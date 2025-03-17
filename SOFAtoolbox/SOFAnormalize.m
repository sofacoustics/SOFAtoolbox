function [norm_S, param_S] = SOFAnormalize(ori_S, param_S)
%SOFAnormalize - Apply normalization to HRTF measurements (SOFA files).
%   Usage: [norm_S, param_S] = SOFAnormalize(ori_S, param_S)
%
%   SOFAnormalize applis a normalization to SOFA files, consisting six steps:
%     1. low-pass filtering to harmonize the high-frequency support
%     2. temporal alignment of the HRIRs
%     3. temporal window to harmonize the HRIR's duration and envelope
%     4. diffuse-field equalization to remove system-specific colorations
%     5. low-frequency extension to harmonize the low-frequency support
%     6. far-field correction to harmonize across measured distances between the sound-source and listeners
%
%   Input parameters:
%     ori_S: original HRTF SOFA file
%     param_S: (optional) structure specifying normalization parameters. Default values are:
%       param_S.do_gain_norm_b    = 0; % normalize by overall max (all dir, all ears)
%       param_S.do_resamp_b       = 0; % frequency resampling to 'default_Fs_f' (Hz)
%       param_S.default_Fs_f      = 48000;
%       param_S.do_lp_b           = 1; % low-pass filtering at cutoff frequency 'cutfreq_f' (Hz)
%       param_S.cutfreq_f         = 18000;
%       param_S.do_talign_b       = 1; % time align frontal HRIR's onset to reference time mark 'talignSec_f' (sec)
%       param_S.talignSec_f       = 0.001;
%       param_S.threshold_f       = 20; % HRIRs onset is detected as the last point prior to and below -threshold_f (dB) of HRIR peak.
%       param_S.do_win_b          = 1; % HRIRs windowing with window length 'windowLengthSec_f' (sec), starting at HRIRs first onset minus
%                                      % 'safety_f' (sec), including half-Hann fade-in of 'fadeIn_f' (sec) and half-Hann fade-out of 'fadeOut_f' (sec).
%       param_S.windowLengthSec_f = 0.0058;
%       param_S.fadeIn_f          = 0.00025;
%       param_S.safety_f          = param_S.fadeIn_f;
%       param_S.fadeOut_f         = 0.001;
%       param_S.do_resize_b       = 0; % set HRIR length to 'targetLengthSec_f' (sec) by zero-padding or cutting the HRIRs
%       param_S.targetLengthSec_f = param_S.windowLengthSec_f*2;
%       param_S.do_eq_b           = 1; % diffuse-field equalization w. option 'eqfiltFlatten_b' (0 or 1) to flatten the filter at high frequencies
%       param_S.eqfiltFlatten_b   = 1;
%       param_S.do_LFext_b        = 1; % low frequency extension: interp HRTF magnitude at 'lowFreq_f' (Hz) down to SHM gain at DC
%       param_S.lowFreq_f         = 250;
%       param_S.do_dist_b         = 1; % far-field correction
%
%   Output parameters:
%     norm_S: normalized HRTF SOFA file
%     param_S: normalization parameters used

% #Author: Helene Bahu: creation of script HRTF_normalization_v1 (11.03.2025)
% #Author: Michael Mihocic: delayseq substituted (dependance on Phased Array Toolbox), adaptions to SOFA Toolbox, header's format updated (12.03.2025)
% #Author: Helene Bahu: functions updated (13.03.2025)
% #Author: Helene Bahu: code improvements and fixes (17.03.2025)

% #Reference:  H. Bahu, T. Carpentier, M. Noisternig, O. Warusfel, J.-M. Jot, M. Mihocic, and P. Majdak. "Towards an improved consistency between databases of head-related transfer functions." JAES, 2025.

% SOFA Toolbox - function SOFAnormalize
% Copyright (C) Helene Bahu, helenebahu(at)gmail.com; Michael Mihocic, Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.



% Load default normalization parameters
default_param_S = local_get_default_norm_param;
param_names_v = fieldnames( default_param_S );

% Check input arguments
if nargin == 2
    input_param_names_v = fieldnames( param_S );
    missingId_v = find( ~ismember( param_names_v, input_param_names_v ));
    % Fill up missing fields with default parameters
    for ii = 1:length( missingId_v )
        param_S.(param_names_v{missingId_v(ii)}) = default_param_S.(param_names_v{missingId_v(ii)});
    end
elseif nargin == 1
    param_S = default_param_S;
else
    error('invalid number of input arguments')
end

% Check that the input SOFA file meets the SOFA convention for HRIRs
assert( strcmp( ori_S.GLOBAL_SOFAConventions, 'SimpleFreeFieldHRIR' ), 'Only HRIRs are supported' )

% Check that parameters are within reasonable bounds
local_check_parameters( param_S );

% Initialize normalized sofa structure
norm_S = ori_S;

% Prepare data coordinates
if strcmp( norm_S.SourcePosition_Type, 'spherical' )
    sphPos_m = norm_S.SourcePosition;
    cartPos_m = SOFAconvertCoordinates( sphPos_m, 'spherical', 'cartesian' );
end
numDir_n = size( norm_S.SourcePosition, 1 );

% Sampling frequency and number of samples
Fs_f = norm_S.Data.SamplingRate;
numSamples_n = size( norm_S.Data.IR, 3 );


%%% 0 - Normalize by overall max

if param_S.do_gain_norm_b

    max_f = max( max( max( abs( norm_S.Data.IR ) )));
    assert( max_f > 0 )
    norm_S.Data.IR = norm_S.Data.IR./max_f;

end


%%% 1 - Frequency resampling to default_Fs_f

if param_S.do_resamp_b

    if Fs_f ~= param_S.default_Fs_f

        % Initialize
        l_hrir_m = squeeze( norm_S.Data.IR( :, 1, : ));
        r_hrir_m = squeeze( norm_S.Data.IR( :, 2, : ));
        norm_S.Data.IR = [];

        % Resample
        for ii = 1 : numDir_n
            l_hrir_v = l_hrir_m( ii, : );
            r_hrir_v = r_hrir_m( ii, : );
            l_hrir_resmp_m( ii, : ) = resample( l_hrir_v, param_S.default_Fs_f, Fs_f );
            r_hrir_resmp_m( ii, : ) = resample( r_hrir_v, param_S.default_Fs_f, Fs_f );
        end

        % Update variables
        Fs_f = param_S.default_Fs_f;
        numSamples_n = size( norm_S.Data.IR, 3 );

        % Save in norm_S
        norm_S.Data.IR( :, 1, : ) = l_hrir_resmp_m;
        norm_S.Data.IR( :, 2, : ) = r_hrir_resmp_m;
        norm_S.Data.SamplingRate = Fs_f;
        norm_S.API.N = numSamples_n;

    end
end


%%% 2 - Low-pass HRIR

if param_S.do_lp_b

    % Initialize
    l_hrir_m = squeeze( norm_S.Data.IR( :, 1, : ));
    r_hrir_m = squeeze( norm_S.Data.IR( :, 2, : ));
    norm_S.Data.IR = [];

    % To guarantee the same filter regardless of Fs, adapt filter length to Fs
    deltaF = 512/96000;
    %firLength_n = round( deltaF*Fs_f );
    firLengthTemp_n = deltaF*Fs_f;
    firLength_n = 2*floor( firLengthTemp_n/2 )+1;% round to nearest odd number to have a centered filter peak

    % Linear-phase FIR design using firls
    cutoff        = param_S.cutfreq_f;
    frequencies_v = [0 cutoff (cutoff+1000) (Fs_f/2)] ./ (Fs_f/2);
    amplitudes_v  = [1 1 0 0];
    fir_v         = firls( firLength_n-1, frequencies_v, amplitudes_v );

    % Estimate filter lag on slope of phase response below 20kHz
    fft_v           = fft( fir_v );
    phase_v         = unwrap( angle( fft_v ));
    freq_v          = linspace( 0, Fs_f/2, firLength_n/2+1 );
    [ ~, iFreqMax ] = min( abs(freq_v-20000));
    f_v             = freq_v( 1:iFreqMax );
    ph_v            = phase_v( 1:iFreqMax );
    P               = polyfit( f_v, ph_v, 1 );
    filtLag_f       = -P(1)/(2*pi);% sec

    % Zero-padd HRIRs by filter lag
    filtLag_n   = floor( filtLag_f*Fs_f );% samples
    l_hrir_zp_m = [ l_hrir_m zeros( numDir_n, filtLag_n ) ];
    r_hrir_zp_m = [ r_hrir_m zeros( numDir_n, filtLag_n ) ];

    % Filter
    l_hrir_filt_m = filter( fir_v, 1, l_hrir_zp_m, [], 2 );
    r_hrir_filt_m = filter( fir_v, 1, r_hrir_zp_m, [], 2 );

    % % Compensate for time lag and zero-padding
    % l_hrir_lc_m = zeros( size( l_hrir_filt_m ) );
    % r_hrir_lc_m = zeros( size( r_hrir_filt_m ) );
    % for ii = 1 : numDir_n
    %     l_hrir_lc_m(ii,:) = delayseq( l_hrir_filt_m(ii,:).', -filtLag_n );
    %     r_hrir_lc_m(ii,:) = delayseq( r_hrir_filt_m(ii,:).', -filtLag_n );
    % end

    % alternative, getting rid of delayseq:
    % Preallocate output matrices
    l_hrir_lc_m = zeros(size(l_hrir_filt_m));
    r_hrir_lc_m = zeros(size(r_hrir_filt_m));

    % Convert delay to integer
    shift = round(-filtLag_n); % Negative to match delayseq behavior

    % Apply integer delay using slicing (vectorized)
    if shift >= 0
        l_hrir_lc_m(:, shift+1:end) = l_hrir_filt_m(:, 1:end-shift);
        r_hrir_lc_m(:, shift+1:end) = r_hrir_filt_m(:, 1:end-shift);
    else
        shift = abs(shift);
        l_hrir_lc_m(:, 1:end-shift) = l_hrir_filt_m(:, shift+1:end);
        r_hrir_lc_m(:, 1:end-shift) = r_hrir_filt_m(:, shift+1:end);
    end
    % end of alternative


    l_hrir_lp_m = l_hrir_lc_m( :, 1:numSamples_n );
    r_hrir_lp_m = r_hrir_lc_m( :, 1:numSamples_n );

    % Save in norm_S
    norm_S.Data.IR( :, 1, : ) = l_hrir_lp_m;
    norm_S.Data.IR( :, 2, : ) = r_hrir_lp_m;

end


%%% 3 - Frontal HRIR alignment to reference time mark

if param_S.do_talign_b

    % Initialize
    l_hrir_m = squeeze( norm_S.Data.IR( :, 1, : ));
    r_hrir_m = squeeze( norm_S.Data.IR( :, 2, : ));
    front_S = norm_S;
    norm_S.Data.IR = [];

    % Reference time mark in samples
    ref_time_smp_n = param_S.talignSec_f*Fs_f;

    % Make sure HRIR is longer than reference time mark + 1ms
    assert( numSamples_n > ref_time_smp_n+(0.001*Fs_f), 'As designed, time alignment will truncate the HRIRs. Zero-padd HRIRs or modify parameter talignSec_f.' )

    % Index of frontal direction
    iFront_n = find( round( sphPos_m(:,1), 1 ) == 0 & round( sphPos_m(:,2), 1 ) == 0 );
    nifront_n = length( iFront_n );
    if( nifront_n == 0 )
        error( 'Frontal direction was not found. Can''t apply time alignement. Set input parameter do_talign_b to 0 to continue.' )
    elseif( nifront_n > 1 )
        warning( 'Several indices found for frontal direction.' )
        iFront_n = iFront_n(1);
    end

    % Create structure with frontal direction
    front_S.SourcePosition = front_S.SourcePosition( iFront_n, : );
    front_S.Data.IR = front_S.Data.IR( iFront_n, :, : );

    % first onset at frontal direction (btw left and right)
    first_onset_front_n = local_detect_first_onset( front_S, param_S.threshold_f );

    % compute time delay (may be positive or negative)
    tshift_n = round( ref_time_smp_n - first_onset_front_n );
    assert( abs( tshift_n ) < numSamples_n )

    % % Apply delay
    % l_hrir_delayed_m = zeros( numDir_n, numSamples_n );
    % r_hrir_delayed_m = zeros( numDir_n, numSamples_n );
    % for ii = 1 : numDir_n
    %     l_hrir_v = l_hrir_m( ii, : );
    %     r_hrir_v = r_hrir_m( ii, : );
    %     l_hrir_delayed_m( ii, : ) = delayseq( l_hrir_v.', tshift_n );% column vector
    %     r_hrir_delayed_m( ii, : ) = delayseq( r_hrir_v.', tshift_n );
    % end

    % alternative to get rid of delayseq
    % Preallocate output matrices
    l_hrir_delayed_m = zeros(size(l_hrir_m));
    r_hrir_delayed_m = zeros(size(r_hrir_m));

    % Convert delay to integer
    shift = round(tshift_n);

    % Apply integer delay using slicing (vectorized)
    if shift >= 0
        l_hrir_delayed_m(:, shift+1:end) = l_hrir_m(:, 1:end-shift);
        r_hrir_delayed_m(:, shift+1:end) = r_hrir_m(:, 1:end-shift);
    else
        shift = abs(shift);
        l_hrir_delayed_m(:, 1:end-shift) = l_hrir_m(:, shift+1:end);
        r_hrir_delayed_m(:, 1:end-shift) = r_hrir_m(:, shift+1:end);
    end


    % end of alternative


    % Save in norm_S
    norm_S.Data.IR( :, 1, : ) = l_hrir_delayed_m;
    norm_S.Data.IR( :, 2, : ) = r_hrir_delayed_m;

end


%%% 4 - Windowing

if param_S.do_win_b

    % Initialize
    l_hrir_m = squeeze( norm_S.Data.IR( :, 1, : ));
    r_hrir_m = squeeze( norm_S.Data.IR( :, 2, : ));

    % Window length in samples
    windowLength_n = round( param_S.windowLengthSec_f*Fs_f );

    % TOA estimation (using threshold method) to identify first onset
    first_onset_n = local_detect_first_onset( norm_S, param_S.threshold_f );

    % Back safety_n points for fade-in
    safety_n = round( param_S.safety_f*Fs_f );
    first_onset_n = max( first_onset_n-safety_n, 1 );

    % Windowing
    fadeInSmp_n     = round( param_S.fadeIn_f*Fs_f );
    fadeOutSmp_n    = round( param_S.fadeOut_f*Fs_f );
    hannWinFadIn_v  = hann( fadeInSmp_n*2 );
    hannWinFadOut_v = hann( fadeOutSmp_n*2 );
    onesLength_n    = windowLength_n-fadeInSmp_n-fadeOutSmp_n;
    window_v        = [ hannWinFadIn_v(1:fadeInSmp_n).' ones( 1, onesLength_n ) hannWinFadOut_v(fadeOutSmp_n+1:end).' ];

    % Apply window
    norm_S.Data.IR = zeros( numDir_n, 2, numSamples_n );
    if numSamples_n >= first_onset_n + windowLength_n
        norm_S.Data.IR( :, 1, first_onset_n : first_onset_n + windowLength_n-1 ) = l_hrir_m( :, first_onset_n : first_onset_n + windowLength_n-1 ).*window_v;
        norm_S.Data.IR( :, 2, first_onset_n : first_onset_n + windowLength_n-1 ) = r_hrir_m( :, first_onset_n : first_onset_n + windowLength_n-1 ).*window_v;
    else
        norm_S.Data.IR( :, 1, first_onset_n : numSamples_n ) = l_hrir_m( :, first_onset_n : numSamples_n ).*window_v(1:numSamples_n-first_onset_n+1);
        norm_S.Data.IR( :, 2, first_onset_n : numSamples_n ) = r_hrir_m( :, first_onset_n : numSamples_n ).*window_v(1:numSamples_n-first_onset_n+1);
    end

    % Zero-padding if HRIRs are shorter than window length
    if numSamples_n < windowLength_n

        % Initialize
        l_hrir_m = squeeze( norm_S.Data.IR(:,1,:));
        r_hrir_m = squeeze( norm_S.Data.IR(:,2,:));
        norm_S.Data.IR = [];

        % Zero-padd
        norm_S.Data.IR(:,1,:) = [ l_hrir_m zeros( numDir_n, windowLength_n-numSamples_n ) ];
        norm_S.Data.IR(:,2,:) = [ r_hrir_m zeros( numDir_n, windowLength_n-numSamples_n ) ];
        assert( size(norm_S.Data.IR,3) == windowLength_n )
        numSamples_n = windowLength_n;
    end
end

%%% Set HRIRs length to targetLengthSec_f (zero-padd or cut)

if param_S.do_resize_b

    % Target length in samples
    targetLength_n = round( param_S.targetLengthSec_f*Fs_f );

    % Initialize
    l_hrir_m = squeeze( norm_S.Data.IR(:,1,:));
    r_hrir_m = squeeze( norm_S.Data.IR(:,2,:));
    norm_S.Data.IR = [];

    if targetLength_n < numSamples_n % cut
        norm_S.Data.IR(:,1,:) = l_hrir_m(:,1:targetLength_n);
        norm_S.Data.IR(:,2,:) = r_hrir_m(:,1:targetLength_n);
    elseif targetLength_n >= numSamples_n % zero-padd
        norm_S.Data.IR(:,1,:) = [ l_hrir_m zeros( numDir_n, targetLength_n-numSamples_n ) ];
        norm_S.Data.IR(:,2,:) = [ r_hrir_m zeros( numDir_n, targetLength_n-numSamples_n ) ];
    end
    assert( size(norm_S.Data.IR,3) == targetLength_n )
    numSamples_n = targetLength_n;

end


%%% 5 - Equalization

if param_S.do_eq_b

    % Initialize
    l_hrir_m = squeeze( norm_S.Data.IR( :, 1, : ));
    r_hrir_m = squeeze( norm_S.Data.IR( :, 2, : ));
    norm_S.Data.IR = [];

    l_hrirVoronoi_m = l_hrir_m;
    r_hrirVoronoi_m = r_hrir_m;

    l_hrtf_eq_m = zeros( numDir_n, numSamples_n );
    r_hrtf_eq_m = zeros( numDir_n, numSamples_n );

    % FFT of HRIR
    l_hrtf_m = fft( l_hrir_m, numSamples_n, 2 );
    r_hrtf_m = fft( r_hrir_m, numSamples_n, 2 );

    % Remove duplicate spatial points (for Voronoi calculation)
    [ ~, iduplicates_v ] = local_removeDuplicatePoints( cartPos_m );
    cartPosVoronoi_m = cartPos_m;
    if isempty( iduplicates_v ) == 0
        cartPosVoronoi_m( iduplicates_v, : ) = [];
        l_hrirVoronoi_m( iduplicates_v, : )  = [];
        r_hrirVoronoi_m( iduplicates_v, : )  = [];
    end

    % Remove points inside/outside sphere (for Voronoi calculation)
    dist_v = sphPos_m(:,3);
    distUniq_v = unique( round( dist_v, 3 ) );
    default_dist_f = mode( dist_v );
    if length( distUniq_v ) > 1
        warning([ 'The input SOFA file contains ' num2str(length( distUniq_v )) ' distances. Only one distance is considered for Voronoi calculation (in diffuse-field equalization).' ])
        [ ~, iOutsideSph_v ] = local_removePointsOutsideSphere( cartPosVoronoi_m, default_dist_f );
        cartPosVoronoi_m( iOutsideSph_v, : ) = [];
        l_hrirVoronoi_m( iOutsideSph_v, : ) = [];
        r_hrirVoronoi_m( iOutsideSph_v, : ) = [];
    end

    % Voronoi weights
    surfaces_v = local_voronoi( cartPosVoronoi_m );
    weights_norm_v = surfaces_v./sum(surfaces_v);% normalize weights

    % Diffuse-field estimate (on HRIR without duplicate points)
    l_hrtfVoronoi_m = fft( l_hrirVoronoi_m, numSamples_n, 2 );
    r_hrtfVoronoi_m = fft( r_hrirVoronoi_m, numSamples_n, 2 );

    l_squared_mag_m = abs( l_hrtfVoronoi_m ).^2;
    l_ctf_mag_v = sqrt( sum( l_squared_mag_m .* ( weights_norm_v*ones( 1, size( l_squared_mag_m, 2 ))), 1 ));

    r_squared_mag_m = abs( r_hrtfVoronoi_m ).^2;
    r_ctf_mag_v = sqrt( sum( r_squared_mag_m .* ( weights_norm_v*ones( 1, size( r_squared_mag_m, 2 ))), 1 ));


    %%% Flatten diffuse-field filter above 16kHz to mean_mag_f
    if param_S.eqfiltFlatten_b

        % Keep positive frequencies
        is_even_b = ~mod( numSamples_n, 2 );
        if is_even_b
            upper_sample_n = numSamples_n/2+1;
        else
            upper_sample_n = ceil(numSamples_n/2);
        end
        l_ctf_mag_v = l_ctf_mag_v( 1:upper_sample_n );
        r_ctf_mag_v = r_ctf_mag_v( 1:upper_sample_n );
        freq_v = linspace( 0, Fs_f/2, upper_sample_n );


        % frequency limits
        lowFreq_f  = 250;
        highFreq_f = 16000;
        highFreqBound_n = 18000; % frequency at which DF eq filt = mean_mag
        [ ~, iLowfreq_n ]  = min( abs( freq_v - lowFreq_f ));
        [ ~, iHighfreq_n ] = min( abs( freq_v - highFreq_f ));
        [ ~, iHighfreqBound_n ] = min( abs( freq_v - highFreqBound_n ));
        freqInd_v          = iHighfreq_n:iHighfreqBound_n;
        numFreqHigh_n      = length( freqInd_v );
        fade_out_v         = linspace( 1, 0, numFreqHigh_n );
        fade_in_v          = linspace( 0, 1, numFreqHigh_n );

        % mean CTF mag
        l_mean_mag_f = mean( l_ctf_mag_v( iLowfreq_n:iHighfreq_n ));% mean mag in [200 16000]
        r_mean_mag_f = mean( r_ctf_mag_v( iLowfreq_n:iHighfreq_n ));

        % Flattening in high frequencies
        l_ctf_mag_v( freqInd_v ) = l_ctf_mag_v( freqInd_v ) .* fade_out_v +  repmat( l_mean_mag_f, 1, numFreqHigh_n ) .* fade_in_v;
        r_ctf_mag_v( freqInd_v ) = r_ctf_mag_v( freqInd_v ) .* fade_out_v +  repmat( r_mean_mag_f, 1, numFreqHigh_n ) .* fade_in_v;

        l_ctf_mag_v( iHighfreqBound_n : upper_sample_n ) = repmat( l_mean_mag_f, 1, length( iHighfreqBound_n : upper_sample_n ) );
        r_ctf_mag_v( iHighfreqBound_n : upper_sample_n ) = repmat( r_mean_mag_f, 1, length( iHighfreqBound_n : upper_sample_n ) );

        % avoid very low values
        l_ctf_mag_v = max( 5*eps, l_ctf_mag_v );
        r_ctf_mag_v = max( 5*eps, r_ctf_mag_v );

        % Mirror the spectrum
        if is_even_b
            l_ctf_mag_v = [ l_ctf_mag_v l_ctf_mag_v(:,upper_sample_n-1:-1:2) ];
            r_ctf_mag_v = [ r_ctf_mag_v r_ctf_mag_v(:,upper_sample_n-1:-1:2) ];
        else
            l_ctf_mag_v = [ l_ctf_mag_v l_ctf_mag_v(:,upper_sample_n:-1:2) ];
            r_ctf_mag_v = [ r_ctf_mag_v r_ctf_mag_v(:,upper_sample_n:-1:2) ];
        end

    end

    % Minimum-phase CTF
    l_min_ph_ctf_v = imag( hilbert( -log( l_ctf_mag_v ) ) );
    l_ctf_v  = l_ctf_mag_v.* exp( 1i * l_min_ph_ctf_v );

    r_min_ph_ctf_v = imag( hilbert( -log( r_ctf_mag_v ) ) );
    r_ctf_v  = r_ctf_mag_v.* exp( 1i * r_min_ph_ctf_v );

    % Diffuse-field equalization
    for ii = 1 : numDir_n
        l_hrtf_eq_m( ii, : ) = l_hrtf_m( ii, : )./l_ctf_v;
        r_hrtf_eq_m( ii, : ) = r_hrtf_m( ii, : )./r_ctf_v;
    end
    l_hrir_eq_m = real( ifft( l_hrtf_eq_m, numSamples_n, 2 ));
    r_hrir_eq_m = real( ifft( r_hrtf_eq_m, numSamples_n, 2 ));

    % Store in norm_S
    norm_S.Data.IR( :, 1, : ) = l_hrir_eq_m;
    norm_S.Data.IR( :, 2, : ) = r_hrir_eq_m;

end


%%% 6 - Low frequency extension

if param_S.do_LFext_b

    % Initialize
    l_hrir_m = squeeze( norm_S.Data.IR( :, 1, : ));
    r_hrir_m = squeeze( norm_S.Data.IR( :, 2, : ));
    norm_S.Data.IR = [];

    % low freq limit
    is_even_b = ~mod( numSamples_n, 2 );
    if is_even_b
        upper_sample_n = numSamples_n/2+1;
    else
        upper_sample_n = ceil(numSamples_n/2);
    end
    freq_v = linspace( 0, Fs_f/2, upper_sample_n );
    [ ~, iLowfreq_n ] = min( abs( freq_v - param_S.lowFreq_f ));

    % FFT
    l_hrtf_m = fft( l_hrir_m, numSamples_n, 2 );
    r_hrtf_m = fft( r_hrir_m, numSamples_n, 2 );

    % magnitude/phase decompo
    [ l_mag_noflat_m, l_excessph_m ] = local_phase_decompo( l_hrtf_m );
    [ r_mag_noflat_m, r_excessph_m ] = local_phase_decompo( r_hrtf_m );

    % initialization
    l_mag_flat_m = l_mag_noflat_m;
    r_mag_flat_m = r_mag_noflat_m;

    % Get SHM gain at LF for each direction
    freqInd_v = 1:iLowfreq_n;
    for ii = 1 : numDir_n

        [ l_hrir_shm_v, r_hrir_shm_v ] = local_get_shm( sphPos_m(ii,:), numSamples_n, Fs_f );
        l_hrtf_shm_v = fft( l_hrir_shm_v, numSamples_n, 2 );
        r_hrtf_shm_v = fft( r_hrir_shm_v, numSamples_n, 2 );

        l_DCtarg_f = abs( l_hrtf_shm_v(2) );% SHM(DC)=1 by default
        r_DCtarg_f = abs( r_hrtf_shm_v(2) );

        % linear interpolation to SHM gain
        l_mag_flat_m( ii, freqInd_v ) = interp1( [ 1, iLowfreq_n ], [ l_DCtarg_f l_mag_noflat_m( ii, iLowfreq_n ) ], freqInd_v ); % 'makima'
        r_mag_flat_m( ii, freqInd_v ) = interp1( [ 1, iLowfreq_n ], [ r_DCtarg_f r_mag_noflat_m( ii, iLowfreq_n ) ], freqInd_v );

    end

    % phase_recompo & back to time domain
    l_hrtf_flat_m = local_phase_recompo( l_mag_flat_m, l_excessph_m, numSamples_n );
    r_hrtf_flat_m = local_phase_recompo( r_mag_flat_m, r_excessph_m, numSamples_n );
    l_hrir_flat_m = real( ifft( l_hrtf_flat_m, numSamples_n, 2 ));
    r_hrir_flat_m = real( ifft( r_hrtf_flat_m, numSamples_n, 2 ));

    % Store in norm_S
    norm_S.Data.IR( :, 1, : ) = l_hrir_flat_m;
    norm_S.Data.IR( :, 2, : ) = r_hrir_flat_m;

end


%%% 7 - Far-field correction

if param_S.do_dist_b

    % Initialize
    l_hrir_m = squeeze( norm_S.Data.IR( :, 1, : ));
    r_hrir_m = squeeze( norm_S.Data.IR( :, 2, : ));
    norm_S.Data.IR = [];

    % FFT
    l_hrtf_m = fft( l_hrir_m, numSamples_n, 2 );
    r_hrtf_m = fft( r_hrir_m, numSamples_n, 2 );

    % HRTF decomposition into magnitude and excess phase
    [ l_mag_m, l_excess_phases_m ] = local_phase_decompo( l_hrtf_m );
    [ r_mag_m, r_excess_phases_m ] = local_phase_decompo( r_hrtf_m );

    % Initialize
    l_mag_ff_corr_m = l_mag_m;
    r_mag_ff_corr_m = r_mag_m;

    % Apply difference filter at each direction
    for ii = 1 : numDir_n

        az = sphPos_m( ii, 1 );
        el = sphPos_m( ii, 2 );
        dist = sphPos_m( ii, 3 );
        [ l_diff_filt_mag_v, r_diff_filt_mag_v ] = local_design_diff_filter( az, el, dist, numSamples_n, Fs_f );

        l_mag_nf_v = l_mag_m(ii,:);
        r_mag_nf_v = r_mag_m(ii,:);

        l_mag_ff_corr_m( ii, : ) = l_mag_nf_v.*l_diff_filt_mag_v;
        r_mag_ff_corr_m( ii, : ) = r_mag_nf_v.*r_diff_filt_mag_v;
    end


    % mag / phase recomposition
    [ l_hrtf_corr_m ] = local_phase_recompo( l_mag_ff_corr_m, l_excess_phases_m, numSamples_n );
    [ r_hrtf_corr_m ] = local_phase_recompo( r_mag_ff_corr_m, r_excess_phases_m, numSamples_n );

    % back to time domain
    l_hrir_corr_m = real( ifft( l_hrtf_corr_m, numSamples_n, 2  ));
    r_hrir_corr_m = real( ifft( r_hrtf_corr_m, numSamples_n, 2  ));

    % Store in norm_S
    norm_S.Data.IR( :, 1, : ) = l_hrir_corr_m;
    norm_S.Data.IR( :, 2, : ) = r_hrir_corr_m;

end


% Add comment in SOFA GLOBAL_History field
norm_S.GLOBAL_History = [ norm_S.GLOBAL_History ' ↵ Normalized according to Bahu et al. (2025)' ];

% Update number of samples in SOFA API field
norm_S.API.N = numSamples_n;

end


%%%%%%%% LOCAL FUNCTIONS %%%%%%%%

function param_S = local_get_default_norm_param()

% Default Normalization Parameters (in order)

param_S.do_gain_norm_b = 0;% normalize by overall max

param_S.do_resamp_b    = 0; param_S.default_Fs_f = 48000; % frequency resampling

param_S.do_lp_b        = 1; param_S.cutfreq_f = 18000; % low-pass filtering

param_S.do_talign_b    = 1; param_S.talignSec_f = 0.001; param_S.threshold_f = 20; % time align frontal HRIR to reference time mark

param_S.do_win_b       = 1; param_S.windowLengthSec_f = 0.0058; param_S.fadeIn_f = 0.00025; param_S.safety_f = param_S.fadeIn_f; param_S.fadeOut_f = 0.001; param_S.threshold_f = 20; % HRIR windowing

param_S.do_resize_b    = 0; param_S.targetLengthSec_f = param_S.windowLengthSec_f*2; % set HRIR length to target length (zero-padd or cut)

param_S.do_eq_b        = 1; param_S.eqfiltFlatten_b = 1; % diffuse-field equalization

param_S.do_LFext_b     = 1; param_S.lowFreq_f = 250; % low frequency extension

param_S.do_dist_b      = 1;% far-field correction


end

function [] = local_check_parameters( param_S )

% Check that parameters in param_S are within reasonable bounds

% Check boolean parameters "do_.._b"
assert( param_S.do_gain_norm_b == 0 | param_S.do_gain_norm_b == 1, 'Parameter do_gain_norm_b should be a boolean.' )
assert( param_S.do_resamp_b == 0 | param_S.do_resamp_b == 1, 'Parameter do_resamp_b should be a boolean.' )
assert( param_S.do_lp_b == 0 | param_S.do_lp_b == 1, 'Parameter do_lp_b should be a boolean.' )
assert( param_S.do_talign_b == 0 | param_S.do_talign_b == 1, 'Parameter do_talign_b should be a boolean.' )
assert( param_S.do_win_b == 0 | param_S.do_win_b == 1, 'Parameter do_win_b should be a boolean.' )
assert( param_S.do_resize_b == 0 | param_S.do_resize_b == 1, 'Parameter do_resize_b should be a boolean.' )
assert( param_S.do_eq_b == 0 | param_S.do_eq_b == 1, 'Parameter do_eq_b should be a boolean.' )
assert( param_S.do_LFext_b == 0 | param_S.do_LFext_b == 1, 'Parameter do_LFext_b should be a boolean.' )
assert( param_S.do_dist_b == 0 | param_S.do_dist_b == 1, 'Parameter do_dist_b should be a boolean.' )
assert( param_S.eqfiltFlatten_b == 0 | param_S.eqfiltFlatten_b == 1, 'Parameter eqfiltFlatten_b should be a boolean.' )

% Define reasonable bounds for normalization parameters
assert( param_S.default_Fs_f >= 8000 & param_S.default_Fs_f <= 192000 , 'Invalid sampling rate default_Fs_f (Hz).' )
assert( param_S.cutfreq_f >= 100 & param_S.cutfreq_f <= param_S.default_Fs_f/2, 'Invalid cutoff frequency cutfreq_f (Hz).' )
assert( param_S.talignSec_f >= 0 & param_S.talignSec_f <= 0.005, 'Invalid reference time mark talignSec_f (sec).' )
assert( param_S.windowLengthSec_f >= 0.001 & param_S.windowLengthSec_f <= 2, 'Invalid window length windowLengthSec_f (sec).' )
assert( param_S.fadeIn_f >= 0 & param_S.fadeIn_f <= param_S.windowLengthSec_f/2, 'Invalid windowing parameter fadeIn_f (sec).' )
assert( param_S.fadeOut_f >= 0 & param_S.fadeOut_f <= param_S.windowLengthSec_f/2, 'Invalid windowing parameter fadeOut_f (sec).' )
assert( param_S.safety_f >= 0 & param_S.safety_f <= param_S.windowLengthSec_f/2, 'Invalid windowing parameter safety_f (sec).' )
assert( param_S.threshold_f >= 0 & param_S.threshold_f <= 100, 'Invalid onset detection parameter threshold_f (sec).' )
assert( param_S.targetLengthSec_f >= 0.001 & param_S.targetLengthSec_f <= 0.02, 'Invalid resize parameter targetLengthSec_f (sec).' )
assert( param_S.lowFreq_f >= 0 & param_S.lowFreq_f <= 2000, 'Invalid low-frequency extension parameter lowFreq_f (Hz).' )
end

function [ first_onset_n, iDirFirstOnset_n, onset_m ] = local_detect_first_onset( struct_S, threshold_f )

% This function corresponds to local function IR_start in SOFA_calculateITD
% Input: - struct_S: SOFA structure containing field .Data.IR of size [ numDir x 2 x numBins ]
%        - threshold_f: (optional) threshold in dB below peak value for onset detection. Default is 20dB
% Output: - first_onset_n: sample of first onset across directions and ears
%         - iDirFirstOnset_n: index of the direction of the first onset
%         - onset_m: onsets for each direction and ear [ numDir x 2 ]

if nargin == 1
    threshold_f = 20;
end

assert( size( struct_S.Data.IR, 2 ) == 2 )

onset_m = zeros( size( struct_S.Data.IR, 1 ), 2 );

for ii = 1 : size( struct_S.Data.IR, 1) % loop on directions

    for ee = 1 : 2 % loop on ears

        clear IR
        IR = struct_S.Data.IR( ii, ee, : );

        % 20210207 - Davi Carvalho, adapted from ita_start_IR.m from https://git.rwth-aachen.de/ita/toolbox/-/blob/master/kernel/DSP/ita_start_IR.m
        threshold_f = -abs(threshold_f);
        IR_square = IR.^2;
        % Max value on IR
        [pk_val, idx_max] = max(IR_square(:));
        abs_dat = 10.*log10(IR_square(1:idx_max)) - 10.*log10(pk_val);

        lastBelowThreshold  = find(abs_dat < threshold_f,1,'last');
        if ~isempty(lastBelowThreshold)
            sampleStart = lastBelowThreshold;
        else
            sampleStart = 1;
        end
        % Check if oscillations exist before the last value below threshold
        % If so, these are part of the RIR and need to be considered.
        idx6dBaboveThreshold = find(abs_dat(1:sampleStart) > threshold_f + 6);
        if ~isempty(idx6dBaboveThreshold)
            tmp = find(abs_dat(1:idx6dBaboveThreshold(1)) < threshold_f, 1 ,'last');
            if isempty(tmp) % without this if, the function would generate an error, if the oscillation persists until the first sample
                sampleStart = 1;
            else
                sampleStart = tmp;
            end
        end
        onset_m( ii, ee ) = sampleStart;
    end
end

[ minOnsetDir_v, iDirFirstOnsetEars_v ] = min( onset_m, [], 1 );
[ first_onset_n, iEar_n ] = min( minOnsetDir_v );
assert( first_onset_n >= 0 )
iDirFirstOnset_n = iDirFirstOnsetEars_v( iEar_n );
end


function [ xyz_unique_m, iDuplicates_v ] = local_removeDuplicatePoints( xyz_m )

% Remove duplicate points in cartesian coordinates matrix
% Input:  - xyz_m is a matrix of size [ numPoints_n x 3 ]
% Output: - xyz_unique_m: matrix with unique points (duplicates removed)
%         - iDuplicates_v: vector with indices of xyz_m (first dim) that are duplicates

assert( nargin == 1, 'Invalid number of input agruments' )
numPoints_n = size( xyz_m, 1 );
if numPoints_n == 3; warning( 'Only 3 points?? Check matrix dimensions' ); end

% Distance tolerance for considering neighboring points as duplicates (m)
tolerance_f = 10^(-10);
kk = 0;
iDupliPairs_m = [ 0, 0 ];

dist_m = zeros( numPoints_n, numPoints_n );
for ii = 1 : numPoints_n
    for jj = 1 : numPoints_n
        % distance between pair of points ii and jj
        dist_m( ii, jj ) = sqrt( ( xyz_m( ii, 1 )-xyz_m( jj, 1 ))^2 + ( xyz_m( ii, 2 )-xyz_m( jj, 2 ))^2 + ( xyz_m( ii, 3 )-xyz_m( jj, 3 ))^2 );
        % if their distance is under tolerance, they are not the same point, duplicate pair not recorded in iDupliPairs_m
        if( dist_m( ii, jj ) < tolerance_f )&&( ii ~= jj )&&( sum( ismember( iDupliPairs_m, [ jj, ii ], 'rows' ) ) == 0 )
            kk = kk+1;
            iDupliPairs_m( kk, 1 ) = ii;
            iDupliPairs_m( kk, 2 ) = jj;
        end
    end
end

% If there are duplicates, remove them
if kk > 0
    iDuplicates_v = unique( iDupliPairs_m( :, 2 ) );
    unique_ind_v = [ 1 : numPoints_n ];
    unique_ind_v( iDuplicates_v )  = [];
    xyz_unique_m = xyz_m( unique_ind_v, : );

else
    xyz_unique_m = xyz_m;
    iDuplicates_v = [];
end

end

function [ xyz_edit_m, iOutsideSphere_v ] = local_removePointsOutsideSphere( xyz_m, default_dist_f )

% Remove points outside of sphere of radius default_dist_f
% Input: - xyz_m is a matrix of size [ numPoints_n x 3 ]
%        - default_dist_f: radius of sampled sphere
% Output: - xyz_edit_m: matrix with points outside sphere removed
%         - iOutsideSphere_v: vector with indices of xyz_m (1st dim) that are outside the sphere of radius default_dist_f

assert( nargin == 2, 'Invalid number of input agruments' )

distFromCenter_v = sqrt( xyz_m(:,1).^2 + xyz_m(:,2).^2 + xyz_m(:,3).^2 );

iOutsideSphere_v = [];
xyz_edit_m = xyz_m;
uniqDist_v = unique( round( distFromCenter_v, 3 ) );
if length( uniqDist_v ) > 1 % there are several distances
    for dd = 1 : length( uniqDist_v )
        if round( uniqDist_v( dd ), 2 ) ~= round( default_dist_f, 2 ) % if distance dd is different from default_dist_f
            clear iDist_v
            [ iDist_v ] = find( round( distFromCenter_v, 3 ) == uniqDist_v(dd) );
            iOutsideSphere_v = [ iOutsideSphere_v; iDist_v ];
        end
    end
end

xyz_edit_m( iOutsideSphere_v, : ) = [];% remove points
end

function [ mag_m, exc_ph_m ] = local_phase_decompo( hrtf_m )

% Decompose HRTF into magnitude and excess phase
% Input:  - hrtf_m: HRTF matrix of size [ numDir_n x numSamples_n ]
% Output: - mag_m, exc_ph_m: magnitudes and excess phases of size [ numDir_n x numPosFreq ]

assert( nargin == 1 )

numSamples_n = size( hrtf_m, 2 );

mag_allBins_m = abs( hrtf_m );

mag_allBins_m = max( 5*eps, mag_allBins_m );

ph_allBins_m = unwrap( angle( hrtf_m ).').';

min_ph_allBins_m = imag( hilbert( -log( mag_allBins_m ).').');

% Indices of positive frequencies
is_even_b = ~mod( numSamples_n, 2 );
if is_even_b
    upper_sample_n = numSamples_n/2+1;
else
    upper_sample_n = ceil( numSamples_n/2 );
end

exc_ph_m = ph_allBins_m(:,1:upper_sample_n) - min_ph_allBins_m(:,1:upper_sample_n);

mag_m = mag_allBins_m(:,1:upper_sample_n);

end

function [ hrtf_m ] = local_phase_recompo( mag_m, exc_ph_m, numSamples_n )

% Recompose HRTF from magnitude and excess phase
% Input: mag_m, exc_ph_m: magnitudes and excess phases of size [ numDir_n x numPosFreq ]
% Output: hrtf_m: complex HRTF of size [ numDir_n x numSamples_n ]

assert( nargin == 3 )

mag_m = max( 5*eps, mag_m );

upper_sample_n = size( mag_m, 2 );

% Mirror the spectrum
is_even_b = ~mod( numSamples_n, 2 );
if is_even_b
    mag_all_bins_m    = [ mag_m mag_m(:,upper_sample_n-1:-1:2) ];
    exc_ph_all_bins_m = [ exc_ph_m -exc_ph_m(:,upper_sample_n-1:-1:2) ];
else
    mag_all_bins_m    = [ mag_m mag_m(:,upper_sample_n:-1:2) ];
    exc_ph_all_bins_m = [ exc_ph_m -exc_ph_m(:,upper_sample_n:-1:2) ];
end

min_ph_all_bins_m = imag( hilbert( -log( mag_all_bins_m ).').');

phase_all_bins_m = min_ph_all_bins_m + exc_ph_all_bins_m;

phase_unwrap_m = unwrap( phase_all_bins_m.' ).';

hrtf_m = mag_all_bins_m .* exp( 1i * phase_unwrap_m );

end

function [ l_hrir_shm_v, r_hrir_shm_v ] = local_get_shm( sphPos_v, numSamples_n, Fs )

% Get left and right Spherical Head Model HRTFs at one direction
% Input:  - sphPos_v: sphercial coordinates of the point where to compute SHM [ az, el, dist ]
%         - numSamples_n: number of samples of SHM impulse responses
%         - Fs: sampling frequency (Hz)
% Output: - l_hrir_shm_v, r_hrir_shm_v: Left and right SHM impulse responses [ 1 x numSamples_n ]

assert( size( sphPos_v, 1 ) == 1, size( sphPos_v, 2 ) == 3 )

% Default parameters of the SHM
symmetric_ears_v = [90 0];
radius_f = 0.087;
Nsh = 100;

% % verify distance is unique
% dist_v = round( sphPos_v(:,3), 3 );
% uniqDist_f = unique( dist_v );
% if ~isscalar( uniqDist_f )
%     warning([ 'The input SOFA file contains ' num2str(length( uniqDist_f )) ' distances. Only one distance is considered for SHM calculation (in low-frequency extension).' ])
% end
% dist_f = mode( dist_v );
% Compute SHM impulse responses
dist_f = sphPos_v(3);
[ hrir_shm_m ] = local_shm( sphPos_v, symmetric_ears_v, radius_f, dist_f, Nsh, numSamples_n, Fs );

% Prepare output vectors
assert( size( hrir_shm_m, 1 ) == numSamples_n & size( hrir_shm_m, 2 ) == 1 & size( hrir_shm_m, 3 ) == 2 )
l_hrir_shm_v = squeeze( hrir_shm_m( :, 1, 1 )).';
r_hrir_shm_v = squeeze( hrir_shm_m( :, 1, 2 )).';

end

function [ l_diff_filt_mag_v, r_diff_filt_mag_v ] = local_design_diff_filter( az, el, measured_dist_f, numSamples_n, Fs )

% Design difference filters (also called Distance Variation Functions) for HRTF far-field correction
% Input:  - az, el, measured_dist_f: spherical coordinates of measured direction
%         - numSamples_n: number of samples the HRTF measurement to be corrected
%         - Fs: sampling frequency (Hz)
% Output: - l_diff_filt_mag_v, r_diff_filt_mag_v: magnitude ofdifference filters for the left and right HRTF
% NB: az 90 deg. = left
% See Kan et al. JASA, 2009

assert( isscalar( az ), 'provide only one direction' )
assert( length(az) == length(el) & length(az) == length(measured_dist_f) )

% Get far-field SHM HRIR
FF_dist = 100;
%h = local_shm( [ az, el, FF_dist ], [90 0], 0.087, FF_dist, 100, numSamples_n, Fs );% h is [ numSamples_n x 1 x 2 ]
[ l_hrir_ff_v, r_hrir_ff_v ] = local_get_shm([ az, el, FF_dist ], numSamples_n, Fs );

%l_hrir_ff_v = squeeze( h( :, 1, 1 )).';
%r_hrir_ff_v = squeeze( h( :, 1, 2 )).';

l_hrtf_ff_v = fft( l_hrir_ff_v );
r_hrtf_ff_v = fft( r_hrir_ff_v );

l_mag_ff_v = local_phase_decompo( l_hrtf_ff_v );
r_mag_ff_v = local_phase_decompo( r_hrtf_ff_v );

l_mag_ff_v(1) = l_mag_ff_v(2);
r_mag_ff_v(1) = r_mag_ff_v(2);

% Get near-field SHM HRIR
[ l_hrir_nf_v, r_hrir_nf_v ] = local_get_shm([ az, el, measured_dist_f ], numSamples_n, Fs );

l_hrtf_nf_v = fft( l_hrir_nf_v );
r_hrtf_nf_v = fft( r_hrir_nf_v );

l_mag_nf_v = local_phase_decompo( l_hrtf_nf_v );
r_mag_nf_v = local_phase_decompo( r_hrtf_nf_v );

l_mag_nf_v(1) = l_mag_nf_v(2);
r_mag_nf_v(1) = r_mag_nf_v(2);

% Difference filter
l_diff_filt_mag_v = l_mag_ff_v./l_mag_nf_v;
r_diff_filt_mag_v = r_mag_ff_v./r_mag_nf_v;

end

function surfaces_v = local_voronoi( xyz_m )

% Compute the spherical voronoi using matlab convhulln
% Author: Thibaut Carpentier
% Input:  - xyz_m: cartesian coordinates of spherical sampling grid [ numDir_n x 3 ]
% Output: - surfaces_v: surface surrounding each point [ numDir_n x 1 ]

assert( size( xyz_m, 2 ) == 3 )

numPoints = size( xyz_m, 1 );
radius_v = sqrt( xyz_m(:,1).^2 + xyz_m(:,2).^2 + xyz_m(:,3).^2 );
radius = unique( round( radius_v, 3 ));
assert( isscalar( radius ), 'radius not unique')

pointIndices_v  = 1:numPoints;
xyz_m           = xyz_m ./ radius;

% compute convex hull facets
facets_m = convhulln(xyz_m);

% compute facets' normals
numFacets = size(facets_m,1);
vertices_m = zeros(numFacets,3);
for facetIndex = 1:numFacets
    a_v = xyz_m(facets_m(facetIndex,3),:)-xyz_m(facets_m(facetIndex,1),:);
    b_v = xyz_m(facets_m(facetIndex,2),:)-xyz_m(facets_m(facetIndex,1),:);
    normal_v = cross(a_v,b_v);
    vertices_m(facetIndex,:) = normal_v/norm(normal_v);
end

% computes voronoi diagram
voronoi_C = {};
for pointIndex = 1:numPoints
    voronoi_C{pointIndex} = [];
end

for pointIndex = 1:numPoints
    % facets that comprise the point
    [facet_indices_v,positions_v] = find(facets_m==pointIndex);
    found_facets_m = facets_m(facet_indices_v,:);
    for facetIndex = 1:length(found_facets_m)
        found_facets_m(facetIndex,:) = ...
            [found_facets_m(facetIndex,positions_v(facetIndex):end) ...
            found_facets_m(facetIndex,1:positions_v(facetIndex)-1)];
    end
    % sorts these facets
    voronoi_v = zeros(1,size(found_facets_m,1));
    searched_point_n = found_facets_m(1,2);
    for voronoi_n = 1:length(voronoi_v)
        voronoi_v(voronoi_n) = find(found_facets_m(:,2)==searched_point_n);
        searched_point_n = found_facets_m(voronoi_v(voronoi_n),3);
    end
    voronoi_C{pointIndices_v(pointIndex)} = facet_indices_v(voronoi_v)';
end

% computes surfaces
surfaces_v = zeros(numPoints,1);
for pointIndex = 1:numPoints
    if ~isempty(voronoi_C{pointIndex})
        voronoi_v = voronoi_C{pointIndex};
        voronoi_v(end+1) = voronoi_v(1);
        vertex1_v = xyz_m(pointIndex,:);
        for voronoi_n = 1:(length(voronoi_v)-1)
            vertex2_v = vertices_m(voronoi_v(voronoi_n+1),:);
            vertex3_v = vertices_m(voronoi_v(voronoi_n),:);
            % calculate solid angle sustended
            surf_f = 2*atan(dot(vertex1_v,cross(vertex2_v,vertex3_v)) ...
                /(1+dot(vertex2_v,vertex3_v)+dot(vertex3_v,vertex1_v)+dot(vertex1_v,vertex2_v)));

            surfaces_v(pointIndex) = surfaces_v(pointIndex) + surf_f;
        end
    end
end

surfaces_v = surfaces_v .* radius^2;

end

function h = local_shm( sg, ear, a, r_0, Nsh, Nsamples, fs )

% Get analytical solution for a Spherical Head Model
% Input: - sg: desired sampling grid in spherical coordinates [ numDir_n x 3 ]
%        - ear: ear positions [ az, el ]
%        - a: radius of SHM (default = 0.0875)
%        - r_0: distance of sound source (default is r_0 = sg(1,3))
%        - Nsh: sphercak harmonics order (default = 100)
%        - Nsamples: number of samples (full spectrum, default = 1024)
%        - fs: sampling frequency (Hz, default = 44100)
% Output: - h: HRIRs of SHM [ numSamples_n x numDir_n x 2 ]
% This code was COPIED/PASTED from AKsphericalHead.m in AKtools

% speed of sound
c = 343;

% check format of ear vector
if numel(ear) == 2
    ear = [ear 360-ear(1) ear(2)];
end

% calculate great circle distances between the sampling grid and the ears
gcd = [acosd( sind(sg(:,2))*sind(ear(2)) + cosd(sg(:,2))*cosd(ear(2)) .* cosd(sg(:,1)-ear(1)) ); ...
    acosd( sind(sg(:,2))*sind(ear(4)) + cosd(sg(:,2))*cosd(ear(4)) .* cosd(sg(:,1)-ear(3)) )];

% get unique list of great circle distances and radii
[GCD, ~, gcdID] = unique([gcd repmat(sg(:,3), 2, 1)], 'rows');
% gcd = reshape(GCD(gcdID), size(gcd));
r   = GCD(:,2);
GCD = GCD(:,1);

% angle of incidence
theta = GCD/180*pi;

% get list of frequencies to be calculated
f = 0:fs/Nsamples:fs/2;


%%% spherical head model according to:
% [1] R. O. Duda and W. L. Martens "Range dependence of the response of a spherical head model." J. Acoust. Soc. Am., 104(5), 3048-3058 (1998).

% allocate space for output (1st dimension: freq., 2nd dimension: angle)
H = zeros(numel(f), numel(theta));

% get unique list of radii
[rUnique, ~, rID] = unique(r);

% normalized distance - Eq. (5) in [1]
rho_0 = r_0     ./ a;
rho   = rUnique ./ a;

% normalized frequency - Eq. (4) in [1]
mu = (2*pi*f*a) / c;

% Calculate H
for i = 1:length(theta)

    % argument for Legendre polynomial in Eq. (3) in [1]
    x = cos(theta(i));

    % initialize the calculation of the Hankel fraction.
    % Appendix A in [1]
    zr = 1./( 1i* mu * rho( rID(i) ) );
    za = 1./(1i * mu);
    Qr2 = zr;
    Qr1 = zr .* (1-zr);
    Qa2 = za;
    Qa1 = za .* (1-za);

    % initialize legendre Polynom for order m=0 (P2) and m=1 (P1)
    P2 = 1;
    P1 = x;

    % initialize the sum - Eq. (A10) in [1]
    sum = 0;

    % calculate the sum for m=0
    term = zr./(za.*(za-1));
    sum = sum + term;

    % calculate sum for m=1
    if Nsh > 0
        term = (3 * x * zr .* (zr-1)) ./ (za .* (2*za.^2 - 2*za+1));
        sum = sum + term;
    end

    % calculate the sum for 2 <= m <= Nsh
    for m = 2:Nsh

        % recursive calculation of the Legendre polynomial of order m
        % (see doc legendreP)
        P = ((2*m-1) * x * P1 - (m-1) * P2) / m;

        % recursive calculation of the Hankel fraction
        Qr = - (2*m-1) * zr .* Qr1 + Qr2;
        Qa = - (2*m-1) * za .* Qa1 + Qa2;

        % update the sum and recursive terms
        term    = ((2*m+1) * P * Qr) ./ ((m+1) * za .* Qa - Qa1);
        id      = ~isnan(term);         % this might become NaN for high SH orders and low frequencies. However, we usually don't need the high orders for low frequencies anyhow...
        sum(id) = sum(id) + term(id);

        Qr2 = Qr1;
        Qr1 = Qr;
        Qa2 = Qa1;
        Qa1 = Qa;
        P2  = P1;
        P1  = P;
    end

    % calculate the pressure - Eq. (A10) in [1]
    H(:,i) = (rho_0 * exp( 1j*( mu*rho(rID(i)) - mu*rho_0 - mu) ) .* sum) ./ (1i*mu);

end

% [1] uses the Fourier convention with the negative exponent for the
% inverse transform - cf. Eq. (13). Since Matlab uses the opposite
% convention H is conjugated
H = conj(H);

% set 0 Hz bin to 1 (0 dB)
H(1,:) = 1;

% make sure bin at fs/2 is real
if f(end) == fs/2
    H(end,:) = abs(H(end,:));
end

% mirror the spectrum
%Hfull1 = AKsingle2bothSidedSpectrum(H, 1-mod(Nsamples, 2));
upper_sample_n = size( H, 1 );
is_even_b = ~mod( Nsamples, 2 );
if is_even_b
    H = [ H; conj( H( upper_sample_n-1:-1:2, : ))];
else
    H = [ H; conj( H( upper_sample_n:-1:2, : ))];
end

% get the impuse responses
hUnique = ifft(H, 'symmetric');

% add delay to shift the pulses away from the very start
hUnique = circshift(hUnique, [round(1.5e-3*fs) 0]);

% resort to match the desired sampling grid
h = zeros(Nsamples, size(sg,1), 2);
h(:,:,1) = hUnique(:, gcdID(1:size(sg,1) )    );
h(:,:,2) = hUnique(:, gcdID(size(sg,1)+1:end) );

end