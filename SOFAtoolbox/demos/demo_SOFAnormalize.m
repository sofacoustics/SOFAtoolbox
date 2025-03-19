%demo_SOFAnormalize - Demo script applying normalization to HRTF SOFA files, using SOFAnormalize.

% #Author: Helene Bahu: creation of test script (11.03.2025)
% #Author: Michael Mihocic: adaptions to SOFA Toolbox, header's format updated (12.03.2025)
% #Author: Helene Bahu: demo script adapted to plot meaningful data (as per Bahu2025 paper) (17.03.2025)
% #Author: Michael Mihocic: minor adaptions to fit SOFA Toolbox (17.03.2025)
% #Author: Michael Mihocic: SS2 SOFA files added to the loop (19.03.2025)
% 
% SOFA Toolbox - demo script
% Copyright (C) Helene Bahu, helenebahu(at)gmail.com; Michael Mihocic, Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Clean up (optional)
% close all;
% clear;
% clear all % 'clear all' can break SOFA functionality
% clc

%% Adapt system parameters here (optional)
% addpath to SOFAstart
% addpath(genpath('D:\Projects\SOFA\Github\SOFAtoolbox-development\SOFAtoolbox'))
% addpath(genpath('/Users/bahu/Documents/MATLAB/SOFA Toolbox/SOFAtoolbox/'))
SofaFiles = {...
    'clubfritz/ClubFritz1.sofa', ...
    'clubfritz/ClubFritz2.sofa',...
    'clubfritz/ClubFritz3.sofa',...
    'clubfritz/ClubFritz4.sofa',...
    'clubfritz/ClubFritz5.sofa',...
    'clubfritz/ClubFritz6.sofa',...
    'clubfritz/ClubFritz7.sofa',...
    'clubfritz/ClubFritz8.sofa',...
    'clubfritz/ClubFritz9.sofa',...
    'clubfritz/ClubFritz10.sofa',...
    'clubfritz/ClubFritz11.sofa',...
    'clubfritz/ClubFritz12.sofa',...
    'bili (hrtf)/IRC_1130_R_HRIR_96000.sofa',...
    'thk/HRIR_L2702_NF050.sofa',...
    'thk/HRIR_L2702_NF150.sofa',...
    'sadie/D1_48K_24bit_256tap_FIR_SOFA.sofa' ...
    'ss2%20(mannequins)/KU100051023_1_processed.sofa' ...
    'ss2%20(mannequins)/KU100051023_2_processed.sofa' ...
    'ss2%20(mannequins)/KU100051023_3_processed.sofa' ...
    'ss2%20(mannequins)/KU100051023_4_processed.sofa' ...
    }; 


%% Process data
% SOFAstart;
% Define a list of file names

% figure
legend_s = [];    
ylim_v = [-35 15];
xlim_v = [90 20000];
% Loop through all SOFA files
for i = 1:length(SofaFiles)
    SofaFile = SofaFiles{i};
    % SOFAload downloads files from SOFA Conventions repository.
    % disp(['Loading file ' i '/' length(SofaFiles) ": " SofaFile]);
    Obj = SOFAload(['db://database/' SofaFile]);
    % Plot original magnitude at frontal direction
    [ l_mag_ori_v, r_mag_ori_v, freq_ori_v ] = local_get_frontal_mag( Obj );
    subplot( 2, 1, 1 )
    semilogx( freq_ori_v, 20*log10( l_mag_ori_v ) )
    hold on
    title('Original')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    grid on
    ylim(ylim_v)
    xlim(xlim_v)
    % plot original files
    % figure('Name',SofaFile);
    % subplot(2,2,1);
    % SOFAplotHRTF(Obj,'ETCHorizontal',1);
    % title('ETC Horizontal')
    % % figure('Name',SofaFile);
    % subplot(2,2,2);
    % SOFAplotHRTF(Obj,'MagMedian',2);
    % title('Mag Median')

    % To modify the normalization parameters, use the following, and add param_S as a 2nd input parameter in HRTF_normalization_v1
    % param_S.do_gain_norm_b = 1;
    % param_S.do_resamp_b    = 0;
    % param_S.do_lp_b        = 0;
    % param_S.do_talign_b    = 0;
    % param_S.do_win_b       = 0;
    % param_S.do_zp_b        = 0;
    % param_S.do_eq_b        = 0;
    % param_S.do_LFext_b     = 0;
    % param_S.do_dist_b      = 0;

    % Apply normalization
    % disp(['Normalizing file ' i '/' length(SofaFiles) ": " SofaFile]);
    param_S.do_resize_b = 1;
    Objnorm_S = SOFAnormalize(Obj,param_S);
    % Plot normalized magnitude at frontal direction
    [ l_mag_norm_v, r_mag_norm_v, freq_norm_v ] = local_get_frontal_mag( Objnorm_S );
    subplot( 2, 1, 2 )
    semilogx( freq_norm_v, 20*log10( l_mag_norm_v ) )
    hold on
    title('Normalized')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    grid on
    ylim(ylim_v)
    xlim(xlim_v)
    legend_s = strvcat( legend_s, num2str(i) );
    % % plot normalized data
    % % figure('Name',[SofaFile ' (normalized)']);
    % subplot(2,2,3);
    % SOFAplotHRTF(Objnorm_S,'ETCHorizontal',1);
    % title('ETC Horizontal, normalized')
    % % figure('Name',[SofaFile ' (normalized)']);
    % subplot(2,2,4);
    % SOFAplotHRTF(Objnorm_S,'MagMedian',2);
    %  title('Mag Median, normalized')
end
legend(legend_s)
% set( gcf, 'Position', [1   917   560   420 ])

function [ l_mag_v, r_mag_v, freq_v ] = local_get_frontal_mag( Objnorm_S )

    % Get left and right magnitudes at frontal direction and frequency
    % vector from SOFA object
    Fs_f = Objnorm_S.Data.SamplingRate;

    % Find index of frontal magnitude
    ind_n = find( round( Objnorm_S.SourcePosition( :, 1 ), 1 ) == 0 & round( Objnorm_S.SourcePosition( :, 2 ), 1 ) == 0 );
    if length(ind_n)>1; ind_n=ind_n(1); end
    if isempty(ind_n); mag_v=[]; freq_v=[]; warning('no frontal direction'); end

    % Indices of positive frequencies
    numSamples_n = size(Objnorm_S.Data.IR,3);
    is_even_b = ~mod( numSamples_n, 2 );
    if is_even_b
        upper_sample_n = numSamples_n/2+1;
    else 
        upper_sample_n = ceil( numSamples_n/2 );
    end

    % Get magntiudes and frequency vector
    l_hrir_v = squeeze( Objnorm_S.Data.IR( ind_n, 1, : )).';
    r_hrir_v = squeeze( Objnorm_S.Data.IR( ind_n, 2, : )).';
    l_hrtf_v = fft( l_hrir_v, numSamples_n );
    r_hrtf_v = fft( r_hrir_v, numSamples_n );
    l_mag_v  = abs( l_hrtf_v( 1:upper_sample_n ) );
    r_mag_v  = abs( r_hrtf_v( 1:upper_sample_n ) );
    freq_v   = linspace( 0, Fs_f/2, numSamples_n/2+1 );

end
