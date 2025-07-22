function [dtf,ctf]=SOFAhrtf2dtf(hrtf,varargin)
%SOFAHRTF2DTF - splits HRTFs into directional transfer functions (DTFs) and common transfer functions (CTFs)
%   Usage:      [dtf,ctf]=SOFAhrtf2dtf(hrtf)
%               [dtf,ctf]=SOFAhrtf2dtf(hrtf,f1,f2)
%
%   Input parameters:
%     hrtf : SOFA object with SimpleFreeFieldHRIR convention
%
%   Output arguments:
%     dtf : SOFA object with the directional transfer functions
%     ctf : SOFA object with the common transfer functions
%
%   `SOFAhrtf2dtf(...)` calculates DTFs using the method from either
%   Majdak et al. (2010; 'log' flag) or Middlebrooks (1999; 'rms' flag).
%   The magnitude spectrum of the CTF is calculated by
%   averaging the (log-)magnitude spectra across all HRTFs for each ear.
%   The phase spectrum of the CTF is the minimum phase
%   corresponding to the magnitude spectrum of the CTF.
%   DTFs result from filtering the HRTF with the inverse complex CTF.
%
%   `SOFAhrtf2dtf` accepts the following key-value pairs:
%
%     'f1',f1     start frequency of the filtering (in Hz; default: 50 Hz)
%     'f2',f2     end frequency of the filtering (in Hz; default: 18 kHz)
%     'atten',a   broadband attenuation in order to avoid clipping (in dB;
%                 default: 20 dB)
%     'weights',w area weights for averaging. Vector of size [M 1]
%                 (M=number of HRIRs), or flase (default: false)
%     'mode','bahu2025' calculate CTF according to Bahu2025, SOFAnormalize.m; Reference: H. Bahu,
%                       T. Carpentier, M. Noisternig, O. Warusfel, J.-M. Jot, M. Mihocic, and P. Majdak.
%                       "Towards an improved consistency between databases of head-related transfer functions." JAES, 2025.
%
%   `SOFAhrtf2dtf` accepts the following flags:
%
%     'log'       evaluate CTF magnitude spectrum by average of log-magnitude
%                 spectra, equivalent to geometric mean of linear filters
%                 (c.f., Majdak et al., 2010; Baumgartner et al., 2014).
%                 This is the default.
%     'rms'       evaluate CTF magnitude spectrum by RMS of linear
%                 magnitude spectra, equivalent to diffuse-field
%                 compensation (c.f., Middlebrooks, 1999; Moller et al.,
%                 1995).
%
%   Restrictions:
%     SOFAhrtf2dtf cannot handle Data.Delay values yet. Issue: https://github.com/sofacoustics/SOFAtoolbox/issues/69
%
% References:
% Baumgartner, R., Majdak, P., & Laback, B. (2014). Modeling sound-source
%  localization in sagittal planes for human listeners. J. Acoust. Soc. Am.
%  136(2), 791-802.
% Majdak, P., Goupell, M. J., & Laback, B. (2010). 3-D localization of
%  virtual sound sources: Effects of visual environment, pointing method,
%  and training. Attention, Perception, & Psychophysics, 72(2), 454-469.
% Middlebrooks, J. C. (1999). Individual differences in external-ear
%  transfer functions reduced by scaling in frequency. J. Acoust. Soc. Am.,
%  106(3), 1480-1492.
% Moller, H., Hammershoi, D., Jensen, C. B., S?rensen, M. F. (1995).
%  Design criteria for headphones. J. Audio Eng. Soc., 43(4), 218-232.

% #Author: Robert Baumgartner (16.01.2014)
% #Author: Fabian Brinkmann: added rms, and weighted averaging (08.09.2016)
% #Author: Piotr Majdak: Octave ifft compatibility (23.02.2018)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: warning if Data.Delay >0; instead of crashing (22.02.2023)
% #Author: Michael Mihocic: method 'bahu2025' implemented (13.03.2025)
% #Author: Michael Mihocic: bug fixed (Octave compatibility) (22.07.2025)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% Check Input
definput.keyvals.f1 = 50;     % Hz
definput.keyvals.f2 = 18000;  % Hz
definput.keyvals.atten = 20;  % dB
definput.keyvals.weights = false;
definput.keyvals.mode = '';   % can also be 'bahu2025'
definput.flags.avg = {'log','rms'};

[flags,kv]=SOFAarghelper({'f1','f2','atten','weights','mode'},definput,varargin);

%% Settings
kv.fs = hrtf.Data.SamplingRate;
hrtfmtx = shiftdim(hrtf.Data.IR,2); % dim 1: time, dim 2: source position, dim 3: receiver/ear
N = size(hrtfmtx,1);
Nfft = N; %2^nextpow2(N);

%% Frequency bounds
df = kv.fs/Nfft;
f = 0:df:kv.fs-df;
idx = f >= kv.f1 & f <= kv.f2;
idx(Nfft/2+2:end) = fliplr(idx(2:Nfft/2));

%% CTF calculation
hrtff=fft(hrtfmtx,Nfft);
if strcmpi(strtrim(kv.mode), 'bahu2025') % method 'Bahu2025'
% if strtrim(lower(convertCharsToStrings(kv.mode))) == lower('bahu2025') % method 'Bahu2025' % not working in Octave

    % Prepare data coordinates
    if strcmp( hrtf.SourcePosition_Type, 'spherical' )
        sphPos_m = hrtf.SourcePosition;
        cartPos_m = SOFAconvertCoordinates( sphPos_m, 'spherical', 'cartesian' );
    end

    cartPosVoronoi_m = cartPos_m;
    surfaces_v = local_voronoi( cartPosVoronoi_m );
    weights_norm_v = surfaces_v./sum(surfaces_v);% normalize weights

    l_hrir_m = squeeze( hrtf.Data.IR( :, 1, : ));
    r_hrir_m = squeeze( hrtf.Data.IR( :, 2, : ));

    l_hrirVoronoi_m = l_hrir_m;
    r_hrirVoronoi_m = r_hrir_m;

    numSamples_n = size( hrtf.Data.IR, 3 );
    l_hrtfVoronoi_m = fft( l_hrirVoronoi_m, numSamples_n, 2 );
    r_hrtfVoronoi_m = fft( r_hrirVoronoi_m, numSamples_n, 2 );

    l_squared_mag_m = abs( l_hrtfVoronoi_m ).^2;
    l_ctf_mag_v = sqrt( sum( l_squared_mag_m .* ( weights_norm_v*ones( 1, size( l_squared_mag_m, 2 ))), 1 ));

    r_squared_mag_m = abs( r_hrtfVoronoi_m ).^2;
    r_ctf_mag_v = sqrt( sum( r_squared_mag_m .* ( weights_norm_v*ones( 1, size( r_squared_mag_m, 2 ))), 1 ));

    ctfmtx(:,:,1)=shiftdim(l_ctf_mag_v,1);
    ctfmtx(:,:,2)=shiftdim(r_ctf_mag_v,1);

    % get complex spectrum
    if flags.do_rms
        ctff = ctfmtx .*exp(1j*angle(ctfmtx));
    else
        ctff = exp(ctfmtx) .*exp(1j*angle(ctfmtx));
    end
else

    % get magnitude response
    if flags.do_rms
        if any(kv.weights)
            kv.weights = squeeze(kv.weights) / sum(kv.weights);
            ctffavg=sqrt(sum( abs(hrtff).^2 .* repmat(kv.weights', [hrtf.API.N 1 hrtf.API.R]), 2 ));
        else
            ctffavg=sqrt(mean(abs(hrtff).^2,2));
        end
    else % flags.do_log
        if any(kv.weights)
            kv.weights = squeeze(kv.weights) / sum(kv.weights);
            ctffavg= sum(log(abs(hrtff)+eps) .* repmat(kv.weights', [hrtf.API.N 1 hrtf.API.R]), 2);
        else
            ctffavg=mean(log(abs(hrtff)+eps),2);
        end
    end

    % Force minimum phase
    ctfflog=mean(log(abs(hrtff)+eps),2);
    ctfcep = ifft(ctfflog,Nfft);
    ctfcep(Nfft/2+2:Nfft,:,:) = 0;    % flip acausal part to causal part or simply multiply
    ctfcep(2:Nfft/2,:,:) = 2*ctfcep(2:Nfft/2,:,:);    % causal part by 2 (due to symmetry)
    ctfflog = fft(ctfcep,Nfft);
    ctfp = exp(ctfflog);

    % get complex spectrum
    if flags.do_rms
        ctff = ctffavg .*exp(1j*angle(ctfp));
    else
        ctff = exp(ctffavg) .*exp(1j*angle(ctfp));
    end
    % get IR
    ctfmtx = real(ifft(ctff, Nfft));
end

%% DTF calculation
dtff = hrtff;
dtff(idx,:,:) = hrtff(idx,:,:)./repmat(ctff(idx,:,:),[1 size(hrtff,2) 1]);
dtfmtx = ifft(dtff,Nfft);

%% Attenuate to avoid clipping
ctfmtx = ctfmtx / 10^(kv.atten/20);
dtfmtx = dtfmtx / 10^(kv.atten/20);

%% Output Objects
dtf = hrtf;
if size(dtfmtx,3)==1 % handle objects with R=1
    dtfmtx(1,1,2)=0;
    dtfmtx = shiftdim(dtfmtx,1);
    dtf.Data.IR = dtfmtx(:,1,:);
else
    dtf.Data.IR = shiftdim(dtfmtx,1);
end
if ~isfield(dtf, 'GLOBAL_Comment'), dtf.GLOBAL_Comment=''; end
dtf.GLOBAL_Comment = [dtf.GLOBAL_Comment '. Directional transfer functions (DTFs) were generated by removing from the HRTFs the direction-independent log-amplitude spectrum for each ear.'];
if max(max(abs(dtf.Data.Delay))) > 0
    warning('SOFAhrtf2dtf does not support ''Data.Delay'' values larger than zero (yet).')
end

ctf = hrtf;
if size(ctfmtx,2)==1 % handle objects with R=1
    ctf.Data.IR = shiftdim(shiftdata([ctfmtx ctfmtx],3),2);
    ctf.Data.IR = ctf.Data.IR(1,:,:);
else
    ctf.Data.IR = shiftdim(ctfmtx,1);
end
ctf.Data.Delay=zeros (1, ctf.API.R); % Data.Delay not supported yet, needs to be fixed

ctf.GLOBAL_Comment = [dtf.GLOBAL_Comment '. Common transfer functions (CTFs) were extracted from the HRTFs in terms of averaging the log-amplitude spectrum for each ear. The phase of a CTF is determined as minimum phase.'];
ctf.SourcePosition = [0 0 0];
% remove non-standard variables
fc=fieldnames(ctf);
fo=fieldnames(SOFAgetConventions(ctf.GLOBAL_SOFAConventions));
f=setdiff(fc,fo); % get non-standard fields
f(strncmpi(f,'GLOBAL_',7)) = []; % variables only
for ii = 1:numel(f)
    if isfield(ctf, f{ii})
        ctf=rmfield(ctf,f{ii}); % remove
    end
    if isfield(ctf.API.Dimensions,f{ii})
        ctf.API.Dimensions=rmfield(ctf.API.Dimensions,f{ii});
    end
end
clear ii
ctf = SOFAupdateDimensions(ctf);


function f=myifftreal(c,N) % thanks goto the LTFAT <http://ltfat.sf.net>
if rem(N,2)==0
    f=[c; flipud(conj(c(2:end-1,:)))];
else
    f=[c; flipud(conj(c(2:end,:)))];
end;
f=real(ifft(f,N,1));


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
