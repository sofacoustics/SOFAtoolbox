function Obj = SOFAresample(Obj, Fs, scale)
%SOFAresample - Modify the sampling rate of an input SOFA object.
%   Usage: Obj = SOFAresample(Obj, Fs)
%          Obj = SOFAresample(Obj, Fs, scale)
%
%   SOFAresample resamples the data in Obj to a specified sampling rate Fs. 
%   If Obj is not in SingleFreeFieldHRIR, it will be converted with the help 
%   of SOFAconvertConventions. 
%
%   Per default, SOFAresample scales the output data to match its energy 
%   to the input data. This corresponds to sampling as a process of energy 
%   preservation in the continuous function and creates amplitude spectra 
%   of the input and output at the same level. 
%   
%   By calling SOFAresample with scale 'noscale', the output data are not
%   scaled such that the output amplitudes match the input amplitudes. This
%   corresponds to sampling as a process of amplitude snapshots of the 
%   continuous function and creates impulse responses of the input and 
%   output at the same level. 
%
%   Input parameters:
%     Obj   : SOFA object 
%     Fs    : Output sampling rate (Hz)
%     scale : 'scale' for scaling the output (default) or 'noscale' otherwise. 
%
%   Output parameters:
%     Obj  : New SOFA object (SimpleFreeFieldHRIR convention)

% #Author: Davi R. Carvalho (09.2022)
% #Author: Michael Mihocic: adapt header for consistency (27.09.2022)
% #Author: Piotr Majdak: clean up and parameter adaptations (4.10.2022)

% SOFA Toolbox
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% parse inputs
if nargin < 3
    do_scale = true;
elseif nargin == 3
    if strcmpi(scale,'scale')
      do_scale=true; 
    else
      do_scale=false; 
    end
end
%% load packages if required and prepare conventions
if exist('OCTAVE_VERSION','builtin'); pkg load signal; end
Obj = SOFAconvertConventions(Obj);
%% Resample
  % prepare data
Fs_in = Obj.Data.SamplingRate;
X = Obj.Data.IR; 
M = Obj.API.M; 
N = Obj.API.N; 
R = Obj.API.R;
[p,q] = rat(Fs / Fs_in, 0.0001); % calculate the rational fraction between the input and output sampling rate
  % pre-calculate the resample filter coefficients for speeding up the resample function
normFc = .965 / max(p,q);
order = 256 * max(p,q);
beta = 12;
lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
lpFilt = lpFilt .* kaiser(order+1,beta)';
lpFilt = lpFilt / sum(lpFilt);
lpFilt = p * lpFilt;
  % Initialize the output matrix
M = size(X, 1);
R = size(X, 2);
N = ceil((Fs/Fs_in) * size(X, 3)); % length after resample
IR=zeros(M, R, N);
  % Resample
for m = 1:M
    for r = 1:R
        IR(m,r,:) = resample(squeeze(X(m,r,:)), p, q, lpFilt);
    end
end
% do scaling
if do_scale, IR = IR.* q/p; end

%% Compile the output structure
Obj.Data.IR = IR;
Obj.Data.SamplingRate = Fs;
Obj = SOFAupdateDimensions(Obj);
end
