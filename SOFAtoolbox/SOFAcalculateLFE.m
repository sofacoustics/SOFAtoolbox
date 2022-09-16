function Obj_lfe = SOFAcalculateLFE(Obj, fmin, fmax)
%SOFAcalculateLFE - Extend HRTFs towards lower frequencies.
%   Usage: Obj = SOFAcalculateLFE(Obj, fmin, fmax)
% 
%   SOFAcalculateLFE extrapolates each HRTF in Obj below fmax down to fmin 
%   by considering a linear extrapolation in the amplitude spectrum
%   and the minimum-phase version of the phase spectrum.
%
%   Input parameters:
%     Obj  : SOFA object (only SimpleFreeFieldHRIR supported)
%     fmin : Minimal frequency to be calculated (default: 15 Hz)
%     fmax : Frequency to be extended down to fmin (default: 500 Hz)
% 
%   Output parameters:
%     Obj  : New SOFA object

% #Author: Davi R. Carvalho, UFSM - Acoustical Engineering (07.04.2021)
% #Author: Michael Mihocic: adapted for SOFA 2.0 release; header documentation updated (20.10.2021)

% SOFA Toolbox - function SOFAcalculateLFE
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.


%% parse inputs
if nargin < 3
    fmax = 500; % max linear frequency
end
if nargin < 2
    fmin = 15; % minimum freq
end

%% Preprocess
fs = Obj.Data.SamplingRate;
IR = shiftdim(Obj.Data.IR, 2);
N = size(IR, 1);
freq_vec = (0:N/2-1)*fs/N;

% Check if defined fmax data exists in raw hrir 
if length(find(fmax > freq_vec)) == 1
    fmax = freq_vec(3);
end


% Minimum length necessary to contain fmin
N_ext = ceil(fs/fmin); 
if N_ext <= N
    N_ext = N;
    freq_vec_ext = freq_vec;
else
    freq_vec_ext = (0:N_ext/2-1)*fs/N_ext;
end
f500Hz = dsearchn(freq_vec_ext.', fmax); % idx at defined linear part of HRTFs

[~, ~, ~, Obj] = SOFAcalculateITD(Obj, 'samples', 'debug', 1);

%% extrap low frequency
ir_interp = zeros(size(IR, 2), size(IR, 3), N_ext);
for k = 1:size(IR, 2)
    for l = 1:size(IR, 3)
        time = [IR(:,k,l); zeros(ceil(N_ext - N), 1)];
        mag = mag2db(abs(fft(time)));
        mag_interp = mag;
        
        % interp 
        x = [1, f500Hz];
        xq = [1:f500Hz];
        y_mag = [mag(f500Hz); mag(f500Hz)];
        if exist('OCTAVE_VERSION','builtin')
          mag_interp(1:f500Hz) = interp1(x, y_mag, xq);
        else
          mag_interp(1:f500Hz) = interp1(x, y_mag, xq, 'makima');
        end
        
        mag_interp = 10.^(mag_interp./20);
        H = mag_interp(1:round(N_ext/2));
        
        % back to time domain
        ir_interp(k,l,:) = circshift(real(ifft(get_min_phase(abs(H)))), Obj.Data.Delay(k,l));
    end
end


%% OUTPUT
Obj_lfe = Obj;
Obj_lfe.Data.IR = ir_interp;
end


function Hmin = get_min_phase(H, varargin)
% Calculate minimum-phase spectrum from magnitude via the Hilbert Transform
%% Preprocess
H = [H; flip(H)]; % back to double sided spectrum

%% Get minimum_phase
phi_min = imag(hilbert(-(log(abs(H)+eps)))); % eps makes avoids log(0) = -inf
% Complex inverse filtering to consider the phase
Hmin = abs(H).*exp(1i*(phi_min));
end