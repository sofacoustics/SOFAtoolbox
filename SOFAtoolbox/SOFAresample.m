function Obj = SOFAresample(Obj, Fs, scale)
%SOFAresample - Modify the sampling rate of an input SOFA object.
%   Usage: Obj = SOFAresample(Obj, Fs)
%          Obj = SOFAresample(Obj, Fs, scale)
%
%   SOFAresample modifies the input SOFA object to a specified sampling rate.
%
%   Input parameters:
%     Obj   : SOFA object (only SimpleFreeFieldHRIR supported)
%     Fs    : Desired output sampling rate
%     scale : Wheter or not to scale the output to have a similar output total energy as the input.
%
%   Output parameters:
%     Obj  : New SOFA object (SimpleFreeFieldHRIR convention)

% #Author: Davi R. Carvalho (09.2022)
% #Author: Michael Mihocic: adapt header for consistency (27.09.2022)

% SOFA Toolbox
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% parse inputs
if nargin < 3
    scale = true;
elseif nargin == 3
    scale = boolean(scale);
end

%% Process
Obj = SOFAconvertConventions(Obj);
Fs_orig = Obj.Data.SamplingRate;
IR = resample_this(Obj.Data.IR, Fs_orig, Fs, scale);

%% Output
Obj.Data.IR = IR;
Obj.Data.SamplingRate = Fs;
Obj = SOFAupdateDimensions(Obj);
end

% --------------------------------------------------------------------------
function IR = resample_this(X, Fs_in, Fs_out, scale)
    [p,q] = rat(Fs_out / Fs_in, 0.0001);
    normFc = .965 / max(p,q);
    order = 256 * max(p,q);
    beta = 12;
    %%% Create a filter via Least-square linear-phase FIR filter design
    if exist('OCTAVE_VERSION','builtin'); pkg load signal; end
    lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
    lpFilt = lpFilt .* kaiser(order+1,beta)';
    lpFilt = lpFilt / sum(lpFilt);
    lpFilt = p * lpFilt;

    % Initialize output matrix
    N_pos = size(X, 1);
    N_ch = size(X, 2);
    N_samples = ceil((Fs_out/Fs_in) * size(X, 3)); % length after resample
    IR=zeros(N_pos, N_ch, N_samples);

    % Actual Resample
    for k = 1:N_pos
        for l = 1:N_ch
            IR(k, l, :) = resample(squeeze(X(k,l,:)), p, q, lpFilt);
        end
    end

    % check scaling
    if scale
        IR = IR.* q/p;
    end

    % make sure signal length is not odd
%     if rem(size(IR,3), 2) ~= 0
%        IR(:,:,end) = [];
%     end
end
