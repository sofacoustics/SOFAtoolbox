function Obj_lfe = SOFAlfe(Obj, fmin, fmax)
% Low frequency extention
%   Description:
%                    This function extrapolates low frequency content of 
%                    SOFA SimpleFreeFieldHRIR objects by considering a
%                    linear behavior to the low frequency content

%   Usage:           Obj = SOFAalfe(Obj, fmin, fmax)
%
%   Input parameters:
%     Obj:           SOFA object with SimpleFreeFieldHRIR convention.
%     fmin:          Minimal frequency to be calculated (default: 15Hz)
%     fmax:          Reference frequency value to be extended until fmin
%                                                       (default: 500Hz)
% 
%   Output arguments:
%     Obj:           SOFA object with 

%   `SOFAgetITD` accepts the following flags:
%     'time'         Output is given in with time as unit. (Default).
%     'samples'      Output is given in samples.


% Author: Davi R. Carvalho, 2021/04/07 @UFSM - Acoustical Engineering


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


%% extrap low frequency
ir_interp = zeros(size(IR, 2), size(IR, 3), N_ext);
for k = 1:size(IR, 2)
    for l = 1:size(IR, 3)
        time = [IR(:,k,l); zeros(ceil(N_ext - N), 1)];
        mag = fft(time);
        mag_interp = mag;
        
        % interp 
        x = [freq_vec_ext(2),    freq_vec_ext(f500Hz:f500Hz+1)];
        xq = freq_vec_ext(2:f500Hz);
        y_mag = [mag(f500Hz); mag(f500Hz:f500Hz+1)];
        mag_interp(2:f500Hz) = interp1(x, y_mag, xq, 'makima');

        % back to time domain
        ir_interp(k,l,:) = real(ifft(mag_interp, N_ext, 'symmetric'));
    end
end


%% OUTPUT
Obj_lfe = Obj;
% "Normalize"
ir_interp = ir_interp./max(abs(ir_interp(:))) .* max(abs(IR(:)));
Obj_lfe.Data.IR = ir_interp;
end