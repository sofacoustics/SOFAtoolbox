function [itd, Obj] = SOFAgetITD(Obj, varargin)
% Calculates Interaural Time Difference (ITD) via the threshold method 
%   Usage:           itd = sofaGetITD(Obj, 'time', 'threshold', 10)
%
%   Input parameters:
%     Obj:           SOFA object with SimpleFreeFieldHRIR convention.
%
%   Output arguments:
%     itd:           Interaural time difference.
%     Obj:           Input object with Obj.Data.Delay filled
%   `SOFAgetITD` accepts the following key-value pairs:
%     'thr',t        Threshold value bellow the IR peak that defines 
%                    the IR onset. (ISO 3382) -- (Default: -10dB).
%     'upsample',v   For output in seconds, the HRIRs are resampled to 
%                    (v) times the original Fs (Default: 16x).

%   `SOFAgetITD` accepts the following flags:
%     'time'         Output is given in with time as unit. (Default).
%     'samples'      Output is given in samples.

% References:
% Brian F. G. Katz and Markus Noisternig. A comparative study of interaural 
%  time delay estimation methods. The Journal of the Acoustical Society of America,
%  page 3530–3540, 2014. doi:10.1121/1.4875714.
% ISO 3382-1:2009 Acoustics — Measurement of room acoustic parameters 

% Author: Davi R. Carvalho, 2021/02/07

%% Parameters
definput.keyvals.thr = 20;  % dB
definput.keyvals.upsample = 16;  % dB
definput.flags.units = {'time', 'samples'};
[flags,kv]=SOFAarghelper({'thr', 'upsample'},definput,varargin);

%% pre-config
itd = zeros(length(Obj.SourcePosition), 1);
delay = zeros(length(Obj.SourcePosition), 2);
IR = shiftdim(Obj.Data.IR, 2);

%% Upsample
if strcmpi(flags.units, 'time')
    fs = Obj.Data.SamplingRate;
    fs_up = fs*kv.upsample;
    Obj_upsample = sofaResample(Obj, fs_up);
    IR = shiftdim(Obj_upsample.Data.IR, 2);
end

%% get ITD
for k = 1:size(IR, 2)
    A = IR(:,k,1); % L
    B = IR(:,k,2); % R    
     % Get onset 
    OnSetL = IR_start(A, kv.thr);
    OnSetR = IR_start(B, kv.thr);        
     % Onset difference
    itd(k) = abs(OnSetL - OnSetR);
    delay(k,:) = [OnSetL, OnSetR];
end

%% Output 
switch flags.units
    case 'time'
        itd = itd./fs_up;
        delay = delay./fs_up;
end
Obj = SOFAaddVariable(Obj,'Data.Delay','MR',delay);
end

%% INTERNAL FUNCTIONS -----------------------------------------------------
function sampleStart = IR_start(IR,threshold)
    % 20210207 - Davi Carvalho, adapted from ita_start_IR.m from https://git.rwth-aachen.de/ita/toolbox/-/blob/master/kernel/DSP/ita_start_IR.m
    threshold = -abs(threshold);
    IR_square = IR.^2; 
    % Max value on IR 
    [pk_val, idx_max] = max(IR_square(:));   
    abs_dat = 10.*log10(IR_square(1:idx_max)) - 10.*log10(pk_val);    
    
    lastBelowThreshold  = find(abs_dat < threshold,1,'last');
    if ~isempty(lastBelowThreshold)
        sampleStart = lastBelowThreshold;
    else
        sampleStart = 1;
    end
    % Check if oscillations exist before the last value below threshold
    % If so, these are part of the RIR and need to be considered.
    idx6dBaboveThreshold = find(abs_dat(1:sampleStart) > threshold + 6);
    if ~isempty(idx6dBaboveThreshold)
         tmp = find(abs_dat(1:idx6dBaboveThreshold(1)) < threshold, 1 ,'last');
         if isempty(tmp) % without this if, the function would generate an error, if the oscillation persists until the first sample
            sampleStart = 1;
        else
            sampleStart = tmp;
         end
    end
end



function Obj = sofaResample(Obj, Fs)
% Change apparent resolution of HRIR in the SOFA object
% Davi R. Carvalho @UFSM - Engenharia Acustica - Junho/2021
%   Input Parameters:
%    Obj:        SimpleFreeFieldHRIR
%    Fs:         Target sample rate
%   Output Parameters:
%     Obj:   SOFA object with the target sampling frequency.
% 
% Matlab R2021a
%% Check if upsampling is necessary first
Fs_sofa = Obj.Data.SamplingRate;
IR = resample_this(Obj.Data.IR, Fs_sofa, Fs);

%% Output
Obj.Data.IR = IR;
% update sampling rate
Obj.Data.SamplingRate = Fs;
Obj = SOFAupdateDimensions(Obj);
end

%--------------------------------------------------------------------------
function IR = resample_this(X, Fs_in, Fs_out)
    [p,q] = rat(Fs_out / Fs_in, 0.0001);
    normFc = .965 / max(p,q);
    order = 256 * max(p,q);
    beta = 12;
    %%% Cria um filtro via Least-square linear-phase FIR filter design
    lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
    lpFilt = lpFilt .* kaiser(order+1,beta)';
    lpFilt = lpFilt / sum(lpFilt);
    lpFilt = p * lpFilt;

    % Initializar matriz
    N_pos = size(X, 1);
    N_ch = size(X, 2);
    N_samples = ceil((Fs_out/Fs_in) * size(X, 3)); % length after resample
    IR=zeros(N_pos, N_ch, N_samples);

    % Actual Resample
    for k = 1:N_pos
        for l = 1:N_ch
            IR(k, l, :) = resample(X(k,l,:), p, q, lpFilt);
        end 
    end
    IR = IR.* q/p; % check scaling

    % make sure signal length is not odd
    if rem(size(IR,3), 2) ~= 0
       IR(:,:,end) = []; 
    end
end