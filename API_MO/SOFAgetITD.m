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
definput.flags.units = {'time', 'samples'};
[flags,kv]=SOFAarghelper({'thr'},definput,varargin);


%% Get ITD 
itd = zeros(length(Obj.SourcePosition), 1);
delay = zeros(length(Obj.SourcePosition), 2);
IR = shiftdim(Obj.Data.IR, 2);

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
        itd = itd./Obj.Data.SamplingRate;
        delay = delay./Obj.Data.SamplingRate;
end
Obj = SOFAaddVariable(Obj,'Data.Delay','MR',delay);
end

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