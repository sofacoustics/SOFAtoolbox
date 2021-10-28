function [toa_diff,toa,IACC,Obj] = SOFAcalculateITD(Obj,varargin)
% SOFAcalculateITD: Estimate ITD from a binaural signal
%   Usage: itd = SOFAcalculateITD(data,mode,threshlvl,lowpass,butterpoly,upper_cutfreq) 
%
%   Input parameters:
% 
%       data:       SOFA object or IR matrix with dimensions: 
%                     emitter x receiver x time
% 
%       fs:         sampling rate, used only if data provided as matrix
% 
%       mode:       (optional) Select one estimation methods
%                   ('Threshold' (default),'Cen_e2','MaxIACCr', 'MaxIACCe',
%                   'CenIACCr','CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD')
%             
%                   For details concerning estimation methods see:
%                   'http://asa.scitation.org/doi/10.1121/1.4996457'
% 
%       lowpass:    (optional) Bandwidth considered. 'lp' for lowpass (default), 'bb' for broadband
%
%       peak:       (optional) Method to find the max, used in Threshold mode only. 
%                   'hp' for max (default), 'fb' for findpeak
% 
%       threshlvl:  (optional) Set threshold level for 'Threshold' mode in dB.        
%                   Default is -10 dB. 
%
%       butterpoly: (optional) Select the order of the polynom
%                   applied in the butterworth filter. ( 2 =< i =< 10 )
%                   Default is 10. 
% 
%       upper_cutfreq: (optional) Set frequency of lowpass cutoff in Hz.
%                      Default is 3000 Hz. 
% 
%       lower_cutfreq: (optional) Set frequency of highpass cutoff in Hz, 
%                      only used in IRGD mode. Default is 1000 Hz.   
%
%       'debug':     (optional, 0/1): Output debug information about calculations.
% 
% 
%   Output parameters:
% 
%       itd:        interaural time difference in seconds
%       toa:        detected activation onsets for left and right channels
%       IACC:       interaural cross-correlation coefficient
%                   Available on when xcorr is used (modes: 'MaxIACCr', 'MaxIACCe',
%                   'CenIACCr','CenIACCe', 'CenIACC2e')
%       Obj:        Input object with Obj.Data.Delay filled
% 
% 
%   Purpose:
%   Estimates the ITD based on biaural impulse responses.
%   Several different estimaton methods can be chosen.
%   MaxIAACe is recommended.
% 
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%
%   Examples:
%   ---------
% 
%   Obj = SOFAload(fullfile(SOFAdbPath,'baumgartner2017','hrtf b_nh15.sofa'));
%   toa_diff = SOFAcalculateITD(Obj,'MaxIACCe','lp','upper_cutfreq',3000)
%   
%   With these settings the estimator uses the MaxIAAce method and applies
%   a lowpass with a cut off frequency of 3 kHz.
%   
%   The output array is structured as the SOFA Data.IR
%   If you would like to select for example only data on the horizontal
%   plane you could:
%
%   plane_idx = find( Obj.SourcePosition(:,2) == 0 );
%   plane_angle = Obj.SourcePosition(plane_idx,1);
%
%   Url: https://amtoolbox.org/amt-1.0.0/doc/common/itdestimator_code.php


% #Author: Piotr Majdak: Copyright (C) 2009-2021 Piotr Majdak and the AMT team.
% #Author: Michael Mihocic: header documentation updated (20.10.2021)
% 
% This file is part of the SOFA_MO 2.0, 
% basing on the function itdestimator in Auditory Modeling Toolbox (AMT) version 1.0.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% AUTHOR: Laurin Steidle
% data matrix option: Robert Baumgartner
% adaption for SOFA_MO 2.0: Michael Mihocic


% ---------------------- SOFAarghelper -------------------------------
% default parameters (adapted from arg_itdestimator)
definput.keyvals.debug = 0;
definput.flags.mode = {'Threshold','Cen_e2','MaxIACCr', 'MaxIACCe', 'CenIACCr', 'CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD'};
definput.flags.lp = {'lp','bb'};
definput.flags.peak = {'hp','fp'};
definput.flags.toaguess = {'noguess','guesstoa'};
definput.keyvals.threshlvl = -10;
definput.keyvals.butterpoly = 10;
definput.keyvals.upper_cutfreq = 3000;
definput.keyvals.lower_cutfreq = 1000;
definput.keyvals.avgtoa = 45;
definput.keyvals.fs = [];
definput.keyvals.thr = 20;  % dB
definput.keyvals.upsample = 16;  % dB
definput.flags.units = {'time', 'samples'};

[flags,kv]=SOFAarghelper({},definput,varargin);

% ---------------------- renaming input parameter ---------------------

if isstruct(Obj)
  pos = Obj.API.M;
  ear = Obj.API.R;
  Ns  = Obj.API.N;
  IR = Obj.Data.IR;
  fs  = Obj.Data.SamplingRate;
else
  pos = size(Obj,1);
  ear = size(Obj,2);
  Ns  = size(Obj,3);
  IR = Obj;
  if isempty(kv.fs)
    error('RB: No sampling rate (fs) provided.')
  end
  fs  = kv.fs;
end

% ---------------------- initialising variables -----------------------

toa = zeros(pos,ear);
toa_diff = zeros(pos,1);
IACC = zeros(pos,1);
delay = zeros(length(Obj.SourcePosition), 2);

if kv.debug == 1; disp('SOFAcalculateITD:'); end

% ---------------------- Applying low-pass ----------------------------

if flags.do_lp
    if kv.debug == 1
        disp('  Applying Butterworth low pass ')
        disp(strcat('  Polynomial order of Butterworth filter: ',num2str(kv.butterpoly)))
        disp(strcat('  Cut-off frequency is: ',num2str(kv.upper_cutfreq),' Hz'))
    end
    if isoctave; pkg load signal; end
    cut_off_freq_norm = kv.upper_cutfreq/(fs/2);
 
    [lp_a,lp_b] = butter(kv.butterpoly,cut_off_freq_norm);
    f_ir = zeros(pos,ear,Ns);
    for ii=1:pos
        for jj=1:ear  
            sir = squeeze( IR(ii,jj,:) );
            f_sir = filter(lp_a,lp_b,sir);
            f_ir(ii,jj,:) = f_sir;
        end
    end

else
    if kv.debug == 1; disp('  No low pass filter is applied'); end
    f_ir = IR;
end

% ---------------------- estimating itd -------------------------------
% ---------------------------------------------------------------------

% ---------------------- Threshold ------------------------------------
switch(flags.mode)
    case 'Threshold'
        if kv.debug == 1
            disp('  Threshold mode')
            disp(strcat('  Threshold level is: ',num2str(kv.threshlvl),'dB'))
        end

        if flags.do_fp
            for ii=1:pos
                 for jj=1:ear
                    indB = 0.5*mag2db(squeeze(f_ir(ii,jj,:)).^2);
                    [~,B] = findpeaks(indB);
                    th_value = indB(B(1)) + kv.threshlvl;
                    toa(ii,jj) = find(indB>th_value,1);
                end
                toa_diff(ii) = toa(ii,1) - toa(ii,2);     
            end

        else
            for ii=1:pos
                for jj=1:ear
                    indB = 0.5*mag2db(squeeze(f_ir(ii,jj,:)).^2);
                    th_value = max(indB) + kv.threshlvl;
                    toa(ii,jj) = find(indB>th_value,1);
                end
                toa_diff(ii) = toa(ii,1) - toa(ii,2);     
            end
        end


% ---------------------- Cross-Correlation ----------------------------        
    case 'Cen_e2'
        if kv.debug == 1; disp('  Cen-e2 mode'); end
        for ii=1:pos
            for jj = 1:ear
                e_sir_sq = abs(hilbert(squeeze(f_ir(ii,jj,:))).^2);
                toa(ii,jj) = centroid(transpose(1:Ns),e_sir_sq);
            end
            toa_diff(ii) = toa(ii,1) - toa(ii,2);     
        end        


    case 'MaxIACCr'
        if kv.debug == 1; disp('  MaxIACCr mode'); end
        for ii=1:pos                
            cc = xcorr(squeeze(f_ir(ii,1,:)),squeeze(f_ir(ii,2,:)));
            [IACC(ii),idx_lag] = max(abs(cc));
            toa_diff(ii) = idx_lag - Ns;                
        end
        if flags.do_guesstoa
            toa = guesstoa(toa_diff,toa, kv.avgtoa);
        end

    case 'MaxIACCe'
        if kv.debug == 1; disp('  MaxIACCe mode'); end
        for ii=1:pos
            e_sir1 = abs(hilbert(squeeze(f_ir(ii,1,:))));
            e_sir2 = abs(hilbert(squeeze(f_ir(ii,2,:))));
            cc = xcorr(e_sir1,e_sir2);
            [IACC(ii),idx_lag] = max(abs(cc));
            toa_diff(ii) = idx_lag - Ns;
        end
        if flags.do_guesstoa
            toa = guesstoa(toa_diff,toa, kv.avgtoa);
        end


    case 'CenIACCr'
        if kv.debug == 1; disp('  CenIACCr mode'); end
        x = transpose(1:(Ns*2-1));
        for ii=1:pos                
            cc = xcorr(squeeze(f_ir(ii,1,:)),squeeze(f_ir(ii,2,:)));
            pos_cc = abs(cc);
            IACC(ii) = max(pos_cc);
            toa_diff(ii) = centroid(x,pos_cc)-Ns;
        end
        if flags.do_guesstoa
            toa = guesstoa(toa_diff,toa, kv.avgtoa);
        end


    case 'CenIACCe'
        if kv.debug == 1; disp('  CenIACCe mode'); end
        x = transpose(1:(Ns*2-1));
        for ii=1:pos
            e_sir1 = abs(hilbert(squeeze(f_ir(ii,1,:))));
            e_sir2 = abs(hilbert(squeeze(f_ir(ii,2,:))));
            cc = xcorr(e_sir1,e_sir2);
            IACC(ii) = max(abs(cc));
            toa_diff(ii) = centroid(x,abs(cc))-Ns;
        end
        if flags.do_guesstoa
            toa = guesstoa(toa_diff,toa, kv.avgtoa);
        end


    case 'CenIACC2e'
        if kv.debug == 1; disp('  CenIACC2e mode'); end
        x = transpose(1:(Ns*2-1));
        for ii=1:pos              
            e_sir1 = abs(hilbert(squeeze(f_ir(ii,1,:))));
            e_sir2 = abs(hilbert(squeeze(f_ir(ii,2,:))));           
            cc = xcorr(e_sir1,e_sir2).^2;
            IACC(ii) = max(abs(cc));
            toa_diff(ii) = centroid(x,abs(cc))-Ns;    
        end
        if flags.do_guesstoa
            toa = guesstoa(toa_diff,toa, kv.avgtoa);
        end


    case 'PhminXcor'
        if kv.debug == 1; disp('  PhminXcor mode'); end
        ir_min=ARI_MinimalPhase(Obj);
        for ii=1:pos
            for jj=1:ear                    
                cc = xcorr(squeeze(IR(ii,jj,:)),squeeze(ir_min(ii,jj,:)));
                [~,toa(ii,jj)] = max(abs(cc));
            end
            toa_diff(ii) = toa(ii,1) - toa(ii,2);
        end


% ---------------------- Groupdelay -----------------------------------
    case 'IRGD'
        if kv.debug == 1; disp('  IRGD mode'); end
        for ii = 1:pos
            for jj = 1:ear
                f_sir = squeeze( f_ir(ii,jj,:) );
                [gd,w] = grpdelay(transpose(double(f_sir)),1,Ns,fs);
                toa(ii,jj)=mean(gd(find(w>kv.lower_cutfreq): ...
                                    find(w>kv.upper_cutfreq)));
            end
            toa_diff(ii) = toa(ii,1) - toa(ii,2);
        end
end
toa_diff = toa_diff/fs;
toa = toa/fs;

% Calculate Delay, add to Obj
for k = 1:size(IR, 1) 
     % Get onset 
    OnSetL = IR_start(IR(k,1,:), kv.thr); % L
    OnSetR = IR_start(IR(k,2,:), kv.thr); % R        
     % Onset difference
    delay(k,:) = [OnSetL, OnSetR];
end

Obj = SOFAaddVariable(Obj,'Data.Delay','MR',delay);



% -------------------------------------------------------------------------
% ---------------------- Functions ----------------------------------------
% -------------------------------------------------------------------------


% ---------------------- Centroid -----------------------------------------
function idx_cent = centroid(x,y)
idx_cent = sum(x.*y)/sum(y);

% ---------------------- guess toa ----------------------------------------
function toa = guesstoa(toa_diff,toa, avgtoa)
toa(:,1) = toa(:,1) + avgtoa + toa_diff/2;
toa(:,2) = toa(:,2) + avgtoa - toa_diff/2;
    
% ---------------------- Create minimal phase -----------------------------
% as used in ziegelwanger2014
function hMmin=ARI_MinimalPhase(Obj)
hM=Obj.Data.IR;
hMmin=hM;

for jj=1:Obj.API.R
    for ii=1:Obj.API.M
        h=squeeze(hM(ii,jj,:));

        amp1=abs(fft(h));
        amp2=amp1;

        an2u=-imag(hilbert(log(amp1)));
        an2u=an2u(1:floor(length(h)/2)+1);

        an3u=[an2u; -flipud(an2u(2:end+mod(length(h),2)-1))];
        an3=an3u-round(an3u/2/pi)*2*pi;

        amp2=amp2(1:floor(length(h)/2)+1);
        amp3=[amp2; flipud(amp2(2:end+mod(length(h),2)-1))];

        h2=real(ifft(amp3.*exp(1i*an3)));
        hMmin(ii,jj,:)=h2(1:Obj.API.N);
    end
end

% -------------------------- sampleStart ----------------------------------
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