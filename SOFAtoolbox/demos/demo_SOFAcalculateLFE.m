%demo_SOFAcalculateLFE - Load HRTF and extends low frequency content.

% #Author: Davi Carvalho
% #Author: Michael Mihocic: bugs fixed (10.2021)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 
% 

clear; % clc; close all;

SOFAfile=fullfile(SOFAdbPath,'database','cipic','subject_003.sofa');
Obj=SOFAload(SOFAfile);


%% Low frequency extension 
fmin = 15;
fmax = 500;
Obj_lfe = SOFAcalculateLFE(Obj, fmin, fmax);


%% general properties
fs = Obj.Data.SamplingRate;
% HRIRs
IR     = shiftdim(Obj.Data.IR, 2);
IR_lfe = shiftdim(Obj_lfe.Data.IR, 2);
% Number of samples
N     = size(IR, 1);
N_lfe = size(IR_lfe, 1);
% Time vector 
tx     = 0:1/fs:(N-1)/fs;
tx_ext = 0:1/fs:(N_lfe-1)/fs;
% Frequency vector
freq     = (0:N/2-1)*fs/N;
freq_lfe = (0:N_lfe/2-1)*fs/N_lfe;


%% PLOTS
ch = 1; % ear
pos = 100;   % position index

%%% Plot time 
figure('Name',mfilename);
plot(tx, IR(:,pos,ch)); hold on
plot(tx_ext(1:N), IR_lfe(1:N, pos, ch), '--','linewidth', 1.3); hold off
title('SimpleFreeFieldHRIR, ch: 1, position: 100,  time domain')
legend('original', 'LFE', 'location', 'best')
xlabel('Time (ms)')
axis tight

%%% Plot freq
figure('Name',mfilename);
ori  = mag2db(abs(fft(IR(:,pos,ch), N_lfe)));
lfe = mag2db(abs(fft(IR_lfe(:,pos,ch))));
semilogx(freq_lfe, ori(1:N_lfe/2)); hold on
semilogx(freq_lfe, lfe(1:N_lfe/2), '--','linewidth', 1.3); hold off
title('SimpleFreeFieldHRIR, ch: 1, position: 100, frequency domain')
legend('original', 'LFE', 'location', 'best')
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')
axis tight

% %% Phase
% t = 0:1/1e3:10;
% fo = 0;
% f1 = 500;
% y = chirp(t,fo,t(end),f1,'linear',0,'complex');
% figure('Name',mfilename);
% semilogx(angle(fft(y)))
% title('phase')

%% 
figure('Name',mfilename);
SOFAplotHRTF(Obj,'maghorizontal');
figure('Name',mfilename);
SOFAplotHRTF(Obj_lfe,'maghorizontal');





 