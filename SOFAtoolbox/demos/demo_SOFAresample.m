%demo_SOFAresample - Demo script to resample SOFA object, using SOFAresample function.

% #Author: Davi Carvalho (09.2022)
% #Author: Michael Mihocic: adapt header for consistency, minor code modifications (27.09.2022)

% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 
% 

% clear;  clc; close all;

SOFAfile=fullfile(SOFAdbPath,'database','cipic','subject_003.sofa');
Obj=SOFAload(SOFAfile);

%% Resample
Fs_target = 96000;
Obj_out = SOFAresample(Obj, Fs_target);

%% general properties
fs = Obj.Data.SamplingRate;
% HRIRs
IR     = shiftdim(Obj.Data.IR, 2);
IR_out = shiftdim(Obj_out.Data.IR, 2);
% Number of samples
N     = size(IR, 1);
N_out = size(IR_out, 1);
% Time vector 
tx     = 0:1/fs:(N-1)/fs;
tx_out = 0:1/Fs_target:(N_out-1)/Fs_target;
% Frequency vector
freq     = (0:N_out/2-1)*fs/N_out;
freq_out = (0:N_out/2-1)*Fs_target/N_out;

%% PLOTS
ch = 1; % ear
pos = 1;   % position index

%%% Plot time 
figure('Name',mfilename);
plot(tx, IR(:,pos,ch)); hold on
plot(tx_out, IR_out(:, pos, ch), '--','linewidth', 1.3); hold off
xlim([0, min(N_out, N)])
legend('original', 'resampled', 'location', 'best')
xlabel('Time (ms)')
axis tight
title('SimpleFreeFieldHRIR: original and resampled')

%%% Plot freq
figure('Name',mfilename);
ori = mag2db(abs(fft(IR(:,pos,ch), N_out)));
out = mag2db(abs(fft(IR_out(:,pos,ch))));
semilogx(freq, ori(1:N_out/2)); hold on
semilogx(freq_out, out(1:N_out/2), '--','linewidth', 1.3); hold off
legend('original', 'resampled', 'location', 'best')
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')
axis tight
title('SimpleFreeFieldHRIR: original and resampled')