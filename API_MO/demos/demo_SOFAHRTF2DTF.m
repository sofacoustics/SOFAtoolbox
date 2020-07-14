% SOFA API - demo script
% Copyright (C) 2012-2014 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 
% 
% load HRTF and plots CTF and average DTF
% 

%% Define parameters
% Subject index of the file to convert
subject=3;
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% load SOFA file
SOFAfn=fullfile(SOFAdbPath, 'database', 'cipic', ['subject_' sprintf('%03d',subject) '.sofa']);
X=SOFAload(SOFAfn);
[D,C]=SOFAhrtf2dtf(X);

% close all;
f = figure;
set(f, 'Position', [50, 400, 1100, 500]); 

subplot(1,2,1);
data=(20*log10(abs(fft(squeeze(C.Data.IR(1,1,:))))));
stepsize=C.Data.SamplingRate/length(data)*2;
plot(1:stepsize:C.Data.SamplingRate,data(1:length(data)/2));
% ax=gca; set(ax,'XScale','log');
hold on;    grid on;

subplot(1,2,2);
% plot(20*log10(abs(fft(squeeze(C.Data.IR(1,1,1:length(C.Data.IR)/2))))));
data=(20*log10(abs(fft(squeeze(C.Data.IR(1,2,:))))));
stepsize=C.Data.SamplingRate/length(data)*2;
plot(1:stepsize:C.Data.SamplingRate,data(1:length(data)/2));
% ax=gca; set(ax,'XScale','log');
hold on;    grid on;

[D,CC]=SOFAhrtf2dtf(D);

subplot(1,2,1);
data=(20*log10(abs(fft(squeeze(CC.Data.IR(1,1,:))))));
stepsize=CC.Data.SamplingRate/length(data)*2;
plot(1:stepsize:CC.Data.SamplingRate,data(1:length(data)/2),'r');
title('left');
xlabel('f (in Hz)'); ylabel('dB');
legend('CTF','Avg DTF','Location','Best')

subplot(1,2,2);
% plot(20*log10(abs(fft(squeeze(CC.Data.IR(1,1,1:length(C.Data.IR)/2))))),'r');
data=(20*log10(abs(fft(squeeze(CC.Data.IR(1,2,:))))));
stepsize=CC.Data.SamplingRate/length(data)*2;
plot(1:stepsize:CC.Data.SamplingRate,data(1:length(data)/2),'r');
title('right');
xlabel('f (in Hz)'); ylabel('dB');
legend('CTF','Avg DTF','Location','Best');
