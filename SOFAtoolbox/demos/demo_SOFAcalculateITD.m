%demo_SOFAcalculateITD - Load HRTF and plots ITD.

% #Author: Piotr Majdak
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

%% Define parameters
% Subject index of the file to convert
subject=3;

%% load SOFA file
SOFAfn=fullfile(SOFAdbPath, 'database', 'cipic', ['subject_' sprintf('%03d',subject) '.sofa']);
Obj=SOFAload(SOFAfn, 'nochecks');

%% Calculate Interaural time delay
[itd_time, ~, ~, Obj_time] = SOFAcalculateITD(Obj, 'time', 'thr', 20);
[itd_samples, ~, ~, Obj_samples] = SOFAcalculateITD(Obj, 'samples', 'thr', 20);

%% Plot results
h = figure('Name',mfilename);
subplot(211)
plot((Obj_time.Data.Delay(:,1) - Obj_time.Data.Delay(:,2))*1e6)
xlabel('Position')
ylabel('Time (\mus)')    
title('ITD (time)')
axis tight

subplot(212)
plot((Obj_samples.Data.Delay(:,1) - Obj_samples.Data.Delay(:,2)))
xlabel('Position')
ylabel(['Samples (Fs:' num2str(Obj.Data.SamplingRate), 'Hz)'])       
title('ITD (samples)')
axis tight

%% Polar plot (not working in Octave)
if ~exist('OCTAVE_VERSION','builtin')
    figure('Name',mfilename);
    SOFAplotHRTF(Obj, 'itdhorizontal');
    title('ITD (time, horizontal plane)')
end

