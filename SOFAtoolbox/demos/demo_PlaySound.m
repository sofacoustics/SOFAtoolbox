%demo_PlaySound - Test audio playback after convolving sound file.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: bugs fixed (10.2021)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% put your information here:
hrtf = SOFAload('path/to_your/HRTF.sofa'); % Load your impulse response into a struct
soundInput = audioread('path/to_your/fancy_audio_file.wav'); % Load your sound file

%% demo script
% Start SOFA
SOFAstart;
% Display some information about the impulse response
SOFAinfo(hrtf);
% Plot a figure with the measurement setup
SOFAplotGeometry(hrtf);
% Have a look at the size of the data
disp(['size [MxRxN]: ' num2str(size(hrtf.Data.IR))])
% Calculate the source position from a listener point of view
apparentSourceVector = SOFAcalculateAPV(hrtf);
% Listen to the HRTF with azimuth of -90°
apparentSourceVector(91, 1)
SOFAplotGeometry(hrtf, 91);
soundOutput = [conv(squeeze(hrtf.Data.IR(91, 1, :)), soundInput) ...
               conv(squeeze(hrtf.Data.IR(91, 2, :)), soundInput)];
sound(soundOutput, hrtf.Data.SamplingRate);