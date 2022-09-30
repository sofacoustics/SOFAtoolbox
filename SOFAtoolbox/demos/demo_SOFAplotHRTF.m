%demo_SOFAplotHRTF - Script demonstrating the usage of SOFAplotHRTF.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: bugs fixed (10.2021)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% 
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% load a SOFA file in SimpleFreeFieldHRIR
SOFAfile=fullfile(SOFAdbPath,'database','ari','dtf_nh2.sofa');
disp(['Loading: ' SOFAfile]);
Obj=SOFAload(SOFAfile);

% plot ETC horizontal plane
figure('Name',mfilename);
SOFAplotHRTF(Obj,'ETCHorizontal',1);
% plot magnitude spectrum in the median plane, channel 2
figure('Name',mfilename);
SOFAplotHRTF(Obj,'MagMedian',2);
% plot non-normalized magnitude spectrum in the median plane, channel 1
figure('Name',mfilename);
SOFAplotHRTF(Obj,'MagMedian','nonormalization');

%% load a GeneralTF SOFA file
SOFAfile=fullfile(SOFAdbPath,'database','ari (sim)','hrtf_nh5_ref.sofa');
Obj=SOFAload(SOFAfile);
% plot magnitude spectrum in the median plane, channel 1
figure('Name',mfilename);;
SOFAplotHRTF(Obj,'MagMedian',1,'conversion2ir');
figure('Name',mfilename);
SOFAplotHRTF(Obj,'MagMedian',1,'noconversion2ir');

