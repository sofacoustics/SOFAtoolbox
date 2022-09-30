%demo_SOFAspat - Demo script showing how to spatialize noise with a SOFA HRTF set.

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


%% Define the filename of the SOFA HRTFs
database='ari';       HRTFfilename='hrtf_nh4.sofa';
% database='cipic';     HRTFfilename='subject_003.sofa';
% database='listen';    HRTFfilename='irc_1002.sofa';
% database='mit';       HRTFfilename='mit_kemar_normal_pinna.sofa';
% database='tu-berlin'; HRTFfilename='qu_kemar_anechoic_0.5m.sofa';
% database='tu-berlin'; HRTFfilename='qu_kemar_anechoic_all.sofa';

%% Define the trajectory
azi=[-45 90 0];	% azimuth angles in degrees. If negative values are found, navigational system (-180;+180) will be used. 
ele=[0 0 -30 90]; %elevation angles in degrees

%% Load the HRTFs
fullfn=fullfile(SOFAdbPath, 'database', database, HRTFfilename);
disp(['Loading ' fullfn]);
Obj=SOFAload(fullfn);

%% Create an input signal
in=randn(5*Obj.Data.SamplingRate,1);	% Five seconds of noise
fade=round(0.02*Obj.Data.SamplingRate); % fade in and out for 20 ms
win=hanning(fade*2);  
in(1:fade)=in(1:fade).*win(1:fade);
in(end-fade+1:end)=in(end-fade+1:end).*win(fade+1:end);
%% Spatialize
[out,azi,ele,idx]=SOFAspat(in,Obj,azi,ele);
disp('Binaural signal rendered');

%% Plot the trajectories
time = (1:length(azi))/Obj.Data.SamplingRate;

figure('Name',mfilename);
subplot(2,1,1);
plot(time,azi); % plot azimuthal trajectory
ylabel('Azimuth (deg)');
title('SOFAspat: Trajectory');

subplot(2,1,2);
plot(time,ele); % plot elevational trajectory
ylabel('Elevation (deg)');
xlabel('Time (s)');

%% Play the sound - use headphones!
if ~exist('dontplay','var'); 
  p=audioplayer(out, Obj.Data.SamplingRate);
  play(p); 
end