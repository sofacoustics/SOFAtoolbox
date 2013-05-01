% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 


%% Define the filename of the SOFA HRTFs
HRTFfilename='ARI NH4 hrtf_M_dtf 256';
% HRTFfilename='CIPIC subject_003 hrir_final';
% HRTFfilename='LISTEN 1002 IRC_1002_C_HRIR';
HRTFfilename='MIT KEMAR normal pinna';
% HRTFfilename='KEMAR horizontal only normal pinna resolution 0.5 deg';
% HRTFfilename='TU-Berlin QU_KEMAR_anechoic_radius 0.5 m';
% HRTFfilename='TU-Berlin QU_KEMAR_anechoic_radius 3 m';
% HRTFfilename='TU-Berlin QU_KEMAR_anechoic_radius 0.5 1 2 3 m';

%% Define the trajectory
azi=[-45 90 0];	% azimuth angles in degrees. If negative values are found, navigational system (-180;+180) will be used. 
ele=[0 0 -30 90]; %elevation angles in degrees

%% Load the HRTFs
% Define the path to the HRTF repository
f=filesep;
% Load the SOFA object
Obj=SOFAload([SOFAdbPath f 'SOFA' f HRTFfilename]);

%% Create an input signal
in=randn(5*Obj.Data.SamplingRate,1);	% Five seconds of noise

%% Spatialize
[out,azi,ele,idx]=SOFAspat(in,Obj,azi,ele);

%% Plot the trajectories
subplot(2,1,1); hold on; box on;
plot(azi); % plot the requested, resampled azimuthal trajectory
plot(Obj.ListenerRotation(idx,1),'rx');
ylabel('Azimuth (deg)');
title('SOFAspat: Trajectory');

subplot(2,1,2); hold on; box on;
plot(ele); 
plot(Obj.ListenerRotation(idx,2),'rx');
ylabel('Elevation (deg)');
xlabel('Time index');
legend({'Requested', 'Actual'},'Location','Best');

%% Play the sound - use headphones!
if ~exist('dontplay','var'); wavplay(out,Obj.Data.SamplingRate); end