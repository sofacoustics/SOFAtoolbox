% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Define the filename of the SOFA HRTFs
if ~exist('HRTFfilename','var'); HRTFfilename='ARI_NH4_hrtf_M_dtf 256'; end;
% HRTFfilename='CIPIC_subject_003_hrir_final';
% HRTFfilename='LISTEN_1002_IRC_1002_C_HRIR';
% HRTFfilename='MIT_KEMAR_normal_pinna';
% HRTFfilename='KEMAR_horizontal_only_normal_pinna_resolution_0.5_deg';
% HRTFfilename='TU-Berlin_QU_KEMAR_anechoic_radius_0.5m';
% HRTFfilename='TU-Berlin_QU_KEMAR_anechoic_radius_3m';
% HRTFfilename='TU-Berlin_QU_KEMAR_anechoic_radius_0.5_1_2_3_m';

%% Define the trajectory
azi=[-45 90 0];	% azimuth angles in degrees. If negative values are found, navigational system (-180;+180) will be used. 
ele=[0 0 -30 90]; %elevation angles in degrees

%% Load the HRTFs
f=filesep;
fullfn=[SOFAdbPath f 'SOFA' f HRTFfilename '.sofa'];
disp(['Loading ' fullfn]);
Obj=SOFAload(fullfn);

%% Create an input signal
in=randn(5*Obj.Data.SamplingRate,1);	% Five seconds of noise

%% Spatialize
[out,azi,ele,idx]=SOFAspat(in,Obj,azi,ele);
disp('Binaural signal rendered');

%% Plot the trajectories
subplot(2,1,1); hold on; box on;
plot(azi); % plot the requested, resampled azimuthal trajectory
plot(Obj.APV(idx,1),'rx');
ylabel('Azimuth (deg)');
title('SOFAspat: Trajectory');

subplot(2,1,2); hold on; box on;
plot(ele); 
plot(Obj.APV(idx,2),'rx');
ylabel('Elevation (deg)');
xlabel('Time index');
legend({'Requested', 'Actual'},'Location','Best');

%% Play the sound - use headphones!
if ~exist('dontplay','var'); wavplay(out,Obj.Data.SamplingRate); end