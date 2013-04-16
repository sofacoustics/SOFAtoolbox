SOFAstart;
%% Define the filename of the SOFA HRTFs
HRTFfilename='ARI NH4 hrtf_M_dtf 256';
% HRTFfilename='CIPIC subject_003 hrir_final';
% HRTFfilename='LISTEN 1002 IRC_1002_C_HRIR';
% HRTFfilename='MIT KEMAR normal pinna';
% HRTFfilename='KEMAR horizontal only normal pinna resolution 0.5 deg';
HRTFfilename='TU-Berlin QU_KEMAR_anechoic_radius 0.5 m';
HRTFfilename='TU-Berlin QU_KEMAR_anechoic_radius 3 m';
% HRTFfilename='TU-Berlin QU_KEMAR_anechoic_radius 0.5 1 2 3 m';

%% Define the trajectory
azi=[-45 90 0];	% azimuth angles in degrees. If negative values are found, navigational system (-180;+180) will be used. 
ele=[0 0 -30 90]; %elevation angles in degrees

%% Load the HRTFs
% Define the path to the HRTF repository
databasepath=[pwd filesep '..' filesep 'HRTFs' filesep 'SOFA'];
% Load the SOFA object
Obj=SOFAload([databasepath filesep HRTFfilename]);

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
wavplay(out,Obj.Data.SamplingRate);