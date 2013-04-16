SOFAstart;
%% Define the filename of the SOFA HRTFs
HRTFfilename='ARI NH4 hrtf_M_dtf 256';
HRTFfilename='CIPIC subject_003 hrir_final';
HRTFfilename='LISTEN 1002 IRC_1002_C_HRIR';
HRTFfilename='MIT KEMAR normal pinna';
HRTFfilename='KEMAR horizontal only normal pinna resolution 0.5 deg';

%% Load
% Define the path to the HRTF repository
databasepath=[pwd filesep '..' filesep 'HRTFs' filesep 'SOFA'];
% Load the SOFA object
Obj=SOFAload([databasepath filesep HRTFfilename]);

%% Create an input signal
in=randn(5*Obj.Data.SamplingRate,1);	% Five seconds of noise

%% Spatialize
[out,azi,ele]=SOFAspat(in,Obj,[-45 45 0 ],[0 0 -30 90]);

%% Output
% Plot the trajectory
subplot(2,1,1); plot(azi);
subplot(2,1,2); plot(ele);
% Play the sound - use headphones!
wavplay(out,Obj.Data.SamplingRate);