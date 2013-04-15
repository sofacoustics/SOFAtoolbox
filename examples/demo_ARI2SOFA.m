%% Define parameters
% Subject index of the file to convert
filename='HRTF ARI NH2';
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% start SOFA
SOFAstart;
%% Load ARI file
ARIfn=[filename '.mat'];
disp(['Loading: ' ARIfn]);
ARI=load(ARIfn);

%% convert
Obj=ARI2SOFA(ARI.hM,ARI.meta,ARI.stimPar);

%% save SOFA file
SOFAfn=['' filename '.sofa'];
disp(['Saving:  ' SOFAfn]);
SOFAsave(SOFAfn, Obj, compression); 