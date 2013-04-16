%% Define parameters
% Subject index of the file to convert
subjectID='NH4';
% File name of the ARI file
ARIfile='hrtf_M_dtf 256';
% path of the database relative to that script
databasepath=[pwd filesep '..' filesep 'HRTFs'];
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% start SOFA
SOFAstart;
%% Load ARI file
ARIfn=[databasepath filesep 'ARI' filesep subjectID filesep ARIfile '.mat'];
disp(['Loading: ' ARIfn]);
ARI=load([ARIfn]);

%% convert
Obj=SOFAconvertARI2SOFA(ARI.hM,ARI.meta,ARI.stimPar);

%% save SOFA file
SOFAfn=[databasepath filesep 'SOFA' filesep 'ARI ' subjectID ' ' ARIfile '.sofa'];
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj, compression); 