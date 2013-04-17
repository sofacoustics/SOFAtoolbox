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

%% convert from ARI to SOFA
Obj=SOFAconvertARI2SOFA(ARI.hM,ARI.meta,ARI.stimPar);

%% convert back from SOFA to ARI
[hM meta stimPar]=SOFAconvertSOFA2ARI(Obj);

%% Calculate the differences
disp(['RMS difference between the new hM and the original ARI.hM: ' num2str(sum(sum(rms(hM-ARI.hM))))]);
disp(['RMS difference between the new meta.pos and the original ARI.meta.pos: ' num2str(sum(rms(meta.pos-ARI.meta.pos)))]);

