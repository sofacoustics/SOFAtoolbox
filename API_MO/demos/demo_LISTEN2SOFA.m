%% Define parameters
% Subject index of the file to convert
subjectID='1002';
% File name of the LISTEN file
LISTENfile=['IRC_' subjectID '_C_HRIR'];
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% start SOFA
[databasepath,f]=SOFAstart;
%% Load LISTEN file
LISTENfn=[databasepath filesep 'LISTEN' filesep 'IRC_' subjectID filesep 'COMPENSATED' filesep 'MAT' filesep 'HRIR' filesep LISTENfile '.mat'];
disp(['Loading: ' LISTENfn]);
LISTEN=load([LISTENfn]);

%% convert
Obj=SOFAconvertLISTEN2SOFA(LISTEN,subjectID);

%% save SOFA file
SOFAfn=[databasepath filesep 'SOFA' filesep 'LISTEN ' subjectID ' ' LISTENfile '.sofa'];
disp(['Saving:  ' SOFAfn]);
SOFAsave(SOFAfn, Obj, compression);