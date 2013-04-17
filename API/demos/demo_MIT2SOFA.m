%% Define parameters
% Two ears are available, normal and large. Select one.
pinna='normal';
pinna='large';

% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

% path of the database relative to that script
databasepath=[pwd filesep '..' filesep 'HRTFs'];

%% start SOFA
SOFAstart;
%% Define directory
MITfn=[databasepath filesep 'MIT KEMAR'];
disp(['Loading: ' MITfn ', pinna:' pinna]);

%% load and convert
Obj=SOFAconvertMIT2SOFA(MITfn,pinna);

%% save SOFA file
SOFAfn=[databasepath filesep 'SOFA' filesep 'MIT KEMAR ' pinna ' pinna.sofa'];
disp(['Saving:  ' SOFAfn]);
SOFAsave(SOFAfn, Obj, compression); 