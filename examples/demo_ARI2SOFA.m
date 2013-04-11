SOFAstart;
% Filename of the file in the ARI format
filename='HRTF ARI NH2';
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

% let's convert
ARI=load([filename '.mat']);
Obj=ARI2SOFA(ARI.hM,ARI.meta,ARI.stimPar);
SOFAsave([filename '.sofa'], Obj, compression); 