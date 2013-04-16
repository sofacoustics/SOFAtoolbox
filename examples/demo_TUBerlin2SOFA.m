%% Define parameters
% Prefix to the files 
TUBfile = 'QU_KEMAR_anechoic_';
% Define radii to be loaded (per default 0.5, 1, 2, and 3 m are available)
radius=[0.5 1 2 3];
% radius=0.5;
% radius=3;

% path of the database relative to that script
databasepath=[pwd filesep '..' filesep 'HRTFs'];
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% start SOFA
SOFAstart;
%% Load and convert the requested TU-Berlin files
for ii=1:length(radius)
	TUBfn=[databasepath filesep 'TU-Berlin KEMAR' filesep TUBfile num2str(radius(ii)) 'm.mat'];
	disp(['Loading: ' TUBfn]);
	TUB=load(TUBfn);
	Obj(ii)=SOFAconvertTUBerlin2SOFA(TUB.irs);
end

%% merge the different radii to a single SOFA object
ObjFull=Obj(1);
if length(radius)>1,
	for ii=2:length(radius)
		ObjFull=SOFAmerge(ObjFull,Obj(ii));
	end
end

%% save the object as a single SOFA file
SOFAfn=[databasepath filesep 'SOFA' filesep 'TU-Berlin ' TUBfile 'radius ' sprintf('%g ',radius) 'm.sofa'];
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, ObjFull, compression);