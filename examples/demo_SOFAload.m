%% Loading the full object
disp('------------------');
disp('Load a full object');
tic;
ObjFull=SOFAload('..\HRTFs\SOFA\TU-Berlin QU_KEMAR_anechoic_radius 0.5 1 2 3 m');
toc;
x=whos('ObjFull');
disp(['Memory requirements: ' num2str(round(x.bytes/1024)) ' kb' 10]);

%% Loading metadata
disp('---------------');
disp('Partial loading');
disp('- Loading metadata');
tic;
Meta=SOFAload('..\HRTFs\SOFA\TU-Berlin QU_KEMAR_anechoic_radius 0.5 1 2 3 m','nodata');
%% Get index of the requested direction
azi=0; ele=0;
idx=find(Meta.ListenerRotation(:,1)==azi & Meta.ListenerRotation(:,2)==ele);
%% Load the objects
disp('- Loading partial data');
for ii=1:length(idx);
	Obj(ii)=SOFAload('..\HRTFs\SOFA\TU-Berlin QU_KEMAR_anechoic_radius 0.5 1 2 3 m',[idx(ii) 1]);
end
%% Merging the objects
disp('- Merging to a single SOFA object');
ObjIdx=Obj(1);
for ii=2:length(idx)
	ObjIdx=SOFAmerge(ObjIdx,Obj(ii));
end
toc;
x=whos('ObjIdx');
disp(['Memory requirements: ' num2str(round(x.bytes/1024)) ' kb']);

%% Plot something
plot(squeeze(ObjIdx.Data.IR(:,1,:))');
legend(num2str(ObjIdx.SourcePosition(:,2)))
title('IR for the left ear with radius as parameter');
xlabel([ObjIdx.N_LongName ' (' ObjIdx.N_Units ')']);
ylabel('Amplitude');