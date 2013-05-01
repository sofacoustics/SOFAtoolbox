% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% Path definitions
f=filesep;
dbpath=[SOFAdbPath f 'SOFA'];
SOFAfile=[dbpath f 'TU-Berlin_QU_KEMAR_anechoic_radius_0.5_1_2_3_m'];

%% Loading the full object
disp('Loading a full object');
tic;
f=filesep; 
ObjFull=SOFAload(SOFAfile);
disp(['  Elapsed time: ' num2str(toc) ' s.']);
x=whos('ObjFull');
disp(['  Memory requirements: ' num2str(round(x.bytes/1024)) ' kb']);

%% Loading metadata
disp('Partial loading');
disp('  Loading metadata only');
tic;
Meta=SOFAload(SOFAfile,'nodata');
%% Get index of the requested direction
azi=0; ele=0;
idx=find(Meta.ListenerRotation(:,1)==azi & Meta.ListenerRotation(:,2)==ele);
%% Load the objects
disp('  Loading partial data only');
clear Obj
for ii=1:length(idx);
	Obj(ii)=SOFAload(SOFAfile,[idx(ii) 1]);
end
%% Merging the objects
disp('  Merging to a single SOFA object');
ObjIdx=Obj(1);
for ii=2:length(idx)
	ObjIdx=SOFAmerge(ObjIdx,Obj(ii));
end
disp(['  Elapsed time: ' num2str(toc) ' s.']);
x=whos('ObjIdx');
disp(['  Memory requirements: ' num2str(round(x.bytes/1024)) ' kb']);

%% Plot something
plot(squeeze(ObjIdx.Data.IR(:,1,:))');
legend(num2str(ObjIdx.SourcePosition(:,2)))
title('IR for the left ear with radius as parameter');
xlabel([ObjIdx.N_LongName ' (' ObjIdx.N_Units ')']);
ylabel('Amplitude');