%demo_SOFAload - Demo showing how to load SOFA files, full or partial data.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: bugs fixed (10.2021)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Path definitions
SOFAfile=fullfile(SOFAdbPath, 'database', 'tu-berlin', 'qu_kemar_anechoic_all.sofa');

%% Loading the full object
disp(['Loading full object: ' SOFAfile]);
tic;
ObjFull=SOFAload(SOFAfile);
disp(['  Elapsed time: ' num2str(toc) ' s.']);
x=whos('ObjFull');
disp(['  Memory requirements: ' num2str(round(x.bytes/1024)) ' kb']);

%% Loading metadata
disp('Loading all metadata and partial data only');
tic;
Meta=SOFAload(SOFAfile,'nodata');

%% Get index of measurements with the same directions
azi=0; ele=0;
idx=find(Meta.SourcePosition(:,1)==azi & Meta.SourcePosition(:,2)==ele);

%% Load the parts of the full objects
disp('Loading partial data only');
clear Obj
for ii=1:length(idx);
	Obj(ii)=SOFAload(SOFAfile,[idx(ii) 1]);
end
disp(['  Elapsed time: ' num2str(toc) ' s.']);
xobj=whos('Obj'); xmeta=whos('Meta');
disp(['  Memory requirements: ' num2str(round((xobj.bytes+xmeta.bytes)/1024)) ' kb']);

%% Load parts of multiple dimensions of the full object
%e.g. only left ear for source positions 0°-90° at distance 1m
tic
idxSTART=find(Meta.SourcePosition(:,1)==0 & Meta.SourcePosition(:,3)==1);
idxEND=find(Meta.SourcePosition(:,1)==90 & Meta.SourcePosition(:,3)==1);
idxCOUNT=idxEND-idxSTART+1;
disp('Loading partial data in multiple dimensions')
ObjPartMultDim=SOFAload(SOFAfile,[idxSTART idxCOUNT],'M',[1 1],'R');
disp(['  Elapsed time: ' num2str(toc) ' s.']);
xobj=whos('ObjPartMultDim'); xmeta=whos('Meta');
disp(['  Memory requirements: ' num2str(round((xobj.bytes+xmeta.bytes)/1024)) ' kb']);

%% Extract and plot the fully loaded data
IRsFull=squeeze(ObjFull.Data.IR(idx,1,:));
legFull=num2str(ObjFull.SourcePosition(idx,3));
figure('Name',mfilename);
subplot(1,2,1);
plot(IRsFull');
legend(legFull);
title(['Fully loaded data'], 'Interpreter', 'none');
xlabel(['Index (Sample Taps), fully loaded']);
ylabel('Amplitude');

%% Extract and plot the partially loaded data
IRs=zeros(length(idx), Obj(1).API.N);
for ii=1:length(idx)
	IRs(ii,:)=squeeze(Obj(ii).Data.IR(:,1,:));
	leg{ii}=num2str(Obj(ii).SourcePosition(:,3));
end
subplot(1,2,2);
plot(IRs');
legend(leg);
title(['IRs for the left ear with radius as parameter' 10 ...
			'Partially loaded data'], 'Interpreter', 'none');
xlabel(['Index (Sample Taps), partially loaded']);
ylabel('Amplitude');
