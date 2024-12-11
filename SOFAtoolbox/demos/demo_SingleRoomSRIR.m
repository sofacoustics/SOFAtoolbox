%demo_SingleRoomSRIR - Demonstrates the usage of the SingleRoomSRIR conventions.

% #Author: Michael Mihocic
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: bug fixed (04.12.2024)
% #Author: Michael Mihocic: bugs fixed (10.12.2024)
% #Author: Michael Mihocic: bugs fixed (11.12.2024)
% 
% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Get an empy conventions structure
clear;
conventions='SingleRoomSRIR';
disp(['Creating SOFA file with ' conventions 'conventions...']);
Obj = SOFAgetConventions(conventions);

%% Fill random data...
Obj.Data.IR=rand(4800,1);
Obj.ListenerPosition=zeros(4800,3); Obj.ListenerPosition(:,1)=1;
Obj.SourcePosition=zeros(4800,3); Obj.SourcePosition(:,2)=1;
Obj = SOFAaddVariable(Obj,'RoomCornerA','IC',[0 0 0]);
Obj = SOFAaddVariable(Obj,'RoomCornerB','IC',[1 1 1]);
Obj = SOFAaddVariable(Obj,'RoomCorners','I',0);
Obj = SOFAaddVariable(Obj,'RoomCorners_Type','S','cartesian');
Obj = SOFAaddVariable(Obj,'RoomCorners_Units','S','metre');

% Add ReceiverDescriptions as string array
str={};
for ii=1:Obj.API.R
  str{ii,1}=['String' num2str(round(rand(1,1)*10000))];
end
Obj = SOFAaddVariable(Obj,'ReceiverDescriptions','RS',str);


%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

%% Fill with attributes
Obj.GLOBAL_ListenerShortName = 'dummy';
Obj.GLOBAL_History = 'created with a demo script';
Obj.GLOBAL_DatabaseName = 'none';
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA Toolbox';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
Obj.GLOBAL_Organization = 'Acoustics Research Institute';
Obj.GLOBAL_AuthorContact = 'michael.mihocic@oeaw.ac.at';

%% save the SOFA file
SOFAfn=fullfile(SOFAdbPath,'sofatoolbox_test',[conventions '_' Obj.GLOBAL_SOFAConventionsVersion '.sofa']);
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj);
