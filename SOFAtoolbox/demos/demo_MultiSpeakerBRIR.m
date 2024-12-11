%demo_MultiSpeakerBRIR - Demonstrates the upgrade from MultiSpeakerBRIR to SingleRoomMIMOSRIR.
% The demo script downloads the file from the SOFA database,
% upgrades the conventions and saves it as a new file in the local database.

% #Author: Michael Mihocic: (23.12.2022)
% #Author: Michael Mihocic: SOFAload: no checks added, source convention outdated (10.07.2023)
% #Author: Michael Mihocic: bugs fixed (10.12.2024)
% #Author: Michael Mihocic: bugs fixed (11.12.2024)
%
% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% Load SOFA file
db='thk';
fn='BRIR_CR1_KU_MICS_L.sofa';
Obj=SOFAload(['db://' fullfile('database',db,fn)],'nochecks');

%% Upgrade convention
disp(['Old conventions: ' Obj.GLOBAL_SOFAConventions ' v' Obj.GLOBAL_SOFAConventionsVersion])
[Obj,modified] = SOFAupgradeConventions(Obj);
disp(['New conventions: ' Obj.GLOBAL_SOFAConventions ' v' Obj.GLOBAL_SOFAConventionsVersion])

%% Fix some fields, add room parameters
Obj.GLOBAL_RoomType = 'shoebox';
Obj = SOFAaddVariable(Obj,'RoomCornerA','IC',[0 0 0]);
Obj = SOFAaddVariable(Obj,'RoomCornerB','IC',[3 3 3]);
Obj = SOFAaddVariable(Obj,'RoomCorners','I',0);
Obj = SOFAaddVariable(Obj,'RoomCorners_Type','S','cartesian');
Obj = SOFAaddVariable(Obj,'RoomCorners_Units','S','metre');

%% save the SOFA file
SOFAfn=fullfile([SOFAdbPath,'\','sofatoolbox_test','\',fn]);
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj);