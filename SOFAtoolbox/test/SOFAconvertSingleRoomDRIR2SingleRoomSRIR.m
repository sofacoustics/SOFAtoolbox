function ObjNew=SOFAconvertSingleRoomDRIR2SingleRoomSRIR(Obj)
%SOFAconvertSingleRoomDRIR2SingleRoomSRIR - converts SOFA objects from SingleRoomDRIR to SingleRoomSRIR convention
%   Usage: Obj=SOFAconvertSingleRoomDRIR2SingleRoomSRIR(Obj)
% 
%   SOFAconvertSOFA2ARI(Obj) converts and upgrades a SOFA object from deprecated SingleRoomDRIR convention to recommended, standardized SingleRoomSRIR convention. 
%     The script was created based on the existing data so far. Deprecated conventions data should not be created anymore.
%
%   Input parameters:
%     Obj : SOFA object (SOFA format) in SingleRoomDRIR convention
% 
%   Output parameters:
%     Obj : SOFA object (SOFA format) in SingleRoomSRIR convention

% #Author: Michael Mihocic: first version of converter (11.2022)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Get an empy conventions structure
ObjNew=SOFAgetConventions('SingleRoomSRIR');

%% Transfer global objects
ObjNew.GLOBAL_ListenerShortName = Obj.GLOBAL_ListenerShortName;          
ObjNew.GLOBAL_ApplicationName = Obj.GLOBAL_ApplicationName;
ObjNew.GLOBAL_ApplicationVersion = Obj.GLOBAL_ApplicationVersion;
ObjNew.GLOBAL_AuthorContact = Obj.GLOBAL_AuthorContact;
ObjNew.GLOBAL_Comment = Obj.GLOBAL_Comment;
% ObjNew.GLOBAL_DataType % always FIR
ObjNew.GLOBAL_History = SOFAappendText(Obj, 'GLOBAL_History', 'Converted from SingleRoomDRIR.'); % append information
ObjNew.GLOBAL_License = Obj.GLOBAL_License;
ObjNew.GLOBAL_Organization = Obj.GLOBAL_Organization;
ObjNew.GLOBAL_References = Obj.GLOBAL_References;
ObjNew.GLOBAL_RoomType = Obj.GLOBAL_RoomType; % shoebox or dae
ObjNew.GLOBAL_Origin = Obj.GLOBAL_Origin;
ObjNew.GLOBAL_Title = Obj.GLOBAL_Title;
ObjNew.GLOBAL_DatabaseName = Obj.GLOBAL_DatabaseName;
ObjNew.GLOBAL_RoomDescription = Obj.GLOBAL_RoomDescription;

%% Transfer data
ObjNew.Data=Obj.Data;

%% Upgrade some objects
ObjNew.ListenerPosition=repmat(Obj.ListenerPosition(1,:),size(ObjNew.Data.IR,1),1);
% for ii=1:size(ObjNew.Data.IR,1)
%     ObjNew.ListenerPosition(ii,:)=Obj.ListenerPosition(1,:);
% end

%% Transfer other objects
ObjNew.ListenerPosition_Type=Obj.ListenerPosition_Type;
ObjNew.ListenerPosition_Units=Obj.ListenerPosition_Units;
ObjNew.ReceiverPosition=Obj.ReceiverPosition;
ObjNew.ReceiverPosition_Type=Obj.ReceiverPosition_Type;
ObjNew.ReceiverPosition_Units=Obj.ReceiverPosition_Units;
ObjNew.SourcePosition=Obj.SourcePosition;
ObjNew.SourcePosition_Type=Obj.SourcePosition_Type;
ObjNew.SourcePosition_Units=Obj.SourcePosition_Units;
ObjNew.EmitterPosition=Obj.EmitterPosition;
ObjNew.EmitterPosition_Type=Obj.EmitterPosition_Type;
ObjNew.EmitterPosition_Units=Obj.EmitterPosition_Units;
ObjNew.ListenerUp=Obj.ListenerUp;
ObjNew.ListenerView=Obj.ListenerView;
ObjNew.ListenerView_Type=Obj.ListenerView_Type;
ObjNew.SourceUp=Obj.SourceUp;
ObjNew.SourceView=Obj.SourceView;
ObjNew.SourceView_Type=Obj.SourceView_Type;
if isfield(Obj,'RoomCornerA')
    if size(Obj.RoomCornerA,1) == 3
        ObjNew.RoomCornerA=Obj.RoomCornerA'; % must be transposed for the data existing so far
    else
        ObjNew.RoomCornerA=Obj.RoomCornerA;
    end
    ObjNew.RoomCornerA_Type=Obj.RoomCornerA_Type;
    ObjNew.RoomCornerA_Units=Obj.RoomCornerA_Units;
end

if isfield(Obj,'RoomCornerB')
    if size(Obj.RoomCornerB,1) == 3
        ObjNew.RoomCornerB=Obj.RoomCornerB'; % must be transposed for the data existing so far
    else
        ObjNew.RoomCornerB=Obj.RoomCornerB;
    end
    ObjNew.RoomCornerB_Type=Obj.RoomCornerB_Type;
    ObjNew.RoomCornerB_Units=Obj.RoomCornerB_Units;
end

%% Update dimensions
ObjNew=SOFAupdateDimensions(ObjNew);

