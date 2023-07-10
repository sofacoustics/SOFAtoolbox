function Obj=SOFAconvertUMA2SOFA(UMA)
%SOFAconvertUMA2SOFA - converts from UMA to SOFA format
%   Usage: OBJ=SOFAconvertARI2SOFA(UMA) 
% 
%   SOFAconvertUMA2SOFA converts audio data in UMA (University of Malaga) mat format and converts it to a SOFA object.
%
%   Input parameters:
%     UMA : Audio data in UMA format
% 
%   Output parameters:
%     Obj : New SOFA object (SOFA format)

% #Author: Michael Mihocic (07.07.2023)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Get an empty conventions structure
Obj = SOFAgetConventions('AnnotatedReceiverAudio');

%% Fill data with data

if isfield(UMA.Data,'Sample') 
    Obj.Data.Receiver=UMA.Data.Sample; % old versions
elseif isfield(UMA.Data,'Receiver')
    Obj.Data.Receiver=UMA.Data.Receiver;
end
Obj.Data.SamplingRate=UMA.Data.SamplingRate;
Obj.Data.SamplingRate_Units=UMA.Data.SamplingRate_Units;

Obj.EmitterPosition=UMA.EmitterPosition;
Obj.EmitterPosition_Type=UMA.EmitterPosition_Type;
Obj.EmitterPosition_Units=UMA.EmitterPosition_Units;
Obj.GLOBAL_ApplicationName=UMA.GLOBAL_ApplicationName;
Obj.GLOBAL_ApplicationVersion=UMA.GLOBAL_ApplicationVersion;
Obj.GLOBAL_AuthorContact=UMA.GLOBAL_AuthorContact;
Obj.GLOBAL_Comment=UMA.GLOBAL_Comment;
Obj.GLOBAL_DateCreated=UMA.GLOBAL_DateCreated;
Obj.GLOBAL_History=UMA.GLOBAL_History;
Obj.GLOBAL_License=UMA.GLOBAL_License;
Obj.GLOBAL_Organization=UMA.GLOBAL_Organization;
Obj.GLOBAL_Origin=UMA.GLOBAL_Origin;
Obj.GLOBAL_References=UMA.GLOBAL_References;
Obj.GLOBAL_RoomType=UMA.GLOBAL_RoomType;
Obj.GLOBAL_Title=UMA.GLOBAL_Title;
Obj.ListenerPosition=UMA.ListenerPosition;
Obj.ListenerPosition_Type=UMA.ListenerPosition_Type;
Obj.ListenerPosition_Units=UMA.ListenerPosition_Units;
Obj.ListenerUp=UMA.ListenerUp;
Obj.ListenerView=UMA.ListenerView;
Obj.ListenerView_Type=UMA.ListenerView_Type;
Obj.ListenerView_Units=UMA.ListenerView_Units;
Obj.M=UMA.M;
Obj.M_LongName=UMA.M_LongName;
Obj.M_Units=UMA.M_Units;
Obj.ReceiverPosition=UMA.ReceiverPosition;
Obj.ReceiverPosition_Type=UMA.ReceiverPosition_Type;
Obj.ReceiverPosition_Units=UMA.ReceiverPosition_Units;
Obj.SourcePosition=UMA.SourcePosition;
Obj.SourcePosition_Type=UMA.SourcePosition_Type;
Obj.SourcePosition_Units=UMA.SourcePosition_Units;

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

% %% Fill with some additional data
% Obj.GLOBAL_History='Converted from the UMA format';
% if size(meta.pos,2)>2, Obj=SOFAaddVariable(Obj,'MeasurementSourceAudioChannel','M',meta.pos(1:size(hM,2),3)); end
% if isfield(meta,'lat'), Obj=SOFAaddVariable(Obj,'MeasurementAudioLatency','MR',meta.lat); end
