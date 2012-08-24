% SOFA API - demo script
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% ARI to Sofa format conversion
clear
% load ARI .mat file
load('HRTF ARI NH4.mat')

Filename = 'NH4';

% convert audio channel to corresponding twist angle (setup at ARI lab)
angles = [-30 -20 -10 0 10 20 30 40 50 60 70 80 -25 -15 -5 5 15 25 35 45 55 65];
for n=1:size(hM,2)
  Twist(n,1) = angles(meta.pos(n,3));
end

 % compose data matrix from hM
for m=1:size(hM,2)
  for r=1:size(hM,3)
    for n=1:size(hM,1)
      data.FIR(m,r,n) = hM(n,m,r);
    end
  end
end

oData = {'Data',data};
oDataType = {'DataType','FIR'};

oSourcePositionType = {'SourcePositionType','cartesian'};
oSourceViewType = {'SourceViewType','cartesian'};
oSourceUpType = {'SourceUpType','cartesian'};
oTransmitterPositionType = {'TransmitterPositionType','cartesian'};
oListenerPositionType = {'ListenerPositionType','cartesian'};
oListenerViewType = {'ListenerViewType','cartesian'};
oListenerUpType = {'ListenerUpType','cartesian'};
oReceiverPositionType = {'ReceiverPositionType','cartesian'};

oSamplingRate = {'SamplingRate',stimPar.SamplingRate};
oSubjectID = {'SubjectID',cellstr(stimPar.SubjectID)};
oApplicationName = {'ApplicationName','test application name'};
oApplicationVersion = {'ApplicationVersion',stimPar.Version};
oSourcePosition = {'SourcePosition',[0 0 0]};
oSourceView = {'SourceView',[1 0 0]};
oSourceUp = {'SourceUp',[0 0 1]};
oSourceRotation = {'SourceRotation',[0 0 0]};
oTransmitterPosition = {'TransmitterPosition',[0 5 0]};
oListenerPosition = {'ListenerPosition',[0 0 0]};

oListenerView = {'ListenerView',[0 -1 0]};
oListenerUp = {'ListenerUp',[0 0 1]};
oListenerRotation = {'ListenerRotation',[meta.pos(:,1) meta.pos(:,2) Twist]};
oReceiverPosition = {'ReceiverPosition',zeros(1,3,2)};
oReceiverPosition{2}(:,:,1) = [0 1 0];
oReceiverPosition{2}(:,:,2) = [0 -1 0];
oMeasurementID = {'MeasurementID',stimPar.ID};
%oMeasurementID = {'MeasurementID',[repmat({stimPar.ID},1550,1) repmat({'dfdsfasdfsad'},1550,1)]};
%oMeasurementID{2}(1,1) = {'short'};
oMeasurementParameterSourceAudioChannel = {'MeasurementParameterSourceAudioChannel',meta.pos(:,3)};
oMeasurementParameterItemIndex = {'MeasurementParameterItemIndex',meta.itemidx};
oMeasurementParameterAudioLatency = {'MeasurementParameterAudioLatency',meta.lat};
oMeasurementParameterSourceAmplitude = {'MeasurementParameterSourceAmplitude',meta.amp};
oRoomType = {'RoomType','free-field'};

Subject.Data = data;

Subject.DataType = 'FIR';
Subject.SourcePositionType = 'cartesian';
Subject.SourceUpType = 'cartesian';
Subject.SourceViewType = 'cartesian';
Subject.TransmitterPositionType = 'cartesian';
Subject.ListenerPositionType = 'cartesian';
Subject.ListenerViewType = 'cartesian';
Subject.ListenerUpType = 'cartesian';
Subject.ReceiverPositionType = 'cartesian';

Subject.SamplingRate = stimPar.SamplingRate;
Subject.SubjectID = stimPar.SubjectID;
Subject.ApplicationName = 'test application name';
Subject.ApplicationVersion = stimPar.Version;
Subject.SourcePosition = [0 0 0];
Subject.SourceView = [1 0 0];
Subject.SourceUp = [0 0 1];
Subject.SourceRotation = [0 0 0];
Subject.TransmitterPosition = [0 5 0];
Subject.ListenerPosition = [0 0 0];
Subject.ListenerView = [0 -1 0];
Subject.ListenerUp = [0 0 1];
Subject.ListenerRotation = [meta.pos(:,1) meta.pos(:,2) Twist];
Subject.ReceiverPosition = zeros(1,3,2);
Subject.ReceiverPosition(:,:,1) = [0 1 0];
Subject.ReceiverPosition(:,:,2) = [0 0 1];
Subject.MeasurementID = stimPar.ID;
Subject.MeasurementParameterSourceAudioChannel = meta.pos(:,3);
Subject.MeasurementParameterItemIndex = meta.itemidx;
Subject.MeasurementParameterAudioLatency = meta.lat;
Subject.MeasurementParameterSourceAmplitude = meta.amp;
Subject.RoomType = 'free-field';

% if the following section is uncommented, Subject will be a cell array:
% clear Subject;
% Subject = {oData,oDataType,oSourcePositionType,oSourceViewType,oSourceUpType, ...
%     oTransmitterPositionType,oListenerPositionType,oListenerViewType, ...
%     oListenerUpType,oReceiverPositionType,oSamplingRate, ...
%     oSubjectID,oApplicationName,oApplicationVersion,oSourcePosition,oSourceView, ...
%     oSourceUp,oSourceRotation,oTransmitterPosition,oListenerPosition,oListenerView, ...
%     oListenerUp,oListenerRotation,oReceiverPosition,oMeasurementID, ...
%     oMeasurementParameterSourceAudioChannel,oMeasurementParameterItemIndex, ...
%     oMeasurementParameterAudioLatency,oMeasurementParameterSourceAmplitude,oRoomType};
  


SOFAsave(Filename,Subject,1); % write data to sofa file

% some very basic examples how the functions might be used:
% load the entire data
results1 = SOFAload(Filename);
% print a list of Metadata
SOFAlistMetadata(Filename);
% get all positions within a range of 80 to 100 degrees azimuth
results3 = SOFAparse({Filename},'ListenerRotation',[90 0 0],'=','TargetValueRange',10);
% get dataset for each position that was found by parse (previous line)
results4 = SOFAget(Filename,results3)