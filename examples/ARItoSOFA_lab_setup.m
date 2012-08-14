% SOFA API - demo script
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% ARI to Sofa format conversion
number = {'2' '4' '5' '30'};
count = 1;
results4{1}{1} = 0;
for n_files=1:4
% load ARI .mat file
load(['NH' number{n_files} ' HRTFs.mat'])

Filename = ['NH' number{n_files}];

% convert audio channel to corresponding twist angle (setup at ARI lab)
angles = [-30 -20 -10 0 10 20 30 40 50 60 70 80 -25 -15 -5 5 15 25 35 45 55 65];
angles = pi*(angles/180); % convert to rad
for n=1:size(hM,2)
  Twist(n,1) = angles(meta.pos(n,3));
end

oData = {'Data',hM};
oDataType = {'DataType','FIR'};

oSourcePositionType = {'SourcePositionType','spherical'};
oSourceViewType = {'SourceViewType','spherical'};
oSourceUpType = {'SourceUpType','spherical'};
oTransmitterPositionType = {'TransmitterPositionType','spherical'};
oListenerPositionType = {'ListenerPositionType','spherical'};
oListenerViewType = {'ListenerViewType','spherical'};
oListenerUpType = {'ListenerUpType','spherical'};
oReceiverPositionType = {'ReceiverPositionType','spherical'};

oSamplingRate = {'SamplingRate',stimPar.SamplingRate};
oSubjectID = {'SubjectID',cellstr(stimPar.SubjectID)};
oApplicationName = {'ApplicationName','test application name'};
oApplicationVersion = {'ApplicationVersion',stimPar.Version};
oSourcePosition = {'SourcePosition',[0 0 0]};
oSourceView = {'SourceView',[0 0 1]};
oSourceUp = {'SourceUp',[0 90 1]};
oSourceRotation = {'SourceRotation',[repmat(0,size(hM,2),1) repmat(0,size(hM,2),1) Twist]};
oTransmitterPosition = {'TransmitterPosition',[90 0 5]}; % azi elev radius
oListenerPosition = {'ListenerPosition',[0 0 0]};
oListenerView = {'ListenerView',[0 0 1]};
oListenerUp = {'ListenerUp',[0 90 1]};
oListenerRotation = {'ListenerRotation',[meta.pos(:,1) meta.pos(:,2) repmat(0,size(hM,2),1)]};
oReceiverPosition = {'ReceiverPosition',zeros(1,3,2)};
oReceiverPosition{2}(:,:,1) = [90 0 2]; % right ear
oReceiverPosition{2}(:,:,2) = [270 0 2]; % left ear
oMeasurementID = {'MeasurementID',stimPar.ID};
oMeasurementParameterSourceAudioChannel = {'MeasurementParameterSourceAudioChannel',meta.pos(:,3)};
oMeasurementParameterItemIndex = {'MeasurementParameterItemIndex',meta.itemidx};
oMeasurementParameterAudioLatency = {'MeasurementParameterAudioLatency',meta.lat};
oMeasurementParameterSourceAmplitude = {'MeasurementParameterSourceAmplitude',meta.amp};
oRoomType = {'RoomType','free-field'};

%%

varargin = {oData,oDataType,oSourcePositionType,oSourceViewType,oSourceUpType, ...
    oTransmitterPositionType,oListenerPositionType,oListenerViewType, ...
    oListenerUpType,oReceiverPositionType,oSamplingRate, ...
    oSubjectID,oApplicationName,oApplicationVersion,oSourcePosition,oSourceView, ...
    oSourceUp,oSourceRotation,oTransmitterPosition,oListenerPosition,oListenerView, ...
    oListenerUp,oListenerRotation,oReceiverPosition,oMeasurementID, ...
    oMeasurementParameterSourceAudioChannel,oMeasurementParameterItemIndex, ...
    oMeasurementParameterAudioLatency,oMeasurementParameterSourceAmplitude,oRoomType};

SOFAsave(Filename,varargin); % write data to sofa file

% load data
SOFAload(Filename,'var');
% get all positions within a range of 80 to 100 degrees azimuth and
% -10 to 10 degrees elevation and roll angle
%results3 = SOFAgetID({Filename},'ListenerRotation',[90 0 0],'=',{'TargetValueRange',10},{'TargetCoordinate',1});
results3 = SOFAgetID({Filename},'ListenerRotation',[90 0 0],'=',{'TargetValueRange',[10 0 0]});

result_temp = SOFAgetData(Filename,results3); % save data from current file
for ii=1:size(result_temp,2)
  if(n_files==1) % for first file a simple assignemnt is enough (no need to merge)
    results4 = result_temp;
  else % merge existing data and data from current file
    results4{ii}{2} = cat(2,results4{ii}{2},result_temp{ii}{2});
  end
end
%count = count - 1; % revert increment of last loop execution
end
for ii=1:size(results4{1}{2},2) % display SubjectID's and ListenerRotations of all results
  SubjectID = results4{12}{2}{ii}
  ListenerRotation = results4{23}{2}{ii} % all within requested value range!
end