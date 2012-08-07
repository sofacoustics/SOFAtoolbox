% SOFA API - demo script
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% ARI to Sofa format conversion
% load ARI .mat file
load('NH30 HRTFs.mat')

Filename = 'NH30';

% convert audio channel to corresponding twist angle (setup at ARI lab)
angles = [-30 -20 -10 0 10 20 30 40 50 60 70 80 -25 -15 -5 5 15 25 35 45 55 65];
angles = pi*(angles/180); % convert to rad
for n=1:size(hM,2)
  Twist(n,1) = angles(meta.pos(n,3));
end

oData = {'Data',hM};
oDataType = {'DataType','FIR'};
oListenerRotationType = {'ListenerRotationType','spherical'};
oPositionType = {'PositionType','Cartesian'};
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
oMeasurementParameterSourceAudioChannel = {'MeasurementParameterSourceAudioChannel',meta.pos(:,3)};
oMeasurementParameterItemIndex = {'MeasurementParameterItemIndex',meta.itemidx};
oMeasurementParameterAudioLatency = {'MeasurementParameterAudioLatency',meta.lat};
oMeasurementParameterSourceAmplitude = {'MeasurementParameterSourceAmplitude',meta.amp};
oRoomType = {'RoomType','free-field'};

varargin = {oData,oDataType,oListenerRotationType,oPositionType,oSamplingRate, ...
    oSubjectID,oApplicationName,oApplicationVersion,oSourcePosition,oSourceView, ...
    oSourceUp,oSourceRotation,oTransmitterPosition,oListenerPosition,oListenerView, ...
    oListenerUp,oListenerRotation,oReceiverPosition,oMeasurementID, ...
    oMeasurementParameterSourceAudioChannel,oMeasurementParameterItemIndex, ...
    oMeasurementParameterAudioLatency,oMeasurementParameterSourceAmplitude,oRoomType};
  
% global ncid; % just for debugging
% try
%   netcdf.close(ncid)
% catch
% end  

SOFAsave(Filename,varargin); % write data to sofa file

% some very basic examples who the functions might be used:
% load the entire data
results1 = SOFAload(Filename);
% get a list of SubjectIDs and corresponding ListenerRotations
results2 = SOFAlistVariables(Filename,'SubjectID','ListenerRotation');
% get all positions within a range of 80 to 100 degrees azimuth
results3 = SOFAgetID({Filename},'ListenerRotation',[90 0 0],'=',{'TargetValueRange',10});
% get data set for each position that was found by getID (previous line)
for ii=1:size(results3)
  results4{ii} = SOFAgetData(Filename,results3(ii));
end