% SOFA API - demo script
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% $LastChangedDate: $
% $LastChangedRevision: $
% $LastChangedBy: $

%% TU to Sofa format conversion
% The files can be downloaded at
% https://dev.qu.tu-berlin.de/projects/measurements/repository/show/2010-11-kemar-anechoic/mat

% === Configuration ===
infile = 'QU_KEMAR_anechoic_3m.mat';
outfile = 'QU_KEMAR_anechoic_3m';

% Load data set
load(infile);

% Get data matrix
M(:,:,1) = irs.left;
M(:,:,2) = irs.right;

oData = {'Data',M};
oDataType = {'DataType','FIR'};
oListenerRotationType = {'ListenerRotationType','spherical'};
oPositionType = {'PositionType','Cartesian'};
oSamplingRate = {'SamplingRate',irs.fs};
oSubjectID = {'SubjectID',irs.head};
oApplicationName = {'ApplicationName','Anechoic measurements, 3m'};
oApplicationVersion = {'ApplicationVersion','1.0'};
oSourcePosition = {'SourcePosition',irs.source_position};
oSourceView = {'SourceView',irs_source_reference};
oSourceUp = {'SourceUp',[0 0 1]};
oSourceRotation = {'SourceRotation',[0 0 0]};
oTransmitterPosition = {'TransmitterPosition',[0 0 0]};
oListenerPosition = {'ListenerPosition',irs.head_position};
oListenerView = {'ListenerView',irs.head_reference};
oListenerUp = {'ListenerUp',[0 0 1]};
oListenerRotation = {'ListenerRotation',[irs.head_azimuth irs.head_elevation 0]};
oReceiverPosition = {'ReceiverPosition',zeros(1,3,2)};  % is this mandatory?
%oMeasurementID = {'MeasurementID',stimPar.ID};
%oMeasurementID = {'MeasurementID',[repmat({stimPar.ID},1,1) repmat({'dfdsfasdfsad'},1,1)]};
%oMeasurementParameterSourceAudioChannel = {'MeasurementParameterSourceAudioChannel',meta.pos(:,3)};
%oMeasurementParameterItemIndex = {'MeasurementParameterItemIndex',meta.itemidx};
%oMeasurementParameterAudioLatency = {'MeasurementParameterAudioLatency',meta.lat};
%oMeasurementParameterSourceAmplitude = {'MeasurementParameterSourceAmplitude',meta.amp};
oRoomType = {'RoomType','free-field'};

% write data to SOFA file
SOFAsave(outfile, ...
    'Data', M, ...
    'DataType','FIR', ...
    'ListenerRotationType','spherical', ...
    'PositionType','Cartesian', ...
    'SamplingRate',irs.fs, ...
    'SubjectID',irs.head, ...
    'ApplicationName','Anechoic measurements, 3m', ...
    'ApplicationVersion','1.0', ...
    'SourcePosition',irs.source_position, ...
    'SourceView',irs_source_reference, ...
    'SourceUp',[0 0 1], ...
    'SourceRotation',[0 0 0], ...
    'TransmitterPosition',[0 0 0], ...
    'ListenerPosition',irs.head_position, ...
    'ListenerView',irs.head_reference, ...
    'ListenerUp',[0 0 1], ...
    'ListenerRotation',[irs.head_azimuth irs.head_elevation 0], ...
    'ReceiverPosition',zeros(1,3,2), ...
    'RoomType','free-field'
    );

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
