deprecated!!!
% SOFA API - demo script
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% ARI to Sofa format conversion
number = {'2' '4' '5' '30'};

for n_files=1:4
% load ARI .mat file
load(['HRTF ARI NH' number{n_files} '.mat'])

Filename = ['NH' number{n_files}];

% convert audio channel to corresponding twist angle (setup at ARI lab)
angles = [-30 -20 -10 0 10 20 30 40 50 60 70 80 -25 -15 -5 5 15 25 35 45 55 65];
for n=1:size(hM,2)
  Twist(n,1) = angles(meta.pos(n,3));
end

for n=1:size(hM,1) % compose data matrix from hM
  for m=1:size(hM,2)
    for r=1:size(hM,3)
      data.FIR(m,r,n) = hM(n,m,r);
    end
  end
end

Subject.Data = data;

Subject.DataType = 'FIR';
Subject.SourcePositionType = 'spherical';
Subject.SourceUpType = 'spherical';
Subject.SourceViewType = 'spherical';
Subject.TransmitterPositionType = 'spherical';
Subject.ListenerPositionType = 'spherical';
Subject.ListenerViewType = 'spherical';
Subject.ListenerUpType = 'spherical';
Subject.ReceiverPositionType = 'spherical';

Subject.SamplingRate = stimPar.SamplingRate;
Subject.SubjectID = stimPar.SubjectID;
Subject.ApplicationName = 'test application name';
Subject.ApplicationVersion = stimPar.Version;
Subject.SourcePosition = [0 0 0];
Subject.SourceView = [0 0 1]; % "x-axis"
Subject.SourceUp = [0 90 1]; % "z-axis"
Subject.SourceRotation = [repmat(0,size(hM,2),1) repmat(0,size(hM,2),1) Twist];
Subject.TransmitterPosition = [90 0 5]; % azi elev radius
Subject.ListenerPosition = [0 0 0];
Subject.ListenerView = [0 0 1];
Subject.ListenerUp = [0 90 1];
Subject.ListenerRotation = [meta.pos(:,1) meta.pos(:,2) repmat(0,size(hM,2),1)];
Subject.ReceiverPosition = zeros(1,3,2);
Subject.ReceiverPosition(:,:,1) = [90 0 2]; % right ear
Subject.ReceiverPosition(:,:,2) = [270 0 2]; % left ear
Subject.MeasurementID = stimPar.ID;
Subject.MeasurementParameterSourceAudioChannel = meta.pos(:,3);
Subject.MeasurementParameterItemIndex = meta.itemidx;
Subject.MeasurementParameterAudioLatency = meta.lat;
Subject.MeasurementParameterSourceAmplitude = meta.amp;
Subject.RoomType = 'free-field';

SOFAsave(Filename,Subject); % write data to sofa file

% get all positions within a range of 80 to 100 degrees azimuth and
% -10 to 10 degrees elevation and roll angle
%Ids = SOFAparse({Filename},'ListenerRotation',[90 0 0],'=',{'TargetValueRange',10},{'TargetCoordinate',1});
Ids = SOFAparse({Filename},'ListenerRotation',[90 0 0],'=','TargetValueRange',[10 0 0]);

result_temp = SOFAget(Filename,Ids); % save data from current file
FieldNames = fieldnames(result_temp);
for ii=1:size(FieldNames,1)
  if(n_files==1) % for first file a simple assignemnt is enough (no need to merge)
    results = result_temp;
  else % merge existing data and data from current file
    CurrentFieldName = FieldNames{ii};
    results.(CurrentFieldName) = cat(2,results.(CurrentFieldName),result_temp.(CurrentFieldName));
  end
end
end % end of loop through files

% display SubjectID's and ListenerRotations of all results
FieldNames = fieldnames(results);
for ii=1:size(results.SubjectID,2) 
  SubjectID = results.SubjectID{ii}
  ListenerRotation = results.ListenerRotation{ii}
end