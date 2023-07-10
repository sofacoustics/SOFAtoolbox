%demo_UMA2SOFA - Load Audio data in UMA format and save as SOFA format.

% #Author: Michael Mihocic (07.07.2023)
% #Author: Michael Mihocic (10.07.2023): Figure with subject tracking data added
% 
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Define parameters
% Data compression (0..uncompressed, 9..most compressed)
% compression=1; % results in a nice compression within a reasonable processing time

%% Load UMA file
UMAfn=fullfile(fileparts(SOFAdbPath), 'UMA','AnnotatedReceiverAudio.mat');

if isfile(UMAfn)
    disp(['Loading: ' UMAfn]);
    UMA=load(UMAfn);
else
    warning(['File not existing: ' UMAfn '  -->  Please get an UMA SOFA file (AnnotatedReceiverAudio convention) from University of Malaga and save it to: ' fullfile(fileparts(SOFAdbPath), 'UMA')]);
    error(['Sorry.... ' mfilename ' cannot complete!']);
end

%% convert UMA to SOFA file
ObjUMA=SOFAconvertUMA2SOFA(UMA);

%% save SOFA file as AnnotatedReceiverAudio
SOFAfn=fullfile(SOFAdbPath,'sofatoolbox_test',['UMA_'  'AnnotatedReceiverAudio.sofa']);
disp(['Saving:  ' SOFAfn]);
ObjUMA=SOFAsave(SOFAfn, ObjUMA); 

%% Plot the subject's rotations

LVazi  = rad2deg(squeeze(ObjUMA.ListenerView(:,1)));
LVele  = rad2deg(squeeze(ObjUMA.ListenerView(:,2)));
% LVroll = squeeze(Obj.ListenerView(:,3));

figure('Name',mfilename);
time = ObjUMA.M; % (1:length(LVazi))/Obj.Data.SamplingRate;

subplot(2,1,1); hold on; 
plot(time,LVazi,'o'); % plot azimuthal trajectory
title('Tracked Subject''s Rotation');
ylabel('Azimuth (deg)');

subplot(2,1,2); hold on; 
plot(time,LVele,'o','Color','red'); % plot elevational trajectory
ylabel('Elevation (deg)');
xlabel('Time (s)');
