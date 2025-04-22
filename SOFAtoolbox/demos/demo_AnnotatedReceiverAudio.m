%demo_AnnotatedReceiverAudio - Script demonstrating the usage of SOFAplotGeometry.

% #Author: Michael Mihocic: ListenerView plotted (22.04.2025)
% 
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% load a SOFA file in AnnotatedReceiverAudio
% clear;
% close all;
% SOFAstart;
database='../examples';       HRTFfilename='AnnotatedReceiverAudio_0.2.sofa';
disp(['Loading ' database ':' HRTFfilename]);
Obj=SOFAload(['db://database/' database '/' HRTFfilename]);
% Obj=SOFAload(['AnnotatedReceiverAudio_0.2.sofa']);

%% convert to spherical coordinates
[azi, ele, ~] = cart2sph(Obj.ListenerView(:,1), Obj.ListenerView(:,2), Obj.ListenerView(:,3));
azi_deg = rad2deg(azi);
ele_deg = rad2deg(ele);

%% Plot ListenerView (azi, ele)
% time = (1:length(azi))/Obj.Data.SamplingRate;
figure('Name',mfilename);
subplot(2,1,1); hold on; 
plot(Obj.M,azi_deg,'LineWidth',2); % plot azimuthal trajectory
ylabel('Azimuth (deg)');
title('AnnotatedReceiverAudio: ListenerView');

% time = (1:length(ele))/Obj.Data.SamplingRate;
subplot(2,1,2); hold on; 
% ele_time = linspace(time(1), time(end), length(ele));  % ergibt [0.01 0.025 0.04]
plot(Obj.M, ele_deg, 'LineWidth', 2);  % plot elevational trajectory
% plot(time,ele,'LineWidth',2); % plot elevational trajectory
ylabel('Elevation (deg)');
xlabel('Time (s)');