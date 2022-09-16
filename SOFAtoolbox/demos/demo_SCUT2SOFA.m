%demo_SCUT2SOFA - Load HRTF in MIT format and save as SOFA format.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Define parameters
% Two ears are available: normal and large. Select one.
pinna='normal';

% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time


%% Define directory
SCUTdata = 'nearfield';
SCUTroot=fullfile(fileparts(SOFAdbPath),'SCUT',SCUTdata);
disp(['Loading: ' SCUTroot]);

%% Define radii for converting
radius=[0.2 0.25 0.3:0.1:1]; 
% radius=0.2;

%% load and convert
Obj=SOFAconvertSCUT2SOFA(SCUTroot,radius);
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA Toolbox';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');

%% save SOFA file
str=sprintf('%g,',radius);
SOFAfn=fullfile(SOFAdbPath,'sofatoolbox_test',['SCUT_KEMAR_radius_' str(1:end-1) '.sofa']);
disp(['Saving:  ' SOFAfn]);
SOFAsave(SOFAfn, Obj, compression); 
