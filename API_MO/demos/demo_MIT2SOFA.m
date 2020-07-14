% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% load HRTF in MIT format and save as SOFA format

%% Define parameters
% Two ears are available: normal and large. Select one.
pinna='normal'; 

% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time


%% Define directory
MITfn=fullfile(fileparts(SOFAdbPath),'MIT KEMAR');
disp(['Loading: ' MITfn ', pinna:' pinna]);

%% load and convert
Obj=SOFAconvertMIT2SOFA(MITfn,pinna);
Obj.GLOBAL_DatabaseName = 'MIT';
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA API';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');

%% save SOFA file
SOFAfn=fullfile(SOFAdbPath, 'sofa_api_mo_test', ['MIT_KEMAR_' pinna '_pinna.sofa']);
disp(['Saving:  ' SOFAfn]);
SOFAsave(SOFAfn, Obj, compression); 