%demo_CIPIC2SOFA - Load HRTF in CIPIC format and save in SOFA format.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Define parameters
% Subject index of the file to convert
subjectID=3;
% File name of the CIPIC file
CIPICfile='hrir_final';
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time


%% load CIPIC file
CIPICfn=fullfile(fileparts(SOFAdbPath), 'CIPIC', ['subject_' sprintf('%03d',subjectID)], [CIPICfile '.mat']);
disp(['Loading: ' CIPICfn]);
CIPIC=load(CIPICfn);

%% convert
Obj=SOFAconvertCIPIC2SOFA(CIPIC);
Obj.GLOBAL_DatabaseName = 'CIPIC';
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA Toolbox';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');

%% save SOFA file
SOFAfn=fullfile(SOFAdbPath,'sofatoolbox_test', ['CIPIC_' 'subject_' sprintf('%03d',subjectID) '_' CIPICfile '.sofa']);
disp(['Saving:  ' SOFAfn])
SOFAsave(SOFAfn, Obj, compression); 