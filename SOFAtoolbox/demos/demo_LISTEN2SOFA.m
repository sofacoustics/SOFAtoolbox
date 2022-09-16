%demo_LISTEN2SOFA - Load HRTF in LISTEN format and save as SOFA format.

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
% Subject index of the file to convert
if ~exist('subjectID','var'); subjectID='1002'; end
% File name of the LISTEN file
LISTENfile=['IRC_' subjectID '_C_HRIR'];
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time


%% Load LISTEN file
LISTENfn=fullfile(fileparts(SOFAdbPath),'LISTEN',['IRC_' subjectID],'COMPENSATED','MAT','HRIR',[LISTENfile '.mat']);
if isfile(LISTENfn)
    disp(['Loading: ' LISTENfn]);
    LISTEN=load(LISTENfn);
else
    warning(['File not existing: ' LISTENfn '  -->  Please download it to: ' fullfile(fileparts(SOFAdbPath),'LISTEN',['IRC_' subjectID],'COMPENSATED','MAT','HRIR')]);
    error(['Sorry.... ' mfilename ' cannot complete! The following file is missing: ' LISTENfile '.mat']);
end

%% convert
Obj=SOFAconvertLISTEN2SOFA(LISTEN,subjectID);
Obj.GLOBAL_DatabaseName = 'LISTEN';
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA Toolbox';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');

%% save SOFA file
SOFAfn=fullfile(SOFAdbPath,'sofatoolbox_test',['LISTEN_' subjectID '_' LISTENfile '.sofa']);
disp(['Saving:  ' SOFAfn]);
SOFAsave(SOFAfn, Obj, compression);