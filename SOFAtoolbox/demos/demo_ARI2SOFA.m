%demo_ARI2SOFA - Load HRTF in ARI format and save as SOFA format.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Define parameters
% Subject index of the file to convert (use default if not defined)
if ~exist('subjectID','var'); subjectID='NH4'; end
% File name of the ARI file
ARIfile='hrtf_M_dtf 256';
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time


%% Load ARI file
ARIfn=fullfile(fileparts(SOFAdbPath), 'ARI', subjectID, [ARIfile '.mat']);

if isfile(ARIfn)
    disp(['Loading: ' ARIfn]);
    ARI=load(ARIfn);
else
    warning(['File not existing: ' ARIfn '  -->  Please download it from http://www.oeaw.ac.at/isf/hrtf and save it to: ' fullfile(fileparts(SOFAdbPath), 'ARI', subjectID)]);
    error(['Sorry.... ' mfilename ' cannot complete!']);
end


%% convert
Obj=SOFAconvertARI2SOFA(ARI.hM,ARI.meta,ARI.stimPar);
Obj.GLOBAL_DatabaseName = 'ARI';
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA Toolbox';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
Obj.GLOBAL_Organization = 'Acoustics Research Institute';
Obj.GLOBAL_AuthorContact = 'piotr@majdak.com';

%% save SOFA file
SOFAfn=fullfile(SOFAdbPath,'sofatoolbox_test',['ARI_' subjectID '_' ARIfile '.sofa']);
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj, compression); 