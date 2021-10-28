% load HRTF in FHK format and save in SOFA format

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% SOFA API - demo script
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% Octave has no support for the miro class of the FHK data
if isoctave
    warning(['demo_FHK2SOFA does not work in Octave.' newline 'Octave is not able to convert FHK to SOFA, use Matlab instead.']);
    return;
end

% load HRTF in FHK format and save as SOFA format

%% Define parameters
% Get a file name of the FHK directory
d=dir(fullfile(fileparts(SOFAdbPath),'FHK','HRIR_*.mat'));
if isempty(d)
    error(['No HRTF files available: ' fullfile(fileparts(SOFAdbPath),'FHK','HRIR_*.mat')]);
end
fn=d(1).name;
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time


%% Load FHK file
FHKfn=fullfile(fileparts(SOFAdbPath), 'FHK', fn);
disp(['Loading: ' FHKfn]);
FHK=load(FHKfn);
FHKvar=fieldnames(FHK);
FHKname=FHKvar{1};
FHKdata=FHK.(FHKname);

%% convert
Obj=SOFAconvertFHK2SOFA(FHKdata);
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA API';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');

%% save SOFA file
SOFAfn=fullfile(SOFAdbPath, 'sofa_api_mo_test', ['FHK_' FHKname '.sofa']);
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj, compression); 
