%demo_FHK2SOFA - Load HRTF in FHK format and save in SOFA format.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% SOFA - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% Octave has no support for the miro class of the FHK data
% if isoctave
if exist('OCTAVE_VERSION','builtin') ~= 0
    warning(['demo_FHK2SOFA does not work in Octave.' newline 'Octave is not able to convert FHK to SOFA, use Matlab instead.']);
    return;
end

%% Check if miro.m class file is available, if not download file from server
% miro class might also be included in other toolboxes, eg. AKtools
if exist('miro','class') ~= 8
    % download miro.m
    disp('Downloading miro.m from TH Köln server...')
    url = 'http://audiogroup.web.th-koeln.de/FILES/miro.m'; 
    basepath=which('SOFAstart');
    basepath=basepath(1:end-12); % Kill the function name from the path.
    target=[basepath filesep 'helpers' filesep 'miro.m'];
    websave (target,url);
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
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA Toolbox';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');

%% save SOFA file
SOFAfn=fullfile(SOFAdbPath, 'sofatoolbox_test', ['FHK_' FHKname '.sofa']);
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj, compression); 
