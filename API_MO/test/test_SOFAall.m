% SOFA API - test script
% Test some of the SOFA API functionality

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

clc;

%% Test converters to SOFA
disp('********************************************');
clear all;
subjectID='NH4';
demo_ARI2SOFA
subjectID='NH2';
demo_ARI2SOFA

clear all;
demo_CIPIC2SOFA;

clear all;
subjectID='1002';
demo_LISTEN2SOFA;

clear all;
pinna='normal';
demo_MIT2SOFA;
pinna='large';
demo_MIT2SOFA;

clear all; 
radius=[0.5 1 2 3];
demo_TUBerlin2SOFA;

clear all;
if ~exist('OCTAVE_VERSION','builtin'), demo_FHK2SOFA; end

clear all;
demo_BTDEI2SOFA;

%% Test converters from SOFA
disp('********************************************');
clear all;
demo_SOFA2ARI;
% SOFAplotGeometry(Obj);

%% Test SOFAmerge and create TU-Berlin KEMAR file with multiple radii
disp('********************************************');
clear all;
demo_SOFAmerge;

%% Test SOFAload
disp('********************************************');
clear all;
demo_SOFAload;

%% Test SOFAspat, but do not play
disp('********************************************');
clear all;
dontplay=1;
demo_SOFAspat;

%% Test SOFAexpand and SOFAcompact
disp('********************************************');
clear all;
demo_SOFAexpandcompact;

%% Test SOFAsave
disp('********************************************');
clear all;
demo_SOFAsave;

%% Test convertions from SimpleFreeFieldHRIR to SimpleFreeFieldTF
disp('********************************************');
clear all;
demo_SimpleFreeFieldHRIR2TF;

%% Test SingleRoomDRIR
disp('********************************************');
clear all
demo_SingleRoomDRIROldenburg;

%% Test variables handling
disp('********************************************');
demo_SOFAvariables

%% Test plotting HRTFs
disp('********************************************');
demo_SOFAplotHRTF

%% Demo of headphone conventions
demo_HpIR
%% Test using string arrays
demo_SOFAstrings