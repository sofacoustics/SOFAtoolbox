% SOFA API - test script
% Test some of the SOFA API functionality

% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

clc; close all;
warning('off','SOFA:upgrade');

%% Test converters TO SOFA
disp('############################################');
disp('#######   TEST CONVERTERS TO SOFA   ########');
disp('############################################');

%% demo_ARI2SOFA
disp('**************  demo_ARI2SOFA  **************');
clear all;
subjectID='NH4';
demo_ARI2SOFA
subjectID='NH2';
demo_ARI2SOFA
disp('*** Finished: demo_ARI2SOFA'); disp('      ');

%% demo_BTDEI2SOFA
disp('**************  demo_BTDEI2SOFA  **************');
clear all;
demo_BTDEI2SOFA;
disp('*** Finished: demo_BTDEI2SOFA'); disp('      ');

%% demo_CIPIC2SOFA
disp('**************  demo_CIPIC2SOFA  **************');
clear all;
demo_CIPIC2SOFA;
disp('*** Finished: demo_CIPIC2SOFA'); disp('      ');

%% demo_FHK2SOFA
clear all;
if ~exist('OCTAVE_VERSION','builtin')
    disp('**************  demo_FHK2SOFA  **************');
    demo_FHK2SOFA;
    disp('*** Finished: demo_FHK2SOFA'); disp('      ');
else
    disp('Skipped: demo_CIPIC2SOFA'); disp('      ');
end

%% demo_FreeFieldDirectivityTF
disp('**************  demo_FreeFieldDirectivityTF  **************');
clear all;
demo_FreeFieldDirectivityTF;
disp('*** Finished: demo_FreeFieldDirectivityTF'); disp('      ');

%% demo_FreeFieldHRTF
disp('**************  demo_FreeFieldHRTF  **************');
clear all;
demo_FreeFieldHRTF;
disp('*** Finished: demo_FreeFieldHRTF'); disp('      ');

%% demo_LISTEN2SOFA
disp('**************  demo_LISTEN2SOFA  **************');
clear all;
subjectID='1002';
demo_LISTEN2SOFA;
disp('*** Finished: demo_LISTEN2SOFA'); disp('      ');

%% demo_MIT2SOFA
disp('**************  demo_MIT2SOFA  **************');
clear all;
pinna='normal';
demo_MIT2SOFA;
pinna='large';
demo_MIT2SOFA;
disp('*** Finished: demo_MIT2SOFA'); disp('      ');

%% demo_SCUT2SOFA
disp('**************  demo_SCUT2SOFA  **************');
clear all;
demo_SCUT2SOFA;
disp('*** Finished: demo_SCUT2SOFA'); disp('      ');

%% demo_SimpleFreeFieldHRIR2TF
% Test convertions from SimpleFreeFieldHRIR to SimpleFreeFieldTF
disp('**************  demo_SimpleFreeFieldHRIR2TF  **************');
clear all;
demo_SimpleFreeFieldHRIR2TF;
disp('*** Finished: demo_SimpleFreeFieldHRIR2TF'); disp('      ');

%% demo_SimpleHeadphoneIR
% old name: demo_HpIR
disp('**************  demo_SimpleHeadphoneIR  **************');
clear all;
demo_SimpleHeadphoneIR;
disp('*** Finished: demo_SimpleHeadphoneIR'); disp('      ');

%% demo_SingleRoomDRIROlcldenburg
% Test SingleRoomDRIR
disp('**************  demo_SingleRoomDRIROldenburg  **************');
clear all
demo_SingleRoomDRIROldenburg;
disp('*** Finished: demo_SingleRoomDRIROldenburg'); disp('      ');

%% demo_SingleRoomMIMOSRIR
disp('**************  demo_SingleRoomMIMOSRIR  **************');
clear all
demo_SingleRoomMIMOSRIR;
disp('*** Finished: demo_SingleRoomMIMOSRIR'); disp('      ');

%% demo_TUBerlin2SOFA
disp('**************  demo_TUBerlin2SOFA  **************');
clear all; 
radius=[0.5 1 2 3];
demo_TUBerlin2SOFA;
disp('*** Finished: demo_TUBerlin2SOFA'); disp('      ');


%% Test converters FROM SOFA
disp('############################################');
disp('######   TEST CONVERTERS FROM SOFA   #######');
disp('############################################');

%% demo_SOFA2ARI
disp('**************  demo_SOFA2ARI  **************');
clear all;
demo_SOFA2ARI;
% SOFAplotGeometry(Obj);
disp('*** Finished: demo_SOFA2ARI'); disp('      ');

%% demo_SOFAexpandcompact
% Test SOFAexpand and SOFAcompact
disp('**************  demo_SOFAexpandcompact  **************');
clear all;
demo_SOFAexpandcompact;
disp('*** Finished: demo_SOFAexpandcompact'); disp('      ');

%% demo_SOFAHRTF2DTF
disp('**************  demo_SOFAHRTF2DTF  **************');
clear all;
demo_SOFAHRTF2DTF;
disp('*** Finished: demo_SOFAHRTF2DTF'); disp('      ');

%% demo_SOFAload
% Test SOFAload
disp('**************  demo_SOFAload  **************');
clear all;
demo_SOFAload;
disp('*** Finished: demo_SOFAload'); disp('      ');

%% demo_SOFAmerge
% Test SOFAmerge and create TU-Berlin KEMAR file with multiple radii
disp('**************  demo_SOFAmerge  **************');
clear all;
demo_SOFAmerge;
disp('*** Finished: demo_SOFAmerge'); disp('      ');

%% demo_SOFAplotGeometry
disp('**************  demo_SOFAplotGeometry  **************');
clear all;
demo_SOFAplotGeometry;
disp('*** Finished: demo_SOFAplotGeometry'); disp('      ');

%% demo_SOFAplotHRTF
% Test plotting HRTFs
disp('**************  demo_SOFAplotHRTF  **************');
demo_SOFAplotHRTF
disp('*** Finished: demo_SOFAplotHRTF'); disp('      ');

%% demo_SOFAsave
disp('**************  SOFAsave  **************');
clear all;
demo_SOFAsave;
disp('*** Finished: SOFAsave'); disp('      ');

%% demo_SOFAspat
% Test SOFAspat, but do not play
disp('**************  demo_SOFAspat  **************');
clear all;
dontplay=1;
demo_SOFAspat;
disp('*** Finished: demo_SOFAspat'); disp('      ');

%% demo_SOFAstrings
% Test using string arrays
disp('**************  demo_SOFAstrings  **************');
demo_SOFAstrings
disp('*** Finished: demo_SOFAstrings'); disp('      ');

%% demo_SOFAvariables
% Test variables handling
disp('**************  demo_SOFAvariables  **************');
demo_SOFAvariables
disp('*** Finished: demo_SOFAvariables'); disp('      ');


%% Prologue
disp('##############################################');
disp('####   COMPLETED ALL DEMOS SUCCESSFULLY   ####');
disp('##############################################');