%test_SOFAall - Test script, running all demos from 'demos' subfolder and testing some SOFA Toolbox functionalities.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: missing demos added, bugs fixed (09-10.2021)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: multiple updates for SOFAtoolbox v2.1 release (2022)
% #Author: Michael Mihocic: new demos added (23.12.2022)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

clc; close all; % clean-up first
tic; % timer
SOFAstart;
warning('off','SOFA:upgrade');
warning('off','SOFA:load');
warning('off','SOFA:save');
warning('off','SOFA:save:API');

%% Test converters TO SOFA
disp('      ');
disp('############################################');
disp('#######   TEST CONVERTERS TO SOFA   ########');
disp('############################################');

disp('      ');
disp('!!! Make sure that all source files are available in the data subdirectories. See individual readme.txt files for more information !!!');
disp('      ');

%% demo_ARI2SOFA
disp('**************  demo_ARI2SOFA  **************');
clear;
% subjectID='NH4'; % default
demo_ARI2SOFA
subjectID='NH2'; % repeat demo, overwrite default subject ID
demo_ARI2SOFA
disp('*** Finished: demo_ARI2SOFA (Output: SOFA-file(s))'); disp('      ');

%% demo_BTDEI2SOFA
disp('**************  demo_BTDEI2SOFA  **************');
clear;
demo_BTDEI2SOFA;
disp('*** Finished: demo_BTDEI2SOFA (Output: SOFA-file(s))'); disp('      ');

%% demo_CIPIC2SOFA
disp('**************  demo_CIPIC2SOFA  **************');
clear;
try
    demo_CIPIC2SOFA;
    disp('*** Finished: demo_CIPIC2SOFA (Output: SOFA-file(s))');
catch
    warning('demo_CIPIC2SOFA cannot finish successfully. Please make sure to save the source files to the \SOFAtoolbox\data\CIPIC\ directory.')
end
disp('      ');

%% demo_FHK2SOFA
clear;
% if ~exist('OCTAVE_VERSION','builtin')
disp('**************  demo_FHK2SOFA  **************');
demo_FHK2SOFA;
disp('*** Finished: demo_FHK2SOFA (Output: SOFA-file(s))'); disp('      ');
% else
%     disp('Skipped: demo_CIPIC2SOFA'); disp('      ');
% end

%% demo_LISTEN2SOFA
disp('**************  demo_LISTEN2SOFA  **************');
clear;
subjectID='1002';
try
    demo_LISTEN2SOFA;
    disp('*** Finished: demo_LISTEN2SOFA (Output: SOFA-file(s))');
catch
    warning('demo_LISTEN2SOFA cannot finish successfully. Please make sure to save the source files to the \SOFAtoolbox\data\LISTEN\ directory.')
end
disp('      ');

%% demo_MIT2SOFA
disp('**************  demo_MIT2SOFA  **************');
clear;
% pinna='normal'; % default value
demo_MIT2SOFA;
pinna='large';
demo_MIT2SOFA;
disp('*** Finished: demo_MIT2SOFA (Output: SOFA-file(s))'); disp('      ');

%% demo_SCUT2SOFA
disp('**************  demo_SCUT2SOFA  **************');
clear;
try
    demo_SCUT2SOFA;
    disp('*** Finished: demo_SCUT2SOFA (Output: SOFA-file(s))');
catch
    warning('demo_SCUT2SOFA cannot finish successfully. Please make sure to save the source files to the \SOFAtoolbox\data\SCUT\ directory.')
end
disp('      ');

%% demo_TUBerlin2SOFA
disp('**************  demo_TUBerlin2SOFA  **************');
clear;
radius=[0.5 1 2 3];
demo_TUBerlin2SOFA;
disp('*** Finished: demo_TUBerlin2SOFA (Output: SOFA-file(s))'); disp('      ');

%% demo_UMA2SOFA
disp('**************  demo_UMA2SOFA  **************');
clear;
demo_UMA2SOFA;
disp('*** Finished: demo_UMA2SOFA (Output: SOFA-file(s), Figure(s))'); disp('      ');

%% Test converters FROM SOFA
disp('############################################');
disp('######   TEST CONVERTERS FROM SOFA   #######');
disp('############################################');

%% demo_SOFA2ARI
disp('**************  demo_SOFA2ARI  **************');
clear;
demo_SOFA2ARI;
% SOFAplotGeometry(Obj);
disp('*** Finished: demo_SOFA2ARI'); disp('      ');

%% demo_SOFAHRTF2DTF
disp('**************  demo_SOFAHRTF2DTF  **************');
clear;
demo_SOFAHRTF2DTF;
disp('*** Finished: demo_SOFAHRTF2DTF (Output: Figure(s))'); disp('      ');


%% Test SOFA functions
disp('############################################');
disp('######      TEST SOFA FUNCTIONS      #######');
disp('############################################');

%% demo_SOFAcalculateITD
disp('**************  demo_SOFAcalculateITD  **************');
clear;
demo_SOFAcalculateITD;
disp('*** Finished: demo_SOFAcalculateITD (Output: Figure(s))'); disp('      ');

%% demo_SOFAcalculateLFE
disp('**************  demo_SOFAcalculateLFE  **************');
clear;
demo_SOFAcalculateLFE;
disp('*** Finished: demo_SOFAcalculateLFE (Output: Figure(s))'); disp('      ');

%% demo_SOFAexpandcompact
% Test SOFAexpand and SOFAcompact
disp('**************  demo_SOFAexpandcompact  **************');
clear;
demo_SOFAexpandcompact;
disp('*** Finished: demo_SOFAexpandcompact'); disp('      ');

%% demo_SOFAload
% Test SOFAload
disp('**************  demo_SOFAload  **************');
clear;
demo_SOFAload;
disp('*** Finished: demo_SOFAload (Output: Figure(s))'); disp('      ');

%% demo_SOFAmerge
% Test SOFAmerge and create TU-Berlin KEMAR file with multiple radii
disp('**************  demo_SOFAmerge  **************');
clear;
demo_SOFAmerge;
disp('*** Finished: demo_SOFAmerge (Output: SOFA-file(s), Figure(s))'); disp('      ');

%% demo_SOFAplotGeometry
disp('**************  demo_SOFAplotGeometry  **************');
clear;
demo_SOFAplotGeometry;
disp('*** Finished: demo_SOFAplotGeometry (Output: Figure(s))'); disp('      ');

% %% demo_plot_trumpet_directivity
% % Test plotting HRTFs
% disp('**************  demo_plot_trumpet_directivity  **************');
% demo_plot_trumpet_directivity
% disp('*** Finished: demo_plot_trumpet_directivity (Output: Figure(s))'); disp('      ');

%% demo_SOFAplotHRTF
% Test plotting HRTFs
disp('**************  demo_SOFAplotHRTF  **************');
demo_SOFAplotHRTF
disp('*** Finished: demo_SOFAplotHRTF (Output: Figure(s))'); disp('      ');

%% demo_SOFAresample
disp('**************  demo_SOFAresample  **************');
if exist('OCTAVE_VERSION','builtin')
    % Octave
    disp('*** Skipped: demo_SOFAresample, might run out of memory in Octave...'); disp('      ');
else
    % Matlab
    clear;
    demo_SOFAresample;
    disp('*** Finished: demo_SOFAresample (Output: Figure(s))'); disp('      ');
end

%% demo_SOFAsave
disp('**************  SOFAsave  **************');
clear;
demo_SOFAsave;
disp('*** Finished: SOFAsave (Output: SOFA-file(s))'); disp('      ');

%% demo_SOFAspat
% Test SOFAspat, but do not play
disp('**************  demo_SOFAspat  **************');
clear;
dontplay=1;
demo_SOFAspat;
disp('*** Finished: demo_SOFAspat (Output: Figure(s))'); disp('      ');

%% demo_SOFAstrings
% Test using string arrays
disp('**************  demo_SOFAstrings  **************');
demo_SOFAstrings
disp('*** Finished: demo_SOFAstrings (Output: SOFA-file(s))'); disp('      ');

%% demo_SOFAvariables
% Test variables handling
disp('**************  demo_SOFAvariables  **************');
demo_SOFAvariables
disp('*** Finished: demo_SOFAvariables (Output: SOFA-file(s))'); disp('      ');


%% Test SOFA conventions
disp('############################################');
disp('######     TEST SOFA CONVENTIONS     #######');
disp('############################################');

%% demo_FreeFieldDirectivityTF
disp('**************  demo_FreeFieldDirectivityTF  **************');
clear;
demo_FreeFieldDirectivityTF;
disp('*** Finished: demo_FreeFieldDirectivityTF (Output: Figure(s))'); disp('      ');

%% demo_FreeFieldHRIR
disp('**************  demo_FreeFieldHRIR  **************');
clear
demo_FreeFieldHRIR;
disp('*** Finished: demo_FreeFieldHRIR (Output: SOFA-file(s))'); disp('      ');

%% demo_FreeFieldHRTF
disp('**************  demo_FreeFieldHRTF  **************');
clear;
demo_FreeFieldHRTF;
disp('*** Finished: demo_FreeFieldHRTF (Output: SOFA-file(s), Figure(s))'); disp('      ');

%% demo_General
disp('**************  demo_General  **************');
clear
demo_General;
disp('*** Finished: General (Output: SOFA-file(s))'); disp('      ');

%% demo_GeneralFIR
disp('**************  demo_GeneralFIR  **************');
clear
demo_GeneralFIR;
disp('*** Finished: demo_GeneralFIR (Output: SOFA-file(s))'); disp('      ');

%% GeneralFIRE: outdated, use GeneralFIR-E instead

%% demo_GeneralFIR-E
% replacing GeneralFIR
disp('**************  demo_GeneralFIR_E  **************');
clear
demo_GeneralFIR_E;
disp('*** Finished: GeneralFIR_E (Output: SOFA-file(s))'); disp('      ');

%% demo_GeneralSOS
disp('**************  demo_GeneralSOS  **************');
clear
demo_GeneralSOS;
disp('*** Finished: GeneralSOS (Output: SOFA-file(s))'); disp('      ');

%% GeneralString: used in demo_SOFAstrings

%% demo_GeneralTF
disp('**************  demo_GeneralTF  **************');
clear
demo_GeneralTF;
disp('*** Finished: GeneralTF (Output: SOFA-file(s))'); disp('      ');

%% demo_GeneralTF-E
disp('**************  demo_GeneralTF_E  **************');
clear
demo_GeneralTF_E;
disp('*** Finished: GeneralTF-E (Output: SOFA-file(s))'); disp('      ');

%% demo_MultiSpeakerBRIR
disp('**************  demo_MultiSpeakerBRIR  **************');
clear
demo_MultiSpeakerBRIR;
disp('*** Finished: MultiSpeakerBRIR (Output: SOFA-file(s))'); disp('      ');

%% demo_SimpleFreeFieldHRIR2TF
% Test conversions from SimpleFreeFieldHRIR to SimpleFreeFieldHRTF
disp('**************  demo_SimpleFreeFieldHRIR2TF  **************');
clear;
demo_SimpleFreeFieldHRIR2TF;
disp('*** Finished: demo_SimpleFreeFieldHRIR2TF (Output: SOFA-file(s))'); disp('      ');

%% demo_SimpleFreeFieldHRSOS
disp('**************  demo_SimpleFreeFieldHRSOS  **************');
clear
demo_SimpleFreeFieldHRSOS;
disp('*** Finished: SimpleFreeFieldHRSOS (Output: SOFA-file(s))'); disp('      ');

%% SimpleFreeFieldHRTF
% used in function demo_SimpleFreeFieldHRIR2TF

%% demo_SimpleFreeFieldSOS
disp('**************  SimpleFreeFieldSOS  **************');
clear
demo_SimpleFreeFieldSOS;
disp('*** Finished: SimpleFreeFieldSOS (Output: SOFA-file(s))'); disp('      ');

%% demo_SimpleHeadphoneIR
% old name: demo_HpIR
disp('**************  demo_SimpleHeadphoneIR  **************');
clear;
demo_SimpleHeadphoneIR;
disp('*** Finished: demo_SimpleHeadphoneIR (Output: Figure(s))'); disp('      ');

%% demo_SingleRoomDRIR
% Test SingleRoomDRIR
disp('**************  demo_SingleRoomDRIR  **************');
clear
demo_SingleRoomDRIR;
disp('*** Finished: demo_SingleRoomDRIR (Output: SOFA-file(s))'); disp('      ');

% %% demo_SingleRoomDRIROlcldenburg (outdated)
% % Test SingleRoomDRIR
% disp('**************  demo_SingleRoomDRIROldenburg  **************');
% clear
% demo_SingleRoomDRIROldenburg;
% disp('*** Finished: demo_SingleRoomDRIROldenburg (Output: SOFA-file(s), Figure(s))'); disp('      ');

%% demo_SingleRoomMIMOSRIR
disp('**************  demo_SingleRoomMIMOSRIR  **************');
clear
demo_SingleRoomMIMOSRIR;
disp('*** Finished: demo_SingleRoomMIMOSRIR (Output: SOFA-file(s))'); disp('      ');

%% demo_SingleRoomSRIR
disp('**************  SingleRoomSRIR  **************');
clear
demo_SingleRoomSRIR;
disp('*** Finished: SingleRoomSRIR (Output: SOFA-file(s))'); disp('      ');

%% Epilogue
disp('##############################################');
disp('####   COMPLETED ALL DEMOS SUCCESSFULLY   ####');
disp('##############################################');
toc; % timer
