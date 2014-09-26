%% Test some of the SOFA API functionality
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
demo_FHK2SOFA;

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

%% Test using string arrays
demo_SOFAstrings