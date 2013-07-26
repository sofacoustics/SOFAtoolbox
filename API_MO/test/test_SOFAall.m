%% Test some of the SOFA API functionality
clc;

%% Test converters to SOFA
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

%% Test converters from SOFA
clear all;
demo_SOFA2ARI;
% SOFAplotGeometry(Obj);

%% Test SOFAmerge and create TU-Berlin KEMAR file with multiple radii
clear all;
demo_SOFAmerge;

%% Test SOFAload
clear all;
demo_SOFAload;

%% Test SOFAmerge, but do not play
clear all;
dontplay=1;
demo_SOFAspat;

%% Test SOFAsave
clear all;
demo_SOFAsave;

%% Test convertions from SimpleFreeFieldHRIR to SimpleFreeFieldTF
clear all;
demo_SimpleFreeFieldHRIR2TF;

%% Test SingleRoomDRIR
clear all
demo_SingleRoomDRIROldenburg;
