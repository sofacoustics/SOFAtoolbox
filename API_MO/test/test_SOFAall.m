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
radius=0.5;
demo_TUBerlin2SOFA;
radius=3;
demo_TUBerlin2SOFA;
radius=[0.5 1 2 3];
demo_TUBerlin2SOFA;

%% Test converters from SOFA
clear all;
demo_SOFA2ARI;

%% Test other things
clear all;
demo_SOFAload;

clear all;
dontplay=1;
demo_SOFAspat;