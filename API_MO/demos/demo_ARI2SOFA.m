% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% load HRTF in ARI format and save as SOFA format

%% Define parameters
% Subject index of the file to convert
subjectID='NH2';
% File name of the ARI file
ARIfile='hrtf_M_dtf 256';
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% start SOFA
SOFAstart;
f=filesep;

%% Load ARI file
ARIfn=[SOFAdbPath f 'ARI' f subjectID f ARIfile '.mat'];
disp(['Loading: ' ARIfn]);
ARI=load(ARIfn);

%% convert
Obj=SOFAconvertARI2SOFA(ARI.hM,ARI.meta,ARI.stimPar);

%% save SOFA file
SOFAfn=[SOFAdbPath f 'SOFA' f 'ARI ' subjectID ' ' ARIfile '.sofa'];
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj, compression); 