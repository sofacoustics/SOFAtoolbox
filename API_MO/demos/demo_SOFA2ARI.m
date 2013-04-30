% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% load HRTF in SOFA format and save as ARI format

%% Define parameters
% Subject index of the file to convert
subjectID='NH4';
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

%% convert from ARI to SOFA
disp('Converting to SOFA...');
Obj=SOFAconvertARI2SOFA(ARI.hM,ARI.meta,ARI.stimPar);

%% convert back from SOFA to ARI
disp('Converting back to ARI (hM, meta, stimPar)...');
[hM, meta, stimPar]=SOFAconvertSOFA2ARI(Obj);

%% Calculate the differences
disp(['RMS difference between the new hM and the original ARI.hM: ' num2str(sum(sum(rms(hM-ARI.hM))))]);
disp(['RMS difference between the new meta.pos and the original ARI.meta.pos: ' num2str(sum(rms(meta.pos-ARI.meta.pos)))]);

