% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% load HRTF in MIT format and save as SOFA format

%% Define parameters
% Two ears are available, normal and large. Select one.
pinna='normal';
pinna='large';

% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% start SOFA
SOFAstart;
f=filesep;

%% Define directory
MITfn=[SOFAdbPath f 'MIT KEMAR'];
disp(['Loading: ' MITfn ', pinna:' pinna]);

%% load and convert
Obj=SOFAconvertMIT2SOFA(MITfn,pinna);

%% save SOFA file
SOFAfn=[SOFAdbPath f 'SOFA' f 'MIT KEMAR ' pinna ' pinna.sofa'];
disp(['Saving:  ' SOFAfn]);
SOFAsave(SOFAfn, Obj, compression); 