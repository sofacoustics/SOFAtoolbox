% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% load HRTF in LISTEN format and save as SOFA format

%% Define parameters
% Subject index of the file to convert
if ~exist('subjectID','var'); subjectID='1002'; end;
% File name of the LISTEN file
LISTENfile=['IRC_' subjectID '_C_HRIR'];
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time


%% Load LISTEN file
f=filesep;
LISTENfn=[SOFAdbPath f 'LISTEN' f 'IRC_' subjectID f 'COMPENSATED' f 'MAT' f 'HRIR' f LISTENfile '.mat'];
disp(['Loading: ' LISTENfn]);
LISTEN=load(LISTENfn);

%% convert
Obj=SOFAconvertLISTEN2SOFA(LISTEN,subjectID);
Obj.GLOBAL_DatabaseName = 'LISTEN';
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA API';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');

%% save SOFA file
SOFAfn=[SOFAdbPath f 'SOFA' f 'LISTEN_' subjectID '_' LISTENfile '.sofa'];
disp(['Saving:  ' SOFAfn]);
SOFAsave(SOFAfn, Obj, compression);