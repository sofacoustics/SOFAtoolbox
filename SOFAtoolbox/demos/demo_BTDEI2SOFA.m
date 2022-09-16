%demo_BTDEI2SOFA - Load HRTF in BT-DEI format and save in SOFA format.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Define parameters
% Headphone index of the files to convert 
hp= 'H010';
% Subject index of the files to convert
subject= 'S115';
% Measurement index of the files to convert
setm= 'Set02'; %Set01 Set02 ... 
% File name of the BTDEI file
BTDEIfold='COMPENSATED'; %RAW %COMPENSATED %EQUALIZED
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time


%% load BTDEI file \ load database structure data
f=filesep;
BTDEI_hp_add=fullfile(fileparts(SOFAdbPath),'BTDEI',hp,'headphones_info.mat');
BTDEI_add=fullfile(fileparts(SOFAdbPath),'BTDEI',hp,subject,setm,BTDEIfold,'MAT',[hp '_' subject '_btdei.mat']);
disp(['Loading BT-DEI data']);

try
    datasheet = load(BTDEI_hp_add);
    BTDEI.hp  = datasheet.hp_specs; 

    switch subject
      case 'S115'
        BTDEI.sbjType = 'dummy head with large pinna';
      case 'S116'
        BTDEI.sbjType = 'dummy head without pinna';
      case 'S117'
        BTDEI.sbjType = 'dummy head without pinna';
      otherwise
        BTDEI.sbjType = 'human';
    end

    container   = load(BTDEI_add);
    BTDEI.specs = container.specs; 
    BTDEI.data  = container.data; 
catch e
	error(['Load BTDEI file - Error message: ' e.message ' Try downloading the BT-DEI database from: http://padva.dei.unipd.it/?page_id=345 to the corresponding directory.']);
end

BTDEI.type    = BTDEIfold;
BTDEI.typeset = setm;

%% convert
Obj = SOFAconvertBTDEI2SOFA(BTDEI);
Obj.GLOBAL_Comment = SOFAappendText(Obj,'GLOBAL_Comment',BTDEIfold);

%% save SOFA file 
SOFAfn=fullfile(SOFAdbPath,'sofatoolbox_test',['BTDEI-hp_' hp '_subj_' subject '-' setm '-' BTDEIfold '.sofa']);
disp(['Saving:  ' SOFAfn])
SOFAsave(SOFAfn, Obj, compression);
