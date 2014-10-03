% SOFA API demo script
% load HRTF in BTdei format and save in SOFA format.

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Define parameters
% Headphone index of the files to convert 
hp= 'H010';
% Subject index of the files to convert
subject= 'S115';
% Measurements index of the files to convert
setm= 'Set02'; %Set01 Set02 ... 
% File name of the BTdei file
BTdeifold='COMPENSATED'; %RAW %COMPENSATED %EQUALIZED
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time


%% load BTdei file \ load database structure data
f=filesep;
BTdei_hp_add=fullfile(SOFAdbPath,'BTdei',hp,'headphones_info.mat');
BTdei_add=fullfile(SOFAdbPath,'BTdei',hp,subject,setm,BTdeifold,'MAT',[hp '-' subject '_btdei.mat']);
disp(['Loading BTdei data']);

try
    datasheet = load(BTdei_hp_add);
    BTdei.hp  = datasheet.hp_specs; 

    switch subject
      case 'S115'
        BTdei.sbjType = 'dummy head with large pinna';
      case 'S116'
        BTdei.sbjType = 'dummy head without pinna';
      case 'S117'
        BTdei.sbjType = 'dummy head without pinna';
      otherwise
        BTdei.sbjType = 'human';
    end

    container   = load(BTdei_add);
    BTdei.specs = container.specs; 
    BTdei.data  = container.data; 
catch e
	error(['Convertion - Error message: ' e.message]);
end

BTdei.type    = BTdeifold;
BTdei.typeset = setm;

%% convert
Obj = SOFAconvertBTdei2SOFA(BTdei);
Obj.GLOBAL_Comment = SOFAappendText(Obj,'GLOBAL_Comment',BTdeifold);

%% save SOFA file 
SOFAfn=fullfile(SOFAdbPath,'SOFA',['BTdei-hp_' hp '-subj_' subject '-' setm '-' BTdeifold '.sofa']);
disp(['Saving:  ' SOFAfn])
SOFAsave(SOFAfn, Obj, compression);
