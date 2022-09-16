%demo_SOFAstrings - Script for testing the string array feature of SOFA.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: updated with variable ReceiverDescriptions instead of Ears (03.08.2022)
% 
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Test Strings as application-specific variable
% Load some arbritrary HRTFs
hrtf = SOFAload(fullfile(SOFAdbPath, 'database','ari','dtf_nh2.sofa'));
% Add a string array
str={};
for ii=1:hrtf.API.M
  str{ii,1}=['String' num2str(round(rand(1,1)*10000))];
end
% SOFAaddVariable(Obj,Name,Dim,Value)
hrtf2 = SOFAaddVariable(hrtf,'Test','MS',str);

% Add a new string with dimensions [RS]
strn={'left ear'; 'right ear'};
hrtf2 = SOFAaddVariable(hrtf2, 'ReceiverDescriptions', 'RS', strn);

% Save as SOFA
SOFAsave('stringtest_applicationvar.sofa',hrtf2);
% Reload the file
hrtf = SOFAload('stringtest_applicationvar.sofa');
% compare the strings
if prod(strcmp(hrtf.Test,hrtf2.Test))
    disp('SimpleFreeFieldHRIR: String Load-Reload: OK');
    delete('stringtest_applicationvar.sofa');
else
    error('String comparison showed differences');
end
clear

%% Test with conventions GeneralString (non-standardized convention, just for testing)
% Create an empty object
Obj = SOFAgetConventions('GeneralString');
% Create numeric data with M=15, R=2, N=10
Obj.Data.Double=rand(15,2,10);
% Create string arrays
str2={}; str={};
for ii=1:15
  id = num2str(round(rand(1,1)*1000000));
  str{ii,1}=['X' id];
  str2{ii,1}=['Left' id];
  str2{ii,2}=['Right' id];
end
Obj.String2 = str2;      % String1=[MRS]
Obj.Data.String1 = str;  % Data.String1=[MS]
Obj.Data.String2 = str2; % Data.String2=[MRS]

% Add a new string with dimensions [RS]
strn={'left ear'; 'right ear'};
Obj = SOFAaddVariable(Obj, 'ReceiverDescriptions', 'RS', strn);

% Update dimensions
Obj = SOFAupdateDimensions(Obj);
% Save as SOFA
SOFAsave('stringtest_generalstring.sofa',Obj);
% Reload the file
Obj2 = SOFAload('stringtest_generalstring.sofa');
% Compare the strings
if ~prod(strcmp(Obj2.Data.String2,Obj.Data.String2))
    error('Data.String2: Comparison showed differences');
end
if ~prod(strcmp(Obj2.String2,Obj.String2))
    error('String2: Comparison showed differences');
end
if ~prod(strcmp(Obj2.Data.String1,Obj.Data.String1))
    error('Data.String1: Comparison showed differences');
end
if ~prod(strcmp(Obj2.ReceiverDescriptions,Obj.ReceiverDescriptions))
    error('ReceiverDescriptions: Comparison showed differences');
end
disp('GeneralString: String1, String2, Data, ReceiverDescriptions: Load-Reload: OK');
clear
delete('stringtest_generalstring.sofa');
