% Script for testing the string array feature of SOFA


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
clear all


%% Test with conventions GeneralString
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
Obj = SOFAaddVariable(Obj, 'Ears', 'RS', strn);
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
if ~prod(strcmp(Obj2.Ears,Obj.Ears))
    error('Ears: Comparison showed differences');
end
disp('GeneralString: String1, String2, Data, Ears: Load-Reload: OK');
clear all
delete('stringtest_generalstring.sofa');
