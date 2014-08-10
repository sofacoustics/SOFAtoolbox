%% Test Strings as application-specific variable
% load any HRTFs
X=SOFAload([SOFAdbPath '\SOFA\ARI_NH2_hrtf_M_dtf 256.sofa']);
% add a string
str={};
for ii=1:X.API.M
  str{ii,1}=['NH' num2str(round(rand(1,1)*10000))];
end
X2=SOFAaddVariable(X,'Test','MS',str);    
% save as SOFA
SOFAsave('stringtest_applicationvar.sofa',X2);
% reload the file
N=SOFAload('stringtest_applicationvar.sofa');
% compare the strings
if prod(strcmp(N.Test,X2.Test)), 
  disp('SimpleFreeFieldHRIR: String Load-Reload: OK');
  delete('stringtest_applicationvar.sofa');
else
  error('String comparison showed differences');
end

%% Test with conventions GeneralString
% create an empty object
X=SOFAgetConventions('GeneralString');
% create a numeric data with M=15, R=2, N=10
X.Data.Double=rand(15,2,10);
% create string arrays
str2={}; str={};
for ii=1:15
  id=num2str(round(rand(1,1)*1000000));
  str{ii,1}=['X' id];
  str2{ii,1}=['L' id];
  str2{ii,2}=['R' id];
end
X.String2=str2; % String1=[MRS]
X.Data.String1=str; % Data.String1=[MS]
X.Data.String2=str2;  % Data.String2=[MRS]
% update dimensions
X2=SOFAupdateDimensions(X);
% save as SOFA
SOFAsave('stringtest_generalstring.sofa',X2);
% reload the file
N=SOFAload('stringtest_generalstring.sofa');
% compare the strings
if ~prod(strcmp(N.Data.String2,X2.Data.String2)), 
  error('Data.String2: Comparison showed differences');
end
if ~prod(strcmp(N.String2,X2.String2)), 
  error('String2: Comparison showed differences');
end
if ~prod(strcmp(N.Data.String1,X2.Data.String1)), 
  error('Data.String1: Comparison showed differences');
end
disp('GeneralString: String1, String2: Load-Reload: OK');
delete('stringtest_generalstring.sofa');


