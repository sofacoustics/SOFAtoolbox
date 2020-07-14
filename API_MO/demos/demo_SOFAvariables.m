% SOFA API - script demonstrating the usage of variables in the API
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% load a SOFA file
SOFAfile=fullfile(SOFAdbPath,'database','ari','hrtf_nh5.sofa');
Obj=SOFAload(SOFAfile);

%% read an API internal variable
M=Obj.API.M;    % read M, the number of measurements
disp(['The number of measurements, M, is ' num2str(M)]);

%% add a user-defined variable and check it
  % create a random variable with the size M
myvar=rand(M,1);
  % add to Obj
Obj=SOFAaddVariable(Obj,'MyVariable','M',myvar);
  % check if it's there and identical to the created one
if ~all((Obj.MyVariable-myvar)<=eps('single')),  error('Error!'); end
  % check if the size is M
if ~all(size(Obj.MyVariable)==[M 1]), error('Error!'); end 
  % check if the dimensions have been correctly stored
if ~strcmp(Obj.API.Dimensions.MyVariable,'M'), error('Error!'); end

%% add a private variable
  % create a random variable with the size 1000 x 10
privatevar=rand(1000,10);
  % add to Obj as private
Obj=SOFAaddVariable(Obj,'MyVariable','Private',privatevar);
  % check if it's there and identical to the created one
if ~all((Obj.PRIVATE.MyVariable-privatevar)<=eps('single')),  error('Error!'); end
  % check if the size is 1000 x 10
if ~all(size(Obj.PRIVATE.MyVariable)==[1000 10]), error('Error!'); end 

%% Save the object
  % create a random file name
fn=[mfilename '_temp_' num2str(rand(1,1)) '.sofa'];
SOFAsave(fn,Obj);

%% Reload the object and remove the temporary file
Obj2=SOFAload(fn);
delete(fn);

%% Check if the user-defined variable is still there
if ~isfield(Obj2,'MyVariable'), error('Error!'); end 
  % check if the size is M
if ~all(size(Obj2.MyVariable)==[M 1]), error('Error!'); end 
  % check if it is identical to the created one
if ~all((Obj2.MyVariable-myvar)<=eps('single')),  error('Error!'); end
  % check if the dimensions have been correctly stored
if ~strcmp(Obj2.API.Dimensions.MyVariable,'M'), error('Error!'); end

%% Make sure that the private variable is not there!
if isfield(Obj2,'PRIVATE'), error('Error!'); end 
