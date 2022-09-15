function Obj = SOFAaddVariable(Obj,Name,Dim,Value)
%SOFAaddVariable - Add a user-defined variable
%
%   Usage: Obj = SOFAaddVariable(Obj, Name, Dim, Value) 
%
%   SOFAaddVariable adds a user-defined variable
%   to the SOFA structure Obj. Name must be a string with the variable name 
%   ('API', 'PRIVATE', or 'GLOBAL' are not allowed). Dim is a string 
%   describing the dimensions of the variable according to SOFA specifications. 
%   Value is the content of the variable with the size of Dim. 
%
%   In Obj, the new variable will be stored as Obj.Name and its
%   dimension will be stored as Obj.API.Dimensions.Name. User-
%   defined variables can be saved in SOFA file and thus remain in the 
%   object when loaded from a SOFA file. 
%
%   Obj = SOFAaddVariable(Obj,Name,'PRIVATE',Value) adds a private variable
%   to Obj. The private variable Name will be stored as Obj.PRIVATE.Name. 
%   Private variables will be not stored in the SOFA file and
%   arbitrary dimensions are allowed.
%
%		Note that adding variables to the Data structure is not supported as
%		user-defined variables in a Data structure are not recommended. 
%		Consider adding a variable at the global level instead, which would be more
%		clear for others.

% #Author: Piotr Majdak
% #Author: Piotr Majdak: dimension is added if not previously found. (9.8.2014)
% #Author: Michael Mihocic: header documentation updated (20.10.2021)
%
% SOFA Toolbox - function SOFAaddVariable
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

Dim=upper(Dim);
switch Dim
  case 'PRIVATE'
    Obj.PRIVATE.(Name)=Value;
  otherwise
    switch Name 
      case {'API','PRIVATE','GLOBAL'}
        error('This variable name is reserved.');
      otherwise
        if strncmp(Name,'Data.',length('Data.'))         
          % add variable to Data
          Name=Name(length('Data.')+1:end);
          Obj.Data.(Name)=Value;
          Obj.API.Dimensions.Data.(Name)=Dim;
        else
          % add variable to root
          Obj.(Name)=Value;
          Obj.API.Dimensions.(Name)=Dim;
        end
        dims=SOFAdefinitions('dimensions');
        for ii=1:length(Dim)  
          if ~isfield(dims,Dim(ii))
            error('Dimension not supported.');
          end
        end
    end
end
