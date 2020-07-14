function Obj = SOFAaddVariable(Obj,Name,Dim,Value)
%SOFAaddVariable
%   Obj = SOFAaddVariable(Obj,Name,Dim,Value) adds a user-defined variable
%   to the SOFA structure OBJ. NAME must be a string with the variable name 
%   ('API', 'PRIVATE', or 'GLOBAL' are not allowed). DIM is a string 
%   describing the dimensions of the variable according to SOFA specifications. 
%   The content of NAME is stored in VALUE which must be of the size DIM. 
%   The used-defined variable NAME will be stored as Obj.NAME and its
%   dimension will be stored as Obj.API.Dimensions.NAME. Note that user-
%   defined variables can be saved in SOFA file and thus remain in the 
%   object when loaded from a SOFA file. 
%
%   Obj = SOFAaddVariable(Obj,Name,'PRIVATE',Value) adds a private variable
%   to OBJ. The private variable NAME will be stored as Obj.PRIVATE.NAME. 
%   Note that the private variables will be not stored in SOFA files and
%   arbitrary dimensions are allowed.
%
%		Note that adding variables to Data is not supported and should not be used
%		as it might be confusing having user-defined variables in a Data structure. 
%		Consider adding a variable at the global level instead, which would be more
%		clear for others.

% 9.8.2014: dimension is added if not previously found.
%
% SOFA API - function SOFAaddVariable
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
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
