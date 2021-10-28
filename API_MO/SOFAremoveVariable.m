function Obj = SOFAremoveVariable(Obj,Name)
%SOFAremoveVariable
%   Obj = SOFAremoveVariable(Obj,Name) removes the user-defined variable
%   from the SOFA structure OBJ. NAME must be a string with the variable name 
%   ('API', 'PRIVATE', or 'GLOBAL' are not allowed). 
%
%

% #Author: Piotr Majdak: adapted from SOFAaddVariable (19.06.2019)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
%
% SOFA API - function SOFAremoveVariable
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


switch Name 
    case {'API','PRIVATE','GLOBAL','PRIVATE','Data'}
        error('This variable name is reserved.');
    otherwise
        if isfield(Obj,Name)
          Obj=rmfield(Obj,Name);
          if isfield(Obj.API.Dimensions,Name)
              Obj.API.Dimensions=rmfield(Obj.API.Dimensions,Name);
          end
       end
    end
end
