function Obj = SOFAremoveVariable(Obj,Name)
%SOFAremoveVariable - Remove a variable, its dimension, and its attributes from a SOFA object
%  Usage: Obj = SOFAremoveVariable(Obj,Name)
%
%  SOFAremoveVariable(Obj,Name) removes the variable Name
%  from the SOFA structure Obj.
%

% #Author: Piotr Majdak: adapted from SOFAaddVariable (19.06.2019)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: removing variable attributes; minor bug fixed (20.09.2022)
%
% SOFA Toolbox - function SOFAremoveVariable
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.


switch Name
    case {'API','PRIVATE','GLOBAL','Data'}
        error('This variable name is reserved.');
    otherwise
        if isfield(Obj,Name)
            Obj=rmfield(Obj,Name);
            if isfield(Obj.API.Dimensions,Name)
                Obj.API.Dimensions=rmfield(Obj.API.Dimensions,Name);
            end
            % check if fields "Name_*" are existing and remove them (variable attributes)
            ObjFields=fieldnames(Obj);

            if exist('OCTAVE_VERSION','builtin')
              % We're in Octave
              IndexC = strfind(ObjFields,[Name '_']);
              IndexContainsName = find(not(cellfun('isempty',IndexC))); % do not replace by 'contains' (not supported in Octave)
            else
              % We're in Matlab
              IndexContainsName = find(contains(ObjFields,[Name '_']));
            end

            for ii=length(IndexContainsName):-1:1 % loop is not entered if empty IndexContainsName; run backwards because of removed indices
                % check if field name starts with "Name_*"
                if startsWith(char(ObjFields(IndexContainsName(ii))),Name)
%                     disp(['removing field: ' char(ObjFields(IndexContainsName(ii)))]) % debug
                    Obj=rmfield(Obj,ObjFields(IndexContainsName(ii)));
                end
            end

        end
end
