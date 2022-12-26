function [Obj] = SOFAsave(filename,Obj,varargin)
%SOFAsave - Save a SOFA object to a SOFA file
%   Usage: Obj = SOFAsave(fn,Obj)
%          Obj = SOFAsave(fn,Obj,compression)
%
%   SOFAsave(fn, Obj) creates a new SOFA file with the name fn and stores 
%   the SOFA object Obj in that file. Before storing, the existence of 
%   mandatory variables in Obj is checked and the dimensions are updated. 
%   Further, DateModified is set to now, DateCreated is created if 
%   not provided in Obj, and all read-only metadata is reset to its standard
%   values. The convention's version is not updated while saving. A matching
%   csv file with convention+version must be existing in the conventions 
%   folder, otherwise an error is thrown.
%
%   SOFAsave(fn,Obj,compression) specifies the amount of compression of the 
%   file as a number between 0 and 9. compression of 0 is no compression and 
%   9 is the largest compression available. For the compression, ZIP algorithm
%   is used.
%   
%   SOFAsave will output a warning if GLOBAL_API and GLOBAL_APIName does not match 
%   the actual API name and version, respectively. To supress these warnings, use
%
%      warning('off','SOFA:save:API');
%
%   before calling SOFAsave(..). 
%   
%   Obj=SOFAsave(..) returns the updated object.
%
 

% #Author: Piotr Majdak
% #Author: Michael Mihocic: doc and header documentation updated (28.10.2021)
%
% SOFA Toolbox - function SOFAsave
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

Def = SOFAdefinitions;

%% check file name
filename=SOFAcheckFilename(filename);

%% Remove private data
if isfield(Obj,'PRIVATE'), Obj=rmfield(Obj,'PRIVATE'); end

%% Check convention: mandatory variables
% try % try first without upgrade
    ObjCheck = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'m',Obj.GLOBAL_SOFAConventionsVersion);
% catch % if convention version not available try to upgrade to latest version
%     disp([Obj.GLOBAL_SOFAConventions ' ' Obj.GLOBAL_SOFAConventionsVersion]);
%     ObjCheck = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'m');
% end

if isempty(ObjCheck), error(['Unknown conventions: ' Obj.GLOBAL_SOFAConventions '. Can''t save.']); end

varNames = fieldnames(ObjCheck);
for ii=1:size(varNames,1);
    if ~isfield(Obj,varNames{ii})
        error(['Mandatory variable/attribute not existing: ' varNames{ii}]);
    end
end

%% Get & set dimensions
Obj=SOFAupdateDimensions(Obj);

%% Check convention: read-only variables
ObjCheck = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'r');
varNames = fieldnames(ObjCheck);

for ii=1:size(varNames,1);
    warningPostfix='';
	if ischar(Obj.(varNames{ii}))
		if ~strcmp(Obj.(varNames{ii}), ObjCheck.(varNames{ii}))
            if exist('OCTAVE_VERSION','builtin')
              % We're in Octave
    %                   if contains(varNames{ii},'GLOBAL_API') > 0; warningPostfix = ':API'; end
              if ismember(varNames{ii},{'GLOBAL_API'}); warningPostfix = ':API'; end
                  % In Octave 'contains' is not available, thus, the list has to be extended manually 
            else
              % We're in Matlab
              if contains(varNames{ii},'GLOBAL_API') > 0; warningPostfix = ':API'; end   
            end
			warning(['SOFA:save' warningPostfix],[varNames{ii} ' is read-only and was reset from ' Obj.(varNames{ii}) ' to '  ObjCheck.(varNames{ii})],0);
%             warning(['SOFA:save:' varNames{ii}], [varNames{ii} ' is read-only and was reset from ' Obj.(varNames{ii}) ' to '  ObjCheck.(varNames{ii})],0);
			Obj.(varNames{ii})=ObjCheck.(varNames{ii});
		end
	else
        if Obj.(varNames{ii}) ~= ObjCheck.(varNames{ii})
            if exist('OCTAVE_VERSION','builtin')
                % We're in Octave
            if ismember(varNames{ii},{'GLOBAL_API'}); warningPostfix = ':API'; end
                  % In Octave 'contains' is not available, thus, the list has to be extended manually 
            else
              % We're in Matlab
              if contains(varNames{ii},'GLOBAL_API') > 0; warningPostfix = ':API'; end   
            end
			warning(['SOFA:save' warningPostfix],[varNames{ii} ' is read-only and was reset from ' Obj.(varNames{ii}) ' to '  ObjCheck.(varNames{ii})],0);
%             warning(['SOFA:save:' varNames{ii}], [varNames{ii} ' is read-only and was reset from ' Obj.(varNames{ii}) ' to '  ObjCheck.(varNames{ii})],0);
            Obj.(varNames{ii})=ObjCheck.(varNames{ii});
        end
	end
end

%% check attributes (syntax, 1-dimensional string)
varNames = fieldnames(Obj);
for ii=1:size(varNames,1);
    
    if size(strfind(varNames{ii},'_'),2) == 1
        if ~ischar(Obj.(varNames{ii}))
            error(['Attribute not a valid string: ' varNames{ii} ' = ' num2str(Obj.(varNames{ii}))]);
        end
    elseif size(strfind(varNames{ii},'_'),2) > 1
        error(['Attribute not valid (only one underscore "_" is allowed in attribute name): ' varNames{ii}]);
    end
    
end

%% check varargin (compression)
if ~isempty(varargin) && isnumeric(varargin{1})
	if isscalar(varargin{1}) && varargin{1}>=0 && varargin{1}<=9
        Compression = varargin{1};
	else
        error('Error: Compression must be a numeric scalar value between 0 and 9.');
	end
else
    Compression = 1; % default
end

%% Set/modify time information
Obj.GLOBAL_DateModified=datestr(now,Def.dateFormat);
if isempty(Obj.GLOBAL_DateCreated), Obj.GLOBAL_DateCreated=Obj.GLOBAL_DateModified; end

%% Save file
NETCDFsave(filename,Obj,Compression);
