function [Obj] = SOFAsave(filename,Obj,varargin)
%SOFASAVE 
%   [] = SOFAsave(filename,Obj,Compression) creates a new SOFA file and
%   writes an entire data set to it.
%
%   filename specifies the name of the SOFA file to which the data is written.
%   Obj is a struct containing the data and meta
%   data to be written to the SOFA file (see below for exact format).
%   Compression is an optional numeric value between 0 and 9 specifying the
%   amount of compression to be applied to the data when writing to the netCDF file.
%   0 is no compression and 9 is the most compression.
% 
%   The existence of mandatory variables will be checked. The dimensions
%   will be updated.

% SOFA API - function SOFAsave
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% check file name
filename=SOFAcheckFilename(filename);

%% Check convention: mandatory variables
ObjCheck = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'m');

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
    if ischar(Obj.(varNames{ii}))
        if ~strcmp(Obj.(varNames{ii}), ObjCheck.(varNames{ii}))
            error(['Read-only variable/attribute was modified: ' varNames{ii}]);
        end
    else
        if Obj.(varNames{ii}) ~= ObjCheck.(varNames{ii})
            error(['Read-only variable/attribute was modified: ' varNames{ii}]);
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
Obj.GLOBAL_DatabaseTimeModified=datestr(now,31);
if isempty(Obj.GLOBAL_DatabaseTimeCreated), Obj.GLOBAL_DatabaseTimeCreated=Obj.GLOBAL_DatabaseTimeModified; end

%% Save file
NETCDFsave(filename,Obj,Compression);
