function [] = SOFAsave(filename,Obj,Var,DataVar,varargin)
%SOFASAVE 
%   [] = SOFAsave(filename,Obj,Compression) creates a new SOFA file and
%   writes an entire data set to it.
%
%   filename specifies the name of the SOFA file to which the data is written.
%   Obj is either a struct or a cell array containing the data and meta
%   data to be written to the SOFA file (see below for exact format).
%   Compression is an optional numeric value between 0 and 9 specifying the
%   amount of compression to be applied to the data when writing to the netCDF file.
%   0 is no compression and 9 is the most compression.
% 
%   If Obj is a struct, it must contain one field called 'Data' for the data
%   and additional fields for each metadata value. The name of these fields
%   are identical to the names of the metadata.
%   If Obj is a cell, it must have the following structure:
%   Obj{x}{y}
%   x ... number of variable
%   y = 1: variable name; y = 2: value
% 
%   In both cases, the existence of mandatory variables will be checked.
%   Coordinate variables are expected to have one of the following
%   dimensions (with numMeasurements being the number of measurements):
%   Source/ListenerPosition, -View, -Up: [1 3], [numMeasurements 3]
%   Transmitter/ReceiverPosition: [1 3 1], [1 3 numReceivers], [numMeasurements 3 1], [numMeasurements 3 numReceivers]
%   (with numReceivers being the number of receivers or transmitters respectively)
%
%   All other meta data variables must have one of the following dimensions:
%   [1 1], [1 x], [numMeasurements 1], [numMeasurements x] (x is arbitary)

% SOFA API - function SOFAsave
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% Last Update: Michael Mihocic, 09.04.2013

% --------------------- check and prepare variables ----------------------
%% check input variables
filename=SOFAcheckFilename(filename);

%% Check convention: mandatory variables
ObjCheck = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'m');
varNames = fieldnames(ObjCheck);

for ii=1:size(varNames,1);
    % check if variable/attribute is existing
    if ~isfield(Obj,varNames{ii})
        error(['Mandatory variable/attribute not existing: ' varNames{ii}]);
    end
end

%% Check convention: read-only variables
ObjCheck = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'rm');
varNames = fieldnames(ObjCheck);

for ii=1:size(varNames,1);
    % check if variable/attribute is existing
%     if Obj.(varNames{ii}) ~= ObjCheck.(varNames{ii})
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
%         disp(['attribute: ' varNames{ii}]);
        if ~ischar(Obj.(varNames{ii}))
            error(['Attribute no valid string: ' varNames{ii} ' = ' num2str(Obj.(varNames{ii}))]);
        end
    elseif size(strfind(varNames{ii},'_'),2) > 1
        error(['Attribute not valid (only one underscore "_" is allowed in attribute name): ' varNames{ii}]);
    else 
%         disp(['variable: ' varNames{ii}]);
    end
    
end

%% check varargin 
if ~isempty(varargin) && isnumeric(varargin{1})
	if isscalar(varargin{1}) && varargin{1}>=0 && varargin{1}<=9
        Compression = varargin{1};
	else
        error('Error: Compression must be a numeric scalar value between 0 and 9.');
	end
else
    % default
    Compression = 1;
end

%% Save file
NETCDFsave(filename,Obj,Var,DataVar,Compression);

end %of function