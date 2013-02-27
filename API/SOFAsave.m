function [] = SOFAsave(filename,Dataset,varargin)
%SOFASAVE 
%   [] = SOFAsave(filename,Dataset,Compression) creates a new SOFA file and
%   writes an entire data set to it.
%
%   filename specifies the name of the SOFA file to which the data is written.
%   Dataset is either a struct or a cell array containing the data and meta
%   data to be written to the SOFA file (see below for exact format).
%   Compression is an optional numeric value between 0 and 9 specifying the
%   amount of compression to be applied to the data when writing to the netCDF file.
%   0 is no compression and 9 is the most compression.
% 
%   If Dataset is a struct, it must contain one field called 'Data' for the data
%   and additional fields for each metadata value. The name of these fields
%   are identical to the names of the metadata.
%   If Dataset is a cell, it must have the following structure:
%   Dataset{x}{y}
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
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% --------------------- check and prepare variables ----------------------
% ------------------------- check input variables -------------------------
filename=SOFAcheckFilename(filename);

% -- check varargin --
if ~isempty(varargin) && isnumeric(varargin{1})
	if isscalar(varargin{1}) && varargin{1}>=0 && varargin{1}<=9
        Compression = varargin{1};
	else
        error('Error: Compression must be a numeric scalar value between 0 and 9.');
	end
else
    % setting a default compression
    % FIXME: (hagenw) what value should be the default compression?
    Compression = 9;
end
% -- check Dataset --
if isstruct(Dataset)
    VarNames = fieldnames(Dataset);
    numVars = size(VarNames,1);%... number of input variables
elseif iscell(Dataset)
    numVars = size(Dataset,2);%... number of input variables
end
Temp = struct; %... temporary variable to create a new struct from cell Dataset
ii = 1;
if iscell(Dataset) % all of this is only necessary if Dataset is a cell
	while ii<=numVars % loop through all input variables
        if mod(size(Dataset{ii},2),2)~=0 % variable name and value are not given correctly
            fprintf(2,'Error: Invalid Arguments.\n');
            return;
        end
        if ~strcmp(cellstr(class(Dataset{ii}{1})),'char') % check type of variable name
            %Dataset{ii}{3}{1}; % TODO???
            fprintf(2,'Error: Invalid Arguments (variable names must be strings).\n');
            return;
        end
        Temp = setfield(Temp,Dataset{ii}{1},Dataset{ii}{2});
        VarNames{ii} = Dataset{ii}{1}; % save variable names for mandatory check
        ii = ii + 1;
    end
    Dataset = Temp; % Dataset is now a struct
    clear Temp; % free memory
elseif ~isstruct(Dataset)
    error('Dataset must be struct or cell array.');
end

% --------------------------- prepare variables ---------------------------
mandatoryVars=SOFAgetVariables('mandatory'); %... mandatory variables
% -- check if mandatroy variables exist --
for ii=1:size(mandatoryVars,2) % go through all mandatory variables
    if ~sum(strcmp(mandatoryVars{ii}{1},VarNames)) % if there's no 1 anywhere
        if ~isempty(mandatoryVars{ii}{2}) % if a default value exists
            Dataset = setfield(Dataset,mandatoryVars{ii}{1},mandatoryVars{ii}{2}); % assign name of missing variable
            numVars = numVars + 1; % update number of input variables
        else % if no default value exists
            fprintf(2,['Error: Mandatory variable ' mandatoryVars{ii}{1} ' is missing.\n']);
            return;
        end
    end
end
Dataset.SOFAVersion = SOFAgetVersion(); % write SOFA version

%% --------------------------- N E T C D F save ---------------------------
NETCDFsave(filename,Dataset,Compression);

end %of function