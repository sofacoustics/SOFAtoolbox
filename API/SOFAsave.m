function [] = SOFAsave(filename,Obj,varargin)
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

%% --------------------- check and prepare variables ----------------------

%% check file name
filename=SOFAcheckFilename(filename);
% if ~ischar(filename)
% 	error('Filename must be a string.');
% end
% idx=strfind(filename,'.');
% if isempty(idx)
%     filename=[filename '.sofa'];
% elseif ~strcmp(filename(idx+1:end),'sofa')
%     error(['SOFA-API does not support *' filename(idx:end) '-files!'])
% end

%% Check convention: mandatory variables
[ObjCheck,Var,DataVar] = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'m');

%%%%%%%%%%%%% LÖSCHEN %%%%%%%%%%%%%%%
Obj.DIMENSIONS=Var;
Obj.DIMENSIONS.Data.IR=DataVar.IR;
Obj.DIMENSIONS.Data.SamplingRate=DataVar.SamplingRate;
clear Var DataVar;
%%%%%%%%%%% ENDE LÖSCHEN %%%%%%%%%%%%

varNames = fieldnames(ObjCheck);
for ii=1:size(varNames,1);
    % check if variable/attribute is existing
    if ~isfield(Obj,varNames{ii})
        error(['Mandatory variable/attribute not existing: ' varNames{ii}]);
    end
end

%% Get & set dimensions
dims={'u';'r';'e';'n';'m';'c';'q'}; % dimensions
varNames = fieldnames(Obj.DIMENSIONS);
dataCount=0; % dimensions in Data
for ii=1:size(varNames,1)
    for jj=1:size(Obj.DIMENSIONS.(varNames{ii}),2)
        if strcmp(varNames{ii},'Data')
        % Data
            dataCount=dataCount+1;
            varNamesData = fieldnames(Obj.DIMENSIONS.Data);
            for ll=1:size(Obj.DIMENSIONS.Data.(varNamesData{dataCount}),2)
                for kk=1:size(dims,1)
                    if strfind(Obj.DIMENSIONS.Data.(varNamesData{dataCount})(ll),dims{kk}) % find lowercase dimension letter
                        Obj.(upper(dims{kk}))=size(Obj.Data.(varNamesData{dataCount}),ll); % save global dimension variable
                        Obj.DIMENSIONS.Data.(varNamesData{dataCount})(ll)=upper(Obj.DIMENSIONS.Data.(varNamesData{dataCount})(ll)); % LC -> UC
                    end
                end
            end
        else
        % DataVar
            for ll=1:size(Obj.DIMENSIONS.(varNames{ii}){jj},2)
                for kk=1:size(dims,1)
                    if strfind(Obj.DIMENSIONS.(varNames{ii}){jj}(ll),dims{kk}) % find lowercase dimension letter
                        Obj.(upper(dims{kk}))=size(Obj.(varNames{ii}),ll); % save global dimension variable
                        Obj.DIMENSIONS.(varNames{ii}){jj}(ll)=upper(Obj.DIMENSIONS.(varNames{ii}){jj}(ll)); % LC -> UC
                    end
                end
            end
        end
    end
    % set dimensions
    for jj=1:size(Obj.DIMENSIONS.(varNames{ii}),2) % set dimension string (instead of array)
        if iscellstr(Obj.DIMENSIONS.(varNames{ii})) 
            if length(size(Obj.(varNames{ii})))==2 % 2-dimensions
                if size(Obj.(varNames{ii}))==[Obj.(Obj.DIMENSIONS.(varNames{ii}){jj}(1)),Obj.(Obj.DIMENSIONS.(varNames{ii}){jj}(2))]
                    disp(varNames{ii});
                    Obj.DIMENSIONS.(varNames{ii})=[Obj.DIMENSIONS.(varNames{ii}){jj}(1) Obj.DIMENSIONS.(varNames{ii}){jj}(2)];
                end
            elseif length(size(Obj.(varNames{ii})))==3  % 3-dimensions
                if size(Obj.(varNames{ii}))==[Obj.(Obj.DIMENSIONS.(varNames{ii}){jj}(1)),Obj.(Obj.DIMENSIONS.(varNames{ii}){jj}(2)),Obj.(Obj.DIMENSIONS.(varNames{ii}){jj}(3))]
                    disp(varNames{ii});
                    Obj.DIMENSIONS.(varNames{ii})=[Obj.DIMENSIONS.(varNames{ii}){jj}(1) Obj.DIMENSIONS.(varNames{ii}){jj}(2) Obj.DIMENSIONS.(varNames{ii}){jj}(3)];
                end
            end
        end
    end
end

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
NETCDFsave(filename,Obj,Compression);

end %of function