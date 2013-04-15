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
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% Global definitions
dims={'i';'r';'e';'n';'m';'c';'q'}; % dimensions

%% check file name
filename=SOFAcheckFilename(filename);

%% Check convention: mandatory variables
ObjCheck = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'m');
Obj.Dimensions=ObjCheck.Dimensions;

varNames = fieldnames(ObjCheck);
for ii=1:size(varNames,1);
    if ~isfield(Obj,varNames{ii})
        error(['Mandatory variable/attribute not existing: ' varNames{ii}]);
    end
end

%% Get & set dimensions
varNames = fieldnames(Obj.Dimensions);
dataCount=0; % dimensions in Data
for ii=1:size(varNames,1)
    for jj=1:size(Obj.Dimensions.(varNames{ii}),2)
        if strcmp(varNames{ii},'Data')
        % Data
            dataCount=dataCount+1;
            varNamesData = fieldnames(Obj.Dimensions.Data);
            for ll=1:size(Obj.Dimensions.Data.(varNamesData{dataCount}),2)
                for kk=1:size(dims,1)
                    if strfind(Obj.Dimensions.Data.(varNamesData{dataCount})(ll),dims{kk}) % find lowercase dimension letter
                        Obj.(upper(dims{kk}))=size(Obj.Data.(varNamesData{dataCount}),ll); % save global dimension variable
                        Obj.Dimensions.Data.(varNamesData{dataCount})(ll)=upper(Obj.Dimensions.Data.(varNamesData{dataCount})(ll)); % LC -> UC
                    end
                end
            end
        else
        % DataVar
%         disp(ii);
        
        
            for ll=1:size(Obj.Dimensions.(varNames{ii}){jj},2)
                for kk=1:size(dims,1)
                    if strfind(Obj.Dimensions.(varNames{ii}){jj}(ll),dims{kk}) % find lowercase dimension letter
                        Obj.(upper(dims{kk}))=size(Obj.(varNames{ii}),ll); % save global dimension variable
                        Obj.Dimensions.(varNames{ii}){jj}(ll)=upper(Obj.Dimensions.(varNames{ii}){jj}(ll)); % LC -> UC
                    end
                end
            end
            
            
            
        end
    end
    % set dimensions
    for jj=1:size(Obj.Dimensions.(varNames{ii}),2) % set dimension string (instead of array)
        if iscellstr(Obj.Dimensions.(varNames{ii})) 
            if length(size(Obj.(varNames{ii})))==2 % 2-dimensions
                if size(Obj.(varNames{ii}))==[Obj.(Obj.Dimensions.(varNames{ii}){jj}(1)),Obj.(Obj.Dimensions.(varNames{ii}){jj}(2))]
%                     disp(varNames{ii});
                    Obj.Dimensions.(varNames{ii})=[Obj.Dimensions.(varNames{ii}){jj}(1) Obj.Dimensions.(varNames{ii}){jj}(2)];
                end
            elseif length(size(Obj.(varNames{ii})))==3  % 3-dimensions
                if size(Obj.(varNames{ii}))==[Obj.(Obj.Dimensions.(varNames{ii}){jj}(1)),Obj.(Obj.Dimensions.(varNames{ii}){jj}(2)),Obj.(Obj.Dimensions.(varNames{ii}){jj}(3))]
%                     disp(varNames{ii});
                    Obj.Dimensions.(varNames{ii})=[Obj.Dimensions.(varNames{ii}){jj}(1) Obj.Dimensions.(varNames{ii}){jj}(2) Obj.Dimensions.(varNames{ii}){jj}(3)];
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
    Compression = 1;
end

%% Set/modify time information

Obj.GLOBAL_DatabaseTimeModified=datestr(now,31);
if isempty(Obj.GLOBAL_DatabaseTimeCreated), Obj.GLOBAL_DatabaseTimeCreated=Obj.GLOBAL_DatabaseTimeModified; end

%% Save file
NETCDFsave(filename,Obj,Compression);
