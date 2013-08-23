function NETCDFsave(filename,Obj,Compression)
%NETCDFSAVE
%   NETCDFsave(filename,Obj,Compression) saves all data and metadata to
%   a SOFA file.

% SOFA API - function octave/NETCDFsave
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% Hagen Wierstorf


% --------------------------- N E T C D F save ---------------------------
try
    % open file
    ncid = netcdf(filename,'c','NETCDF4 with classical model');

    % loop through fileds in Obj
    fields = fieldnames(Obj);
    for ii=1:length(fields)
        % store current field name and its value
        fieldName = fields{ii};
        fieldVal = Obj.(fieldName);

        % ----- GLOBAL ATTRIBUTES ----------------------------------------
        if ~isempty(strfind(fieldName,'GLOBAL'))
            % strip "GLOBAL_" from the name
            fieldName = fieldName(strfind(fieldName,'_')+1:end);
            % store as global attribute
            ncid.(fieldName) = fieldVal;
        end
    end


    % ----- DIMENSIONS ---------------------------------------------------
    % get dimensions
    dims=cell2mat(fieldnames(rmfield(Obj.API,'Dimensions'))');
    for ii=1:length(dims)
        dimName = dims(ii);
        dimVal = Obj.API.(dimName);
        % create dimension
        ncid(dimName) = dimVal;
    end


    % ----- VARIABLES and ATTRIBUTES -------------------------------------
    fields = fieldnames(rmfield(Obj,{'Data','API'}));
    for ii=1:length(fields)
        % store current field name and its dimension
        fieldName = fields{ii};
        fieldVal = Obj.(fieldName);
        % --- Variables ---
        if isempty(strfind(fieldName,'_'))
            fieldDim = Obj.API.Dimensions.(fieldName);
            if length(fieldDim)==3
                ncid{fieldName} = ncfloat(fieldDim(1), ...
                                          fieldDim(2), ...
                                          fieldDim(3));
            elseif length(fieldDim)==2
                ncid{fieldName} = ncfloat(fieldDim(1), ...
                                          fieldDim(2));
            else
                ncid{fieldName} = ncfloat(fieldDim);
            end
            % store variable
            ncid{fieldName}(:) = fieldVal;
        % --- Attributes ---
        elseif isempty(strfind(fieldName,'GLOBAL_'))
            % here we store descriptive attributes of variables, for example
            % units or long names
            % split "VariableName_Attribute" into "VariableName" and "Attribute"
            fieldNameBase = fieldName(1:strfind(fieldName,'_')-1);
            fieldNameAttr = fieldName(strfind(fieldName,'_')+1:end);
            % store field
            ncid{fieldNameBase}.(fieldNameAttr) = fieldVal;
        end
    end

    % ----- DATA ---------------------------------------------------------
    fields = fieldnames(Obj.Data);
    for ii=1:length(fields)
        fieldName = fields{ii};
        fieldVal = Obj.Data.(fieldName);
        if isempty(strfind(fieldName,'_')) % skip all attributes
            fieldDim = Obj.API.Dimensions.Data.(fieldName);
            if length(fieldDim)==3
                ncid{['Data.' fieldName]} = ncfloat(fieldDim(1), ...
                                                    fieldDim(2), ...
                                                    fieldDim(3));
            elseif length(fieldDim)==2
                ncid{['Data.' fieldName]} = ncfloat(fieldDim(1), ...
                                                    fieldDim(2));
            else
                ncid{['Data.' fieldName]} = ncfloat(fieldDim);
            end
            % store variable
            ncid{['Data.' fieldName]}(:) = fieldVal;
        else
            % split "VariableName_Attribute" into "VariableName" and "Attribute"
            fieldNameBase = fieldName(1:strfind(fieldName,'_')-1);
            fieldNameAttr = fieldName(strfind(fieldName,'_')+1:end);
            % store field
            ncid{['Data.' fieldNameBase]}.(fieldNameAttr) = fieldVal;
        end
    end
catch
    if exist('ncid','var') && ~isempty(ncid)
        close(ncid);
    end
    error(['An error occured during writing the SOFA file: ' ...
        lasterror.message lasterror.stack]);
end_try_catch

% close file
close(ncid);

% Move files to get the desired compression and a newer NetCDF version
% FIXME: this doesn't work with old versions of nccopy
%if isunix
%    unix(['nccopy -k 4 -d ' num2str(Compression) ' ' filename ' temp_netCDF3_to_netCDF4.nc']);
%    unix(['cp temp_netCDF3_to_netCDF4.nc ' filename]);
%    unix('rm temp_netCDF3_to_netCDF4.nc');
%else
%    dos(['nccopy -k 4 -d ' num2str(Compression) ' ' filename ' temp_netCDF3_to_netCDF4.nc']);
%    dos(['cp temp_netCDF3_to_netCDF4.nc ' filename]);
%    dos('rm temp_netCDF3_to_netCDF4.nc');
%end

end %of function
