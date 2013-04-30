function NETCDFsave(filename,Obj,Compression)
%NETCDFSAVE
%   NETCDFsave(filename,Obj,Compression) saves all data and metadata to
%   a SOFA file.

% SOFA API - function octave/NETCDFsave
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 


%% ===== Check and prepare variables =====================================
%varNames = fieldnames(Dataset);
%numVars = size(varNames,1);
%dimNames = SOFAgetDimensions();
%SourceListenerVars=SOFAgetVariables('sourcelistener');
%TransmitterReceiverVars=SOFAgetVariables('transmitterreceiver');
%[numMeasurements,numReceivers,numSamples]=SOFAcheckDimensions(Dataset);

Def = SOFAdefinitions;
dims = Def.dimensions;

%% --------------------------- N E T C D F save ---------------------------
% create file
ncid = netcdf(filename,'c','NETCDF4 with classical model');

%% ===== Loop through all fields in Obj -1- ==============================
% remove Dimension definition fields before the loop
Obj1 = rmfield(Obj,'Dimensions');
fields = fieldnames(Obj1);
for ii=1:length(fields)
    % store current field name and its value
    fieldName = fields{ii};
    fieldVal = Obj.(fieldName);

    % ----- Dimensions ---------------------------------------------------
    if any(strcmp(struct2cell(dims),fieldName))
        % FIXME: this is not working at the moment, because we can also have
        % [R C M] for receiver positions, but M is only allowed as the first
        % dimension, when it should be able to grow dynamically.
        %if strcmp('M',fieldName)
        %    % store the M dimension as the one that can grow
        %    ncid(fieldName) = 0;
        %    ncid{fieldName} = ncdouble(fieldName);
        %    ncid{fieldName}(1:fieldVal) = 1:fieldVal;
        %else
            % create dimension
            ncid(fieldName) = fieldVal;
            ncid{fieldName} = ncdouble(fieldName);
            ncid{fieldName}(:) = 1:fieldVal;
        %end
    elseif strcmp('RoomCorner',fieldName)
        ncid('NumberOfRoomCorners') = 2;
        ncid{'NumberOfRoomCorners'} = ncdouble('NumberOfRoomCorners');
        ncid{'NumberOfRoomCorners'}(:) = 1:2;
    end
end


%% ===== Loop through all fields in Obj -2- ==============================
% remove dimension fields before the loop
Obj2 = rmfield(Obj1,struct2cell(dims));
fields = fieldnames(Obj2);
for ii=1:length(fields);
    % store the current field name and its value
    fieldName = fields{ii};
    fieldVal = Obj.(fieldName);


    % ----- DATA ----------------------------------------------------------
    if strcmp(fieldName,'Data')

        dataFields = fieldnames(fieldVal);
        for jj=1:length(dataFields)
            % store current Data field name and its value
            dataFieldName = dataFields{jj};
            dataFieldVal = fieldVal.(dataFieldName);
            
            if isempty(strfind(dataFieldName,'_')) % skip attributes

                % FIXME: find a better way than the checking of the length!
                % Should the data matrix be allowed to have every number of
                % dimension they want?
                if ndims(dataFieldVal)==2 && ...
                    size(dataFieldVal)==[1 1]                       % [I]
                    ncid{['Data.' dataFieldName]} = ncdouble(dims.I);
                elseif ndims(dataFieldVal)==2 && ...
                    strcmp('TOAModel',dataFieldName)                % [M 5]
                    ncid('TOAModelParameter') = 5;
                    ncid{['Data.' dataFieldName]} = ncdouble(dims.M, ...
                        'TOAModelParameter');
                elseif ndims(dataFieldVal)==2 && ...
                    size(dataFieldVal)==[Obj.M Obj.R]               % [M R]
                    ncid{['Data.' dataFieldName]} = ncdouble(dims.M,dims.R);
                elseif ndims(dataFieldVal)==3 && ...
                    size(dataFieldVal)==[Obj.M Obj.R Obj.N]         % [M R N]
                    ncid{['Data.' dataFieldName]} = ncdouble(dims.M, ...
                                                             dims.R, ...
                                                             dims.N);
                end
                % store data
                ncid{['Data.' dataFieldName]}(:) = dataFieldVal;
            end
        end
        % store attributes of the variables
        for jj=1:length(dataFields)
            if ~isempty(strfind(dataFields{jj},'_')) % only attributes
                % store current Field field name and its value
                dataFieldName = dataFields{jj};
                dataFieldVal = fieldVal.(dataFieldName);
                dataFieldNameBase = ...
                    dataFieldName(1:strfind(dataFieldName,'_')-1);
                dataFieldNameAttr = ...
                    dataFieldName(strfind(dataFieldName,'_')+1:end);
                % store attribute
                ncid{['Data.' dataFieldNameBase]}.(dataFieldNameAttr) = dataFieldVal;
            end
        end


    % ----- GLOBAL ATTRIBUTES --------------------------------------------
    elseif ~isempty(strfind(fieldName,'GLOBAL'))
        % strip "GLOBAL_" from the name
        fieldName = fieldName(strfind(fieldName,'_')+1:end);
        % store as global attribute
        ncid.(fieldName) = fieldVal;


    % ----- LISTENER and SOURCE -------------------------------------------
    elseif (~isempty(strfind(fieldName,'Listener')) || ...
           ~isempty(strfind(fieldName,'Source')) ) && ...
            isempty(strfind(fieldName,'_'))

        if size(fieldVal)==[Obj.I Obj.C]                            % [C]
            ncid{fieldName} = ncfloat(dims.C);
        elseif size(fieldVal)==[Obj.M Obj.C]                        % [M C]
            ncid{fieldName} = ncfloat(dims.M,dims.C);
        end

        % store variable
        ncid{fieldName}(:) = fieldVal;


    % ----- RECEIVER -----------------------------------------------------
    elseif ~isempty(strfind(fieldName,'Receiver')) && ...
            isempty(strfind(fieldName,'_'))

        if ndims(fieldVal)==2 && size(fieldVal)==[Obj.R Obj.C]      % [R C]
            ncid{fieldName} = ncfloat(dims.R,dims.C);
        elseif ndims(fieldVal)==3 && size(fieldVal)==[Obj.R Obj.C Obj.M] % [R C M]
            ncid{fieldName} = ncfloat(dims.R,dims.C,dims.M);
        end

        % store variable
        ncid{fieldName}(:) = fieldVal;


    % ----- EMITTER ------------------------------------------------------
    elseif ~isempty(strfind(fieldName,'Emitter')) && ...
            isempty(strfind(fieldName,'_'))

        if ndims(fieldVal)==2 && size(fieldVal)==[Obj.E Obj.C]      % [E C]
            ncid{fieldName} = ncfloat(dims.E,dims.C);
        elseif ndims(fieldVal)==3 && size(fieldVal)==[Obj.E Obj.C Obj.M] % [E C M]
            ncid{fieldName} = ncfloat(dims.E,dims.C,dims.M);
        end

        % store variable
        ncid{fieldName}(:) = fieldVal;


    % ----- ROOM CORNERS --------------------------------------------------
    elseif ~isempty(strfind(fieldName,'RoomCorner')) && ...
            isempty(strfind(fieldName,'_'))

        if ndims(fieldVal)==2 && size(fieldVal)==[2 Obj.C]          % [2 C]
            ncid{fieldName} = ncfloat('NumberOfRoomCorners',dims.C);
        elseif ndims(fieldVal)==3 && size(fieldVal)==[2 Obj.C Obj.M] % [2 C M]
            ncid{fieldName} = ncfloat('NumberOfRoomCorners',dims.C,dims.M);
        end

        % store variable
        ncid{fieldName}(:) = fieldVal;


    % ----- NUMERIC VARIABLES ---------------------------------------------
    % store other numeric vectors
    elseif isnumeric(fieldVal)

        if ndims(fieldVal)==2
            ncid{fieldName} = ncfloat([fieldName 'Dimension1'], ...
                [fieldName 'Dimension2']);
        elseif ndims(fieldVal)==3
            ncid{fieldName} = ncfloat([fieldName 'Dimension1'], ...
                [fieldName 'Dimension2'], ...
                [fieldName 'Dimension3']);
        end

        % store variable
        ncid{fieldName}(:) = fieldVal;

        
    % ----- ATTRIBUTES ---------------------------------------------------
    % here we store descriptive attributes of variables, for example units or
    % long names
    elseif ~isempty(strfind(fieldName,'_'))
        % split "VariableName_Attribute" into "VariableName" and "Attribute"
        fieldNameBase = ...
            fieldName(1:strfind(fieldName,'_')-1);
        fieldNameAttr = ...
            fieldName(strfind(fieldName,'_')+1:end);
        % store field
        ncid{fieldNameBase}.(fieldNameAttr) = fieldVal;
    
    else
        error('%s: your variable type is not supported by SOFA',upper(mfilename));
    end
end

close(ncid);

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
