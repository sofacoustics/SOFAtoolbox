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

%% --------------------------- N E T C D F save ---------------------------
% create file
ncid = netcdf(filename,'c','NETCDF4 with classical model');%'NETCDF4 with classical model'

% Define dimensions
% M - number of measurements
ncid('M') = Obj.M;
% R - number of receivers
ncid('R') = Obj.R;
% N - number of data samples per measurement
ncid('N') = Obj.N;
% E - number of emitters
ncid('E') = Obj.E;
% C - coordinate dimension
ncid('C') = 3;
% Q - quaternions (optional) dimension
% TODO: ?
% I - singleton dimension
ncid('I') = 1;


%% ===== Loop through all fields in Obj ==================================
attributes = fieldnames(Obj);
for ii=1:length(attributes)
    % store the current field name and its value
    attributeName = attributes{ii};
    attributeVal = Obj.(attributeName);


    % ----- DATA ----------------------------------------------------------
    if strcmp(attributeName,'Data')

        if strcmp(Obj.DataType,'FIR')
            % define dimensions
            ncid{'Data.FIR'} = ncdouble(Obj.N, ... % Samples
                                        Obj.R, ... % Receivers
                                        Obj.M);    % Measurements
            % store data
            ncid{'Data.FIR'}(:) = permute(attributeVal.FIR,[3 2 1]);

        elseif strcmp(Obj.DataType,'SpectraMagnitudePhase')
            % define dimensions
            ncid{'Data.Mag'} =   ncdouble(Obj.N, ... % Samples
                                          Obj.R, ... % Receivers
                                          Obj.M);    % Measurements
            ncid{'Data.Phase'} = ncdouble(Obj.N, ... % Samples
                                          Obj.R, ... % Receivers
                                          Obj.M);    % Measurements
            % store data
            ncid{'Data.Mag'}(:) = permute(attributeVal.Mag,[3 2 1]);
            ncid{'Data.Phase'}(:) = permute(attributeVal.Phase,[3 2 1]);
        end


    % ----- GLOBAL ATTRIBUTES --------------------------------------------
    elseif ~isempty(strfind(attributeName,'GLOBAL'))
        % store as global attribute
        ncid.(attributeName) = attributeVal;


    % ----- STRING VARIABLES ----------------------------------------------
    elseif ischar(attributeName)
        % define dimension
        [attributeName 'DIM']
        ncid([attributeName 'DIM']) = size(attributeVal,2);
        ncid{attributeName} = ncchar([attributeName 'DIM'],Obj.I);
        % store string
        ncid{attributeName}(:) = attributeVal;


    % ----- NUMERIC VARIABLES ---------------------------------------------
    elseif ndims(attributeVal)==2
        attributeVal = transpose(attributeVal);
        % define dimensions
        if size(attributeVal)==[3 1]                                     % [3 1]
            ncid{attributeName} = ncfloat(Obj.C,Obj.I);
        elseif size(attributeVal)==[3 Obj.M]                             % [3 M]
            ncid{attributeName} = ncfloat(Obj.C,Obj.M);
        % "normal" numeric variables
        elseif size(attributeVal)==[1 1]                                 % [1 1]
            ncid{attributeName} = ncfloat(Obj.I,Obj.I);
        elseif size(attributeVal)==[1 M]                                 % [M 1]
            ncid{attributeName} = ncfloat(Obj.I,Obj.M);
        elseif size(attributeVal,1)==1                                   % [1 x]
            ncid(attributeName) = size(attributeVal,1);
            ncid{attributeName} = ncfloat(attributeName,Obj.I);
        elseif size(attributeVal,1)==M                                   % [M x]
            ncid(attributeName) = size(attributeVal,1);
            ncid{attributeName} = ncfloat(attributeName,Obj.M);
        end
        % store variable
        ncid{atributeName}(:) = attributeVal;
    elseif ndims(attributeVal)==3
        attributeVal = permute(attributeVal,[3 2 1]);
        % define dimensions
        % receiver/transmitter positions
        if size(attributeVal)==[1 3]                                   % [1 3 1]
            ncid{attributeName} = ncfloat(Obj.I,Obj.C,Obj.I); 
        elseif size(attributeVal)==[R 3]                               % [1 3 R]
            ncid{attributeName} = ncfloat(Obj.R,Obj.C,Obj.I);
        elseif size(attributeVal)==[1 3 M]                             % [M 3 1]
            ncid{attributeName} = ncfloat(Obj.I,Obj.C,Obj.M);
        elseif size(attributeVal)==[R 3 M]                             % [M 3 R]
            ncid{attributeName} = ncfloat(Obj.R,Obj.C,Obj.M);
        end
        % store variable
        ncid{attributeName}(:) = attributeVal; 
    else
        error('%s: your variable type is not supported by SOFA',upper(mfilename));
    end
end

close(ncid);

if isunix
    unix(['nccopy -k 4 -d ' num2str(Compression) ' ' filename ' temp_netCDF3_to_netCDF4.nc']);
    unix(['cp temp_netCDF3_to_netCDF4.nc ' filename]);
    unix('rm temp_netCDF3_to_netCDF4.nc');
else
    dos(['nccopy -k 4 -d ' num2str(Compression) ' ' filename ' temp_netCDF3_to_netCDF4.nc']);
    dos(['cp temp_netCDF3_to_netCDF4.nc ' filename]);
    dos('rm temp_netCDF3_to_netCDF4.nc');
end

end %of function
