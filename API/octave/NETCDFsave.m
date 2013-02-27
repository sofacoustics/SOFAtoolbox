function NETCDFsave(filename,Dataset,Compression)
%NETCDFSAVE
%   NETCDFsave(filename,Dataset,Compression) saves all data and metadata to
%   a SOFA file.

% SOFA API - function octave/NETCDFsave
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% --------------------- check and prepare variables ----------------------
varNames = fieldnames(Dataset);
numVars = size(varNames,1);
dimNames = SOFAgetDimensions();
SourceListenerVars=SOFAgetVariables('sourcelistener');
TransmitterReceiverVars=SOFAgetVariables('transmitterreceiver');
[numMeasurements,numReceivers,numSamples]=SOFAcheckDimensions(Dataset);

%% --------------------------- N E T C D F save ---------------------------
% create file
ncid = netcdf(filename,'c','NETCDF4');

% define some fixed dimension sizes
ncid(dimNames.Scalar) = 1;
ncid(dimNames.Coordinates) = 3;
% FIXME: this is not working at themoment
%ncid(dimNames.String) = 0; % unlimited length
ncid(dimNames.Measurements)=numMeasurements;
ncid(dimNames.Receivers)=numReceivers;
ncid(dimNames.Samples)=numSamples;



for ii=1:numVars % loop through all input variables
    currentVarName = varNames{ii};
    currentVarValue = getfield(Dataset,varNames{ii});

    % ----- STRING VARIABLES ----------------------------------------------
    if ischar(currentVarValue) % -- if currentVarValue is a string
        currentVarValue=transpose(currentVarValue);
        if size(currentVarValue,1)==1
            ncid(currentVarName) = size(currentVarValue,2);      % define dimension
            ncid{currentVarName} = ncchar(currentVarName);  % allocate dim
        else
            ncid(['x' currentVarName]) = size(currentVarValue,1);
            ncid(['y' currentVarName]) = size(currentVarValue,2);
            ncid{currentVarName} = ncchar(['x' currentVarName], ...
                ['y' currentVarName]);
        end
        ncid{currentVarName}(:) = currentVarValue;                     % store variable

    % ----- DATA MATRIX ---------------------------------------------------
    elseif strcmp(currentVarName,'Data')  % data [measurements receiver samples]
        if strcmp(Dataset.DataType,'FIR')
            ncid{'Data.FIR'} = ncdouble(dimNames.Samples, ...
                                        dimNames.Receivers, ...
                                        dimNames.Measurements);
            ncid{'Data.FIR'}(:) = permute(currentVarValue.FIR,[3 2 1]);
        elseif strcmp(Dataset.dataType,'SpectraMagnitudePhase')
            ncid{'Data.Mag'} = ncdouble(dimNames.Samples, ...
                                        dimNames.Receivers, ...
                                        dimNames.Measurements);
            ncid{'Data.Phase'} = ncdouble(dimNames.Samples, ...
                                        dimNames.Receivers, ...
                                        dimNames.Measurements);

            ncid{'Data.Mag'}(:) = permute(currentVarValue.Mag,[3 2 1]);
            ncid{'Data.Phase'}(:) = permute(currentVarValue.Phase,[3 2 1]);
        end

    % ----- NUMERIC VARIABLES ---------------------------------------------
    elseif ndims(currentVarValue)==2
        currentVarValue=transpose(currentVarValue);
        % positions and vectors
        if size(currentVarValue)==[3 1]                                     % [1 3]
            ncid{currentVarName} = ncfloat(dimNames.Coordinates,dimNames.Scalar);
        elseif size(currentVarValue)==[3 numMeasurements]                                 % [numMeasurements 3]
            ncid{currentVarName} = ncfloat(dimNames.Coordinates,dimNames.Measurements);
        % "normal" numeric variables
        elseif size(currentVarValue)==[1 1]                                 % [1 1]
            ncid{currentVarName} = ncfloat(dimNames.Scalar,dimNames.Scalar);
        elseif size(currentVarValue)==[1 numMeasurements]                                 % [numMeasurements 1]
            ncid{currentVarName} = ncfloat(dimNames.Scalar,dimNames.Measurements);
        elseif size(currentVarValue,1)==1                                   % [1 x]
            ncid(currentVarName) = size(currentVarValue,1);
            ncid{currentVarName} = ncfloat(currentVarName,dimNames.Scalar);
        elseif size(currentVarValue,1)==numMeasurements                                   % [numMeasurements x]
            ncid(currentVarName) = size(currentVarValue,1);
            ncid{currentVarName} = ncfloat(currentVarName,dimNames.Measurements);
        end
        % store variable
        ncid{currentVarName}(:) = currentVarValue;
    elseif ndims(currentVarValue)==3
        currentVarValue=permute(currentVarValue,[3 2 1]);
        % receiver/transmitter positions
        if size(currentVarValue)==[1 3]                                   % [1 3 1]
            ncid{currentVarName} = ...
                ncfloat(dimNames.Scalar,dimNames.Coordinates,dimNames.Scalar);
        elseif size(currentVarValue)==[numReceivers 3]                               % [1 3 numReceivers]
            ncid{currentVarName} = ...
                ncfloat(dimNames.Receivers,dimNames.Coordinates,dimNames.Scalar);
        elseif size(currentVarValue)==[1 3 numMeasurements]                               % [numMeasurements 3 1]
            ncid{currentVarName} = ...
                ncfloat(dimNames.Scalar,dimNames.Coordinates,dimNames.Measurements);
        elseif size(currentVarValue)==[numReceivers 3 numMeasurements]                               % [numMeasurements 3 numReceivers]
            ncid{currentVarName} = ...
                ncfloat(dimNames.Receivers,dimNames.Coordinateso,dimNames.Measurements);
        end
        % store variable
        ncid{currentVarName}(:) = currentVarValue; 
    else
        error('%s: your variable type is not supported by SOFA',upper(mfilename));
    end
end

close(ncid);

end %of function