function NETCDFsave(Filename,Dataset,Compression)
%NETCDFSAVE
%   NETCDFsave(Filename,Dataset,Compression) saves all data and metadata to
%   a SOFA file.

% SOFA API - function SOFAloadData
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
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
ncid = netcdf(Filename,'c','64bit-offset');

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

  % ----- STRING VARIABLES -----------------------------------------------
  if ischar(currentVarValue) % -- if currentVarValue is a string
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

  % ----- DATA MATRIX ----------------------------------------------------
  elseif strcmp(currentVarName,'Data')  % data [measurements receiver samples]
    numMeasurements = size(currentVarValue,1);
    numReceivers = size(currentVarValue,2);
    ncid(dimNames.Measurements) = numMeasurements;
    ncid(dimNames.Receivers) = numReceivers;
    if strcmp(Dataset.DataType,'FIR')
      numSamples = length(currentVarValue(1,1).FIR);
      ncid(dimNames.Samples) = numSamples;
      ncid{'Data.FIR'} = ncdouble(dimNames.Measurements, ...
                                dimNames.Receivers, ...
                                dimNames.Samples);
      ncid{'Data.FIR'}(:) = currentVarValue.FIR;
    elseif strcmp(Dataset.dataType,'SpectraMagnitudePhase')
      numSamples = length(currentVarValue(1,1).Mag);
      ncid(dimNames.Samples) = numSamples;
      ncid{'Data.Mag'} = ncdouble(dimNames.Measurements, ...
                                dimNames.Receivers, ...
                                dimNames.Samples);
      ncid{'Data.Phase'} = ncdouble(dimNames.Measurements, ...
                                  dimNames.Receivers, ...
                                  dimNames.Samples);

      ncid{'Data.Mag'}(:) = currentVarValue.Mag;
      ncid{'Data.Phase'}(:) = currentVarValue.Phase;
    end

  % ----- NUMERIC VARIABLES ----------------------------------------------
  elseif ndims(currentVarValue)==2
    % positions and vectors
    if size(currentVarValue)==[1 3]                                     % [1 3]
      ncid{currentVarName} = ncfloat(dimNames.Scalar,dimNames.Coordinates);
    elseif size(currentVarValue)==[numMeasurements 3]                                 % [numMeasurements 3]
      ncid{currentVarName} = ncfloat(dimNames.Measurements,dimNames.Coordinates);
    % "normal" numeric variables
    elseif size(currentVarValue)==[1 1]                                 % [1 1]
      ncid{currentVarName} = ncfloat(dimNames.Scalar,dimNames.Scalar);
    elseif size(currentVarValue)==[numMeasurements 1]                                 % [numMeasurements 1]
      ncid{currentVarName} = ncfloat(dimNames.Measurements,dimNames.Scalar);
    elseif size(currentVarValue,1)==1                                   % [1 x]
      ncid(currentVarName) = size(currentVarValue,2);
      ncid{currentVarName} = ncfloat(dimNames.Scalar,currentVarName);
    elseif size(currentVarValue,1)==numMeasurements                                   % [numMeasurements x]
      ncid(currentVarName) = size(currentVarValue,2);
      ncid{currentVarName} = ncfloat(dimNames.Measurements,currentVarName);
    end
    % store variable
    ncid{currentVarName}(:) = currentVarValue;
  elseif ndims(currentVarValue)==3
    % receiver/transmitter positions
    if size(currentVarValue)==[1 3 1]                                   % [1 3 1]
      ncid{currentVarName} = ...
        ncfloat(dimNames.Scalar,dimNames.Coordinates,dimNames.Scalar);
    elseif size(currentVarValue)==[1 3 numReceivers]                               % [1 3 numReceivers]
      ncid{currentVarName} = ...
        ncfloat(dimNames.Scalar,dimNames.Coordinates,dimNames.Receivers);
    elseif size(currentVarValue)==[numMeasurements 3 1]                               % [numMeasurements 3 1]
      ncid{currentVarName} = ...
        ncfloat(dimNames.Measurements,dimNames.Coordinates,dimNames.Scalar);
    elseif size(currentVarValue)==[numMeasurements 3 numReceivers]                               % [numMeasurements 3 numReceivers]
      ncid{currentVarName} = ...
        ncfloat(dimNames.Measurements,dimNames.Coordinateso,dimNames.Receivers);
    end
    % store variable
    ncid{currentVarName}(:) = currentVarValue;
  
  
  else
    error('%s: your variable type is not supported by SOFA',upper(mfilename));
  end

end

close(ncid);

end %of function