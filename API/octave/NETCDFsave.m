function NETCDFsave(Filename,Dataset,Compression)
% NETCDFsave savesthe given Dataset to the Filename using netCDF
%
%   Usage: NETCDFsave(Filename,Dataset,Compression)
%

VarNames = fieldnames(Dataset); % update VarNames
numVar = size(VarNames,1); % update number of variables


% get names for the dimensions, this is put in an extra function in order to
% change the names consistently at one place
DimNames = SOFAdimensions();


%% ----- N E T C D F save ------------------------------------------------

% create file
nc = netcdf(Filename,'c','64bit-offset');
%nc = netcdf(Filename,'c');

% define some fixed dimension sizes
nc(DimNames.Scalar) = 1;
nc(DimNames.Coordinates) = 3;
% FIXME: this is not working at themoment
%nc(DimNames.String) = 0; % unlimited length

% loop through all input variables
for ii=1:numVar
  
  CurrentVarName = VarNames{ii};
  CurrentVar = getfield(Dataset,VarNames{ii});

  % ----- STRING VARIABLES -----------------------------------------------
  if ischar(CurrentVar) % -- if CurrentVar is a string
    if size(CurrentVar,1)==1
      nc(CurrentVarName) = size(CurrentVar,2);      % define dimension
      nc{CurrentVarName} = ncchar(CurrentVarName);  % allocate dim
    else
      nc(['x' CurrentVarName]) = size(CurrentVar,1);
      nc(['y' CurrentVarName]) = size(CurrentVar,2);
      nc{CurrentVarName} = ncchar(['x' CurrentVarName], ...
        ['y' CurrentVarName]);
    end
    nc{CurrentVarName}(:) = CurrentVar;                     % store variable

  % ----- DATA MATRIX ----------------------------------------------------
  elseif strcmp(CurrentVarName,'Data')  % data [measurements receiver samples]
    M = size(CurrentVar,1);
    R = size(CurrentVar,2);
    nc(DimNames.Measurements) = M;
    nc(DimNames.Receivers) = R;
    if strcmp(Dataset.DataType,'FIR')
      N = length(CurrentVar(1,1).FIR);
      nc(DimNames.Samples) = N;
      nc{'Data.FIR'} = ncdouble(DimNames.Measurements, ...
                                DimNames.Receivers, ...
                                DimNames.Samples);
      Temp = zeros(M,R,N); % preallocating memory
      for m=1:M
        for r=1:R
          Temp(m,r,:) = CurrentVar(m,r).FIR;
        end
      end
      nc{'Data.FIR'}(:) = Temp;
      clear Temp;
    elseif strcmp(Dataset.dataType,'SpectraMagnitudePhase')
      N = length(CurrentVar(1,1).Mag);
      nc(DimNames.Samples) = N;
      nc{'Data.Mag'} = ncdouble(DimNames.Measurements, ...
                                DimNames.Receivers, ...
                                DimNames.Samples);
      nc{'Data.Phase'} = ncdouble(DimNames.Measurements, ...
                                  DimNames.Receivers, ...
                                  DimNames.Samples);
      Temp1 = zeros(M,R,N); % preallocating memory
      Temp2 = zeros(M,R,N); % preallocating memory
      for m=1:M
        for r=1:R
          Temp1(m,r,:) = CurrentVar(m,r).Mag;
          Temp2(m,r,:) = CurrentVar(m,r).Phase;
        end
      end
      nc{'Data.Mag'}(:) = Temp1;
      nc{'Data.Phase'}(:) = Temp2;
      clear Temp1;
      clear Temp2;
    end

  % ----- NUMERIC VARIABLES ----------------------------------------------
  elseif ndims(CurrentVar)==2
    % positions and vectors
    if size(CurrentVar)==[1 3]                                     % [1 3]
      nc{CurrentVarName} = ncfloat(DimNames.Scalar,DimNames.Coordinates);
    elseif size(CurrentVar)==[M 3]                                 % [M 3]
      nc{CurrentVarName} = ncfloat(DimNames.Measurements,DimNames.Coordinates);
    % "normal" numeric variables
    elseif size(CurrentVar)==[1 1]                                 % [1 1]
      nc{CurrentVarName} = ncfloat(DimNames.Scalar,DimNames.Scalar);
    elseif size(CurrentVar)==[M 1]                                 % [M 1]
      nc{CurrentVarName} = ncfloat(DimNames.Measurements,DimNames.Scalar);
    elseif size(CurrentVar,1)==1                                   % [1 x]
      nc(CurrentVarName) = size(CurrentVar,2);
      nc{CurrentVarName} = ncfloat(DimNames.Scalar,CurrentVarName);
    elseif size(CurrentVar,1)==M                                   % [M x]
      nc(CurrentVarName) = size(CurrentVar,2);
      nc{CurrentVarName} = ncfloat(DimNames.Measurements,CurrentVarName);
    end
    % store variable
    nc{CurrentVarName}(:) = CurrentVar;
  elseif ndims(CurrentVar)==3
    % receiver/transmitter positions
    if size(CurrentVar)==[1 3 1]                                   % [1 3 1]
      nc{CurrentVarName} = ...
        ncfloat(DimNames.Scalar,DimNames.Coordinates,DimNames.Scalar);
    elseif size(CurrentVar)==[1 3 R]                               % [1 3 R]
      nc{CurrentVarName} = ...
        ncfloat(DimNames.Scalar,DimNames.Coordinates,DimNames.Receivers);
    elseif size(CurrentVar)==[M 3 1]                               % [M 3 1]
      nc{CurrentVarName} = ...
        ncfloat(DimNames.Measurements,DimNames.Coordinates,DimNames.Scalar);
    elseif size(CurrentVar)==[M 3 R]                               % [M 3 R]
      nc{CurrentVarName} = ...
        ncfloat(DimNames.Measurements,DimNames.Coordinateso,DimNames.Receivers);
    end
    % store variable
    nc{CurrentVarName}(:) = CurrentVar;
  
  
  else
    error('%s: your variable type is not supported by SOFA',upper(mfilename));
  end

end

close(nc);

% vim:sw=2:ts=2
