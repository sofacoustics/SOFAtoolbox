% SOFAsave: Creates a new SOFA file and writes an entire data set to it.
% [] = SOFAsave(Filename,Dataset,Compression)
% Filename specifies the name of the SOFA file to which the data is written.
% Dataset is either a struct or a cell array containing the data and meta
% data to be written to the SOFA file (see below for exact format).
% Compression is an optional numeric value between 0 and 9 specifying the
% amount of compression to be applied to the data when writing to the netCDF file.
% 0 is no compression and 9 is the most compression.
% 
% If Dataset is a struct, it must contain one field called 'Data' for the data
% and additional fields for each metadata value. The name of these fields
% are identical to the names of the metadata.
% If Dataset is a cell, it must have the following structure:
% Dataset{x}{y}
% x ... number of variable
% y = 1: variable name; y = 2: value
% 
% In both cases, the existence of mandatory variables will be checked.
% Coordinate variables are expected to have one of the following
% dimensions (with M being the number of measurements):
% Source/ListenerPosition, -View, -Up: [1 3], [M 3]
% Transmitter/ReceiverPosition: [1 3 1], [1 3 R], [M 3 1], [M 3 R]
% (with R being the number of receivers or transmitters respectively)
%
% All other meta data variables must have one of the following dimensions:
% [1 1], [1 x], [M 1], [M x] (x is arbitary)

% SOFA API - function SOFAsave
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

function [] = SOFAsave(Filename,Dataset,varargin)
%% -- check format and validity of input variables

% V ... number of input variables
if(strcmp(class(Dataset),'struct'))
  VarNames = fieldnames(Dataset);
  V = size(VarNames,1);
elseif(strcmp(class(Dataset),'cell'))
  V = size(Dataset,2);
end

% -- list of mandatory variables and default values
% first entry of cell is variable name
% second entry of cell is default value. empty string -> mandatory variable
MandatoryVars = {{'Data',''},{'DataType','FIR'},{'SourcePositionType',''}, ...
  {'SourceViewType',''},{'SourceUpType',''},{'TransmitterPositionType',''}, ...
  {'ListenerPositionType',''},{'ListenerUpType',''},{'ReceiverPositionType',''}, ...
  {'RoomType',''},{'SamplingRate',''},{'ListenerPosition',''},{'ListenerView',''}, ...
  {'ListenerUp',''},{'ListenerRotation',[0 0 0]},{'ReceiverPosition',[0 0 0]}, ...
  {'SourcePosition',''},{'SourceView',''},{'SourceUp',''},{'SourceRotation',[0 0 0]}, ...
  {'TransmitterPosition',[0 0 0]}};

% source/listener variables
SourceListenerVars = {'ListenerPosition','ListenerView','ListenerUp','ListenerRotation', ...
                      'SourcePosition','SourceView','SourceUp','SourceRotation'};
                    
% transmitter/receiver variables
TransmitterReceiverVars = {'ReceiverPosition','TransmitterPosition'};
         

% ----------- check input variable types -------------
if(~strcmp(cellstr(class(Filename)),'char')) % check type of filename
  error('Filename must be a string.');
end

if(~isempty(varargin) && isnumeric(varargin{1}))
  if(isscalar(varargin{1}) && varargin{1}>=0 && varargin{1}<=9)
    Compression = varargin{1};
  else
    error('Error: Compression must be a numeric scalar value between 0 and 9.');
  end
end

Temp = struct; % temporary variable to create a new struct from cell Dataset
ii = 1;
if(iscell(Dataset)) % -- all of this is only necessary if Dataset is a cell
  while ii<=V % loop through all input variables
    if(~(mod(size(Dataset{ii},2),2)==0)) % variable name and value are not given correctly
      fprintf(2,'Error: Invalid Arguments.\n');
      return;
    end
    if(~strcmp(cellstr(class(Dataset{ii}{1})),'char')) % check type of variable name
      %Dataset{ii}{3}{1}; % TODO???
      fprintf(2,'Error: Invalid Arguments (variable names must be strings).\n');
      return;
    end
    Temp = setfield(Temp,Dataset{ii}{1},Dataset{ii}{2});
    VarNames{ii} = Dataset{ii}{1}; % save variable names for mandatory check
    ii = ii + 1;
  end % -- end of while loop
  Dataset = Temp; % Dataset is now a struct
  clear Temp; % free memory

% if dataset is neither cell nor struct (cell has already been checked)
elseif(~isstruct(Dataset))
  error('Dataset must be struct or cell array.');
end

Dataset.SOFAVersion = SOFAgetVersion(); % write SOFA version

% -------- check if mandatroy variables exist --------
for ii=1:size(MandatoryVars,2) % go through all mandatory variables
  if(~sum(strcmp(MandatoryVars{ii}{1},VarNames))) % if there's no 1 anywhere
    if(~isempty(MandatoryVars{ii}{2})) % if a default value exists
      Dataset = setfield(Dataset,MandatoryVars{ii}{1},MandatoryVars{ii}{2}); % assign name of missing variable
      V = V + 1; % update number of input variables
    else % if no default value exists
      fprintf(2,['Error: Mandatory variable ' MandatoryVars{ii}{1} ' is missing.\n']);
      return;
    end
  end
end

% ----------- get some variable dimensions -----------
% M ... number of measurements (always encoded in rows, except for data)
% N ... number of samples
% R ... number of receivers
% T ... number of transmitters

VarNames = fieldnames(Dataset); % update VarNames
V = size(VarNames,1); % update number of variables
if(strcmp(Dataset.DataType,'FIR')) [M R N] = size(Dataset.Data.FIR); % retrieve size of data array
elseif(strcmp(Dataset.DataType,'SpectralMagnitudePhase')) [M R N] = size(Dataset.Data.Mag);
% >->-><> add other Data types here in the future
else
  error('Unkown DataType.');
end
T = size(Dataset.TransmitterPosition,3); % retrieve number of transmitters

% -------- check matrix dimensions --------
for ii=1:V
  CurrentValue = Dataset.(VarNames{ii});
  if(strcmp(VarNames{ii},'Data')) % do nothing (Data sets dimensions)
  elseif(sum(strcmp(SourceListenerVars,VarNames{ii}))) % Source/ListenerVars
    if(~((all(size(CurrentValue)==[M 3])) || (all(size(CurrentValue)==[1 3]))))
     error('Dimensions of coordinate variables are not valid.');
    end
  elseif(strcmp(VarNames{ii},'ReceiverPosition')) % ReceiverPosition
    if(~((all(size(CurrentValue)==[1 3 1])) || (all(size(CurrentValue)==[M 3 1]) || ...
         (all(size(CurrentValue)==[1 3 R])) || (all(size(CurrentValue)==[M 3 R])))))
      error('Dimensions of ReceiverPosition are not valid.');
    end
  elseif(strcmp(VarNames{ii},'TransmitterPosition')) % TransmitterPosition
    % T doesn't need to be checked, as it is defined by the size of TransmitterPosition
    if(~((all(size(CurrentValue)==[1 3])) || (all(size(CurrentValue)==[M 3]))))
      error('Dimensions of TransmitterPosition are not valid.');
    end
  elseif(~(size(size(CurrentValue),2)>2))
    % if matrix is not more than 2D -> check dimensions: [1 1], [M 1], [1 x], [M x]
    if(~(all(size(CurrentValue)==[1 1]) || all(size(CurrentValue)==[M 1]) || ...
        (size(CurrentValue,1)==1 && size(CurrentValue,2)>1) || ...
        (size(CurrentValue,1)==M && size(CurrentValue,2)>1)))
      error('Invalid matrix dimensions.');
    end
  elseif((size(size(CurrentValue),2)>2)) % if matrix is more than 2D
    error('Invalid matrix dimensions.');
  end
end
  %% -- N E T C D F save
ncid = netcdf.create([Filename '.sofa'],'NETCDF4');
try

% ----------------------- dimensions ---------------------------
% DimId ... vector which contains dimension Ids for every string
% float ... netcdf.getConstant('NC_FLOAT')
% 'ID' might be part of meta data name; 'Id' is a netcdf ID

float = netcdf.getConstant('NC_FLOAT');

MDimId = netcdf.defDim(ncid,'MDim',M);
RDimId = netcdf.defDim(ncid,'RDim',R);
NDimId = netcdf.defDim(ncid,'NDim',N);
TDimId = netcdf.defDim(ncid,'TDim',T);

% >->-><> define additional dimensions for future data types here

ScalarDimId = netcdf.defDim(ncid,'ScalarDim',1);
CoordDimId = netcdf.defDim(ncid,'CoordDim',3);
%UnlimitedDimId= netcdf.defDim(ncid,'UnlimitedDim',netcdf.getConstant('NC_UNLIMITED'));

netcdf.endDef(ncid);

% ------- L O O P ---------
for ii=1:V % loop through all input variables
  VarId = 0; % reset VarId (otherwise it becomes a vector!?)
  DimId = 0;
  CurrentVarName = VarNames{ii};
  CurrentVar = Dataset.(VarNames{ii});
  % --------------- check and prepare variables ------------------
  % -- convert all strings to cells
  if(~isnumeric(CurrentVar) && ~isstruct(CurrentVar)) % if CurrentVar is a string
    CurrentVar = cellstr(CurrentVar);
  end
  
  % dimensions (for length of string) if a cell only contains one string
  if(~isnumeric(CurrentVar) && ~isstruct(CurrentVar)) % -- if CurrentVar is a string
    if(size(CurrentVar,1) == 1 && size(CurrentVar,2) == 1) % [1 1]
      DimId = netcdf.defDim(ncid,[CurrentVarName 'DIM'],length(CurrentVar{1}));
    end
    if(size(CurrentVar,1) == 1 && size(CurrentVar,2) > 1) % [1 x]
      xDimId = netcdf.defDim(ncid,[CurrentVarName 'xDIM'],size(CurrentVar,2));
      for n=1:size(CurrentVar,2) % go through all strings up to x
        lengths(n) = length(CurrentVar{n}); % store all string lengths
      end
      % length of dimension is maximum of all string lengths
      DimId = netcdf.defDim(ncid,[CurrentVarName 'DIM'],max(lengths));
    end
    if(size(CurrentVar,1) == M && size(CurrentVar,2) == 1) % [M 1]
      for n=1:M % go through all strings up to x
        lengths(n) = length(CurrentVar{n});
      end
      DimId = netcdf.defDim(ncid,[CurrentVarName 'DIM'],max(lengths));
    end
    if(size(CurrentVar,1) == M && size(CurrentVar,2) > 1) % [M x]
      xDimId = netcdf.defDim(ncid,[CurrentVarName 'xDIM'],size(CurrentVar,2));
      for n=1:M % go through all strings up to M and x (2D)
        for m=1:size(CurrentVar,2)
          lengths(n,m) = length(CurrentVar{n,m});
        end
      end
      DimId = netcdf.defDim(ncid,[CurrentVarName 'DIM'],max(max(lengths)));
    end
  % dimensions of length x for normal, numeric variables, [1 x] or [M x]
    elseif(~(strcmp(CurrentVarName,'Data') | sum(strcmp(CurrentVarName,SourceListenerVars)) | ...
             sum(strcmp(CurrentVarName,TransmitterReceiverVars))))
      if(size(CurrentVar,2) > 1) DimId = netcdf.defDim(ncid,[CurrentVarName 'xDIM'],size(CurrentVar,2)); end
  end

  % ------------------------ variables ---------------------------
  if(~isnumeric(CurrentVar) && ~isstruct(CurrentVar)) % --- define string variables ---
  % string variable, [1 1], [1 x], [M 1], [M x]
    if(size(CurrentVar,1) == 1 && size(CurrentVar,2) == 1) % [1 1]
      VarId = netcdf.defVar(ncid,CurrentVarName,2,[ScalarDimId DimId]);
    end
    if(size(CurrentVar,1) == 1 && size(CurrentVar,2) > 1) % [1 x]
      VarId = netcdf.defVar(ncid,CurrentVarName,2,[ScalarDimId xDimId DimId]);
    end
    if(size(CurrentVar,1) == M && size(CurrentVar,2) == 1) % [M 1]
      VarId = netcdf.defVar(ncid,CurrentVarName,2,[MDimId ScalarDimId DimId]);
    end
    if(size(CurrentVar,1) == M && size(CurrentVar,2) > 1) % [M x]
      VarId = netcdf.defVar(ncid,CurrentVarName,2,[MDimId xDimId DimId]);
    end
  
  else % --- define numeric variables ---
    if(strcmp(CurrentVarName,'Data')) % -- Data, float, [M R N]
      if(strcmp(Dataset.DataType,'FIR'))
        VarId = netcdf.defVar(ncid,'Data.FIR','double',[MDimId RDimId NDimId]);
      elseif(strcmp(Dataset.DataType,'SpectralMagnitudePhase'))
        VarIdMag = netcdf.defVar(ncid,'Data.Mag','double',[MDimId RDimId NDimId]);
        VarIdPhase = netcdf.defVar(ncid,'Data.Phase','double',[MDimId RDimId NDimId]);
      
      % >->-><> define additional variables for future data types here
      %         (variable names must be 'Data.xxx'!)
      end    
    elseif(sum(strcmp(CurrentVarName,SourceListenerVars)))
      % -- positions and vectors, float, [1 3] or [M 3]
      if(size(CurrentVar,1) > 1) VarId = netcdf.defVar(ncid,CurrentVarName,float,[MDimId CoordDimId]);
      else VarId = netcdf.defVar(ncid,CurrentVarName,float,[ScalarDimId CoordDimId]); end
      
    elseif(sum(strcmp(CurrentVarName,TransmitterReceiverVars)))
       % receiver/transmitter position, float, [1 3 1], [1 3 R], [M 3 1] or [M 3 R]
      if((size(CurrentVar,1) == 1) && (size(CurrentVar,3) == 1))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[ScalarDimId CoordDimId ScalarDimId]); end
      if((size(CurrentVar,1) == 1) && (size(CurrentVar,3) > 1))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[ScalarDimId CoordDimId RDimId]); end
      if((size(CurrentVar,1) > 1) && (size(CurrentVar,3) == 1))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[MDimId CoordDimId ScalarDimId]); end
      if((size(CurrentVar,1) > 1) && (size(CurrentVar,3) > 1) && strcmp(CurrentVarName,'TransmitterPosition'))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[MDimId CoordDimId TDimId]); end      
      if((size(CurrentVar,1) > 1) && (size(CurrentVar,3) > 1) && strcmp(CurrentVarName,'ReceiverPosition'))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[MDimId CoordDimId RDimId]); end
    else % "normal" numeric variables, float, [1 1], [M 1], [1 x], [M x]
      if((size(CurrentVar,1) == 1) && (size(CurrentVar,2) == 1))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,ScalarDimId); end
      if((size(CurrentVar,1) == M) && (size(CurrentVar,2) == 1))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[MDimId ScalarDimId]); end
      if((size(CurrentVar,1) == 1) && (size(CurrentVar,2) > 1))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[ScalarDimId DimId]); end
      if((size(CurrentVar,1) == M) && (size(CurrentVar,2) > 1))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[MDimId DimId]); end
    end
    if(exist('Compression','var'))
      netcdf.defVarDeflate(ncid,VarId,true,true,Compression);
    end
  end
  % ------------------- write values to variables -----------------
  if(~isnumeric(CurrentVar) && ~isstruct(CurrentVar)) % string variables
    if(size(CurrentVar,1) == 1 && size(CurrentVar,2) == 1) % [1 1]
      netcdf.putVar(ncid,VarId,char(CurrentVar));
    end
    if(size(CurrentVar,1) == 1 && size(CurrentVar,2) > 1) % [1 x]
      for n=1:size(CurrentVar,2) % write elements of cell to variable one-by-one
        netcdf.putVar(ncid,VarId,[0 n-1 0],[1 1 length(CurrentVar{n})],CurrentVar{n});
      end
    end
    if(size(CurrentVar,1) == M && size(CurrentVar,2) == 1) % [M 1]
      for n=1:M % write elements of cell to variable one-by-one
        netcdf.putVar(ncid,VarId,[n-1 0 0],[1 1 length(CurrentVar{n})],CurrentVar{n});
      end
    end
    if(size(CurrentVar,1) == M && size(CurrentVar,2) > 1) % [M x]
      for n=1:M % write elements of cell to variable one-by-one
        for m=1:size(CurrentVar,2)
          netcdf.putVar(ncid,VarId,[n-1 m-1 0],[1 1 length(CurrentVar{n,m})],CurrentVar{n,m});
        end
      end
    end
  elseif(strcmp(CurrentVarName,'Data')) % write data variables
    if(strcmp(Dataset.DataType,'FIR'))
      netcdf.putVar(ncid,VarId,CurrentVar.FIR);
    elseif(strcmp(Dataset.DataType,'SpectralMagnitudePhase'))
      netcdf.putVar(ncid,VarIdMag,CurrentVar.Mag);
      netcdf.putVar(ncid,VarIdPhase,CurrentVar.Phase);
      
    % >->-><> write values of data for future data types here
    end
    
  else % numeric variables
    netcdf.putVar(ncid,VarId,CurrentVar);
  end
end
catch
  if(exist('ncid','var') && ~isempty(ncid)) netcdf.close(ncid); end
  error(['An error occured during reading the SOFA file: ' lasterr()]);
  % TODO lasterr() should not be used any more...
end
netcdf.close(ncid);
end % of function