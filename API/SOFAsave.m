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
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Piotr Majdak; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

function [] = SOFAsave(Filename,Dataset,varargin)
%global ncid;
%try
%  netcdf.close(ncid)
%catch
%end


%% ----- check format and validity of input variables --------------------
isargchar(Filename)
% V ... number of input variables
if isstruct(Dataset)
  VarNames = fieldnames(Dataset);
  numVar = size(VarNames,1);
elseif iscell(Dataset)
  numVar = size(Dataset,2);
else
  error('Dataset must be struct or cell array.');
end
if(~isempty(varargin) && isnumeric(varargin{1}))
  if(isscalar(varargin{1}) && varargin{1}>=0 && varargin{1}<=9)
    Compression = varargin{1};
  else
    error('Error: Compression must be a numeric scalar value between 0 and 9.');
  end
else
  % setting a default compression
  % FIXME: (hagenw) what value should be the default compression?
  Compression = 9;
end

% ensure Dataset is a struct, by converting the cell
Temp = struct; % temporary variable to create a new struct from cell Dataset
ii = 1;
if(iscell(Dataset)) % -- all of this is only necessary if Dataset is a cell
  while ii<=numVar % loop through all input variables
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
end


% ----- list of mandatory variables and default values -------------------
% first entry of cell is variable name
% second entry of cell is default value. empty string -> mandatory variable
MandatoryVars = {...
  {'Data',''}, ...
  {'DataType','FIR'}, ...
  {'SourcePositionType',''}, ...
  {'SourceViewType',''}, ...
  {'SourceUpType',''}, ...
  {'TransmitterPositionType',''}, ...
  {'ListenerPositionType',''}, ...
  {'ListenerUpType',''}, ...
  {'ReceiverPositionType',''}, ...
  {'RoomType',''}, ...
  {'SamplingRate',''}, ...
  {'ListenerPosition',''}, ...
  {'ListenerView',''}, ...
  {'ListenerUp',''}, ...
  {'ListenerRotation',[0 0 0]}, ...
  {'ReceiverPosition',[0 0 0]}, ...
  {'SourcePosition',''}, ...
  {'SourceView',''}, ...
  {'SourceUp',''}, ...
  {'SourceRotation',[0 0 0]}, ...
  {'TransmitterPosition',[0 0 0]} ...
};

% source/listener variables
SourceListenerVars = {'ListenerPosition','ListenerView','ListenerUp','ListenerRotation', ...
                      'SourcePosition','SourceView','SourceUp','SourceRotation'};

% transmitter/receiver variables
TransmitterReceiverVars = {'ReceiverPosition','TransmitterPosition'};
         

Dataset.SOFAVersion = SOFAgetVersion(); % write SOFA version

% ------ check if mandatroy variables exist ------------------------------
for ii=1:size(MandatoryVars,2) % go through all mandatory variables
  if(~sum(strcmp(MandatoryVars{ii}{1},VarNames))) % if there's no 1 anywhere
    if(~isempty(MandatoryVars{ii}{2})) % if a default value exists
      Dataset = setfield(Dataset,MandatoryVars{ii}{1},MandatoryVars{ii}{2}); % assign name of missing variable
      numVar = numVar + 1; % update number of input variables
    else % if no default value exists
      fprintf(2,['Error: Mandatory variable ' MandatoryVars{ii}{1} ' is missing.\n']);
      return;
    end
  end
end

% ----- get some variable dimensions -------------------------------------
% M ... number of measurements (always encoded in rows, except for data)
% N ... number of samples
% R ... number of receivers
% T ... number of transmitters

VarNames = fieldnames(Dataset); % update VarNames
numVar = size(VarNames,1); % update number of variables
 % length of one piece of data for different data types:
switch Dataset.DataType
	case 'FIR'
%     N = length(Dataset.Data(1,1).FIR);
		[M R N] = size(Dataset.Data.FIR); % retrieve size of data array
	case 'SpectralMagnitudePhase'
%     N = length(Dataset.Data(1,1).Mag);
		[M R N] = size(Dataset.Data.Mag); % retrieve size of data array
	otherwise
		error('Unknown DataType');
end

T = size(Dataset.TransmitterPosition,3); % retrieve number of transmitters

% -------- check matrix dimensions --------
for ii=1:numVar
  CurrentValue = getfield(Dataset,VarNames{ii});
  if(strcmp(VarNames{ii},'Data')) % do nothing (Data sets dimensions)
  elseif(sum(strcmp(SourceListenerVars,VarNames{ii}))) % Source/ListenerVars
    if(~((all(size(CurrentValue)==[M 3])) || (all(size(CurrentValue)==[1 3]))))
     error('%s: Dimensions of coordinate variables are not valid.',VarNames{ii});
    end
  elseif(strcmp(VarNames{ii},'ReceiverPosition')) % ReceiverPosition
    if(~((all(size(CurrentValue)==[1 3 1])) || (all(size(CurrentValue)==[M 3 1]) || ...
         (all(size(CurrentValue)==[1 3 R])) || (all(size(CurrentValue)==[M 3 R])))))
      error('Dimensions of ReceiverPosition (%d %d %d) are not valid.', size(CurrentValue,1), size(CurrentValue,2), size(CurrentValue,3));
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



%% ----- N E T C D F save
NETCDFsave(Filename,Dataset,Compression);


end % of function

% vim:sw=2:ts=2
