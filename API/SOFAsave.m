% SOFAsave: Creates a new SOFA file and writes an entire data set to it.
% results = SOFAsave(Filename,{var_name_1,value},var_name_2,value,...)
% Filename specifies the name of the SOFA file to which the data is written.
% The variable names and data can either be given as cell arrays or as
% consecutive arguments (as indicated above). The existence of mandatory
% variables will be checked.
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

function results = SOFAsave(Filename,varargin)
global ncid;
try
  netcdf.close(ncid)
catch
end
%% -- check format and validity of input variables
results = 0;
varargin = varargin{:}; % make "column cell"

% V ... number of input variables
V = size(varargin,2);

% list of mandatory variables
% TODO complete list: orientation and position types, room?
MandatoryVars = {'Data','SamplingRate','ListenerPosition','ListenerView', ...
                 'ListenerUp','ListenerRotation','SourcePosition','SourceView', ...
                 'SourceUp','SourceRotation','ReceiverPosition','TransmitterPosition'};

% source/listener variables
SourceListenerVars = {'ListenerPosition','ListenerView','ListenerUp','ListenerRotation', ...
                      'SourcePosition','SourceView','SourceUp','SourceRotation'};
                    
% transmitter/receiver variables
TransmitterReceiverVars = {'ReceiverPosition','TransmitterPosition'};
         

% ----------- check input variable types -------------
if(~strcmp(cellstr(class(Filename)),'char')) % check type of filename
  fprintf(2,'Error: Filename must be a string.\n');
  return;
end
ii = 1;
while ii<=V % loop through all input variables
  if(~iscell(varargin{ii})) % pack all arguments into cells
    varargin{ii} = {varargin{ii},varargin{ii+1}};
    varargin(ii+1) = []; % delete superfluous entry
    V = size(varargin,2); % reset V
  end
  if(~(mod(size(varargin{ii},2),2)==0)) % variable name and value are not given correctly
    fprintf(2,'Error: Invalid Arguments.\n');
    return;
  end
  if(~strcmp(cellstr(class(varargin{ii}{1})),'char')) % check type of variable name
    varargin{ii}{3}{1};
    fprintf(2,'Error: Invalid Arguments (variable names must be strings).\n');
    return;
  end
  VarNames{ii} = varargin{ii}{1}; % save variable names for mandatory check
  ii = ii + 1;
end

% -------- check if mandatroy variables exist --------
for ii=1:size(MandatoryVars,2) % go through all mandatory variables
  if(~sum(strcmp(MandatoryVars{ii},VarNames))) % if there's no 1 anywhere
    fprintf(2,['Error: Mandatory variable ' MandatoryVars{ii} ' is missing.\n']);
    return;
  end
end

% ----------- get some variable dimensions -----------
% M ... number of measurements (always encoded in rows, except for data)
% N ... number of samples
% R ... number of receivers
% T ... number of transmitters

for ii=1:V % loop through all input variables
  if(strcmp(varargin{ii}{1},'Data'))
    [N M R] = size(varargin{ii}{2}); % retrieve size of data array
  end
  if(strcmp(varargin{ii}{1},'TransmitterPosition'))
   T = size(varargin{ii}{2},3); % retrieve number of transmitters
  end
end

% -------- check matrix dimensions --------
for ii=1:V
  if(strcmp(varargin{ii}{1},'Data')) % do nothing (Data sets dimensions)
  elseif(sum(strcmp(SourceListenerVars,varargin{ii}{1}))) % Source/ListenerVars
    if(~((all(size(varargin{ii}{2})==[M 3])) | (all(size(varargin{ii}{2})==[1 3]))))
      fprintf(2,'Error: Dimensions of coordinate variables are not valid.\n');
      return;
    end
  elseif(strcmp(varargin{ii}{1},'ReceiverPosition')) % ReceiverPosition
    if(~((all(size(varargin{ii}{2})==[1 3 1])) | (all(size(varargin{ii}{2})==[M 3 1]) | ...
         (all(size(varargin{ii}{2})==[1 3 R])) | (all(size(varargin{ii}{2})==[M 3 R])))))
      fprintf(2,'Error: Dimensions of ReceiverPosition are not valid.\n');
      return;
    end
  elseif(strcmp(varargin{ii}{1},'TransmitterPosition')) % TransmitterPosition
    % T doesn't to be checked, as it is defined by the size of TransmitterPosition
    if(~((all(size(varargin{ii}{2})==[1 3])) | (all(size(varargin{ii}{2})==[M 3]))))
      fprintf(2,'Error: Dimensions of TransmitterPosition are not valid.\n');
      return;
    end
  elseif(~(size(size(varargin{ii}{2}),2)>2))
    % if matrix is not more than 2D -> check dimensions: [1 1], [M 1], [1 x], [M x]
    if(~(all(size(varargin{ii}{2})==[1 1]) | all(size(varargin{ii}{2})==[M 1]) | ...
        (size(varargin{ii}{2},1)==1 & size(varargin{ii}{2},2)>1) | ...
        (size(varargin{ii}{2},1)==M & size(varargin{ii}{2},2)>1)))
      fprintf(2,'Error: Invalid matrix dimensions.\n');
      return;
    end
  elseif((size(size(varargin{ii}{2}),2)>2)) % if matrix is more than 2D
    fprintf(2,'Error: Invalid matrix dimensions.\n');
    return;
  end
end
  %% -- N E T C D F save

  ncid = netcdf.create([Filename '.sofa'],'NETCDF4');

% ----------------------- dimensions ---------------------------
% DimId ... vector which contains dimension Ids for every string
% float ... netcdf.getConstant('NC_FLOAT')
% 'ID' might be part of meta data name; 'Id' is a netcdf ID

float = netcdf.getConstant('NC_FLOAT');

MDimId = netcdf.defDim(ncid,'MDim',M);
NDimId = netcdf.defDim(ncid,'NDim',N);
RDimId = netcdf.defDim(ncid,'RDim',R);

ScalarDimId = netcdf.defDim(ncid,'ScalarDim',1);
CoordDimId = netcdf.defDim(ncid,'CoordDim',3);
UnlimitedDimId= netcdf.defDim(ncid,'UnlimitedDim',netcdf.getConstant('NC_UNLIMITED'));

netcdf.endDef(ncid);

% ------- L O O P ---------
for ii=1:V % loop through all input variables
  VarId = 0; % reset VarId (otherwise it becomes a vector!?)
  DimId = 0;
  CurrentVarName = varargin{ii}{1};
  CurrentVar = varargin{ii}{2};
  % --------------- check and prepare variables ------------------
  % -- convert all strings to cells
  if(~isnumeric(CurrentVar)) % if CurrentVar is a string
    CurrentVar = cellstr(CurrentVar);
  end
  
  % dimensions (for length of string) if a cell only contains one string
  if(~isnumeric(CurrentVar)) % -- if CurrentVar is a string
    if(size(CurrentVar,1) == 1 & size(CurrentVar,2) == 1) % [1 1]
      DimId = netcdf.defDim(ncid,[CurrentVarName 'DimId'],length(CurrentVar{1}));
    end
    if(size(CurrentVar,1) == 1 & size(CurrentVar,2) > 1) % [1 x]
      xDimId = netcdf.defDim(ncid,[CurrentVarName 'xDimId'],size(CurrentVar,2));
      for n=1:size(CurrentVar,2) % go through all strings up to x
        lengths(n) = length(CurrentVar{n}); % store all string lengths
      end
      % length of dimension is maximum of all string lengths
      DimId = netcdf.defDim(ncid,[CurrentVarName 'DimId'],max(lengths));
    end
    if(size(CurrentVar,1) == M & size(CurrentVar,2) == 1) % [M 1]
      for n=1:M % go through all strings up to x
        lengths(n) = length(CurrentVar{n});
      end
      DimId = netcdf.defDim(ncid,[CurrentVarName 'DimId'],max(lengths));
    end
    if(size(CurrentVar,1) == M & size(CurrentVar,2) > 1) % [M x]
      xDimId = netcdf.defDim(ncid,[CurrentVarName 'xDimId'],size(CurrentVar,2));
      for n=1:M % go through all strings up to M and x (2D)
        for m=1:size(CurrentVar,2)
          lengths(n,m) = length(CurrentVar{n,m});
        end
      end
      DimId = netcdf.defDim(ncid,[CurrentVarName 'DimId'],max(max(lengths)));
    end
  % dimensions of length x for normal, numeric variables, [1 x] or [M x]
  elseif(~(strcmp(CurrentVarName,'Data') | ... % TODO find more elegant solution...
     strcmp(CurrentVarName,'ListenerPosition') | strcmp(CurrentVarName,'ListenerView') | ...
     strcmp(CurrentVarName,'ListenerUp') | strcmp(CurrentVarName,'ListenerRotation') | ...
     strcmp(CurrentVarName,'SourcePosition') | strcmp(CurrentVarName,'SourceView') | ...
     strcmp(CurrentVarName,'SourceUp') | strcmp(CurrentVarName,'SourceRotation') | ...
     strcmp(CurrentVarName,'ReceiverPosition') | strcmp(CurrentVarName,'TransmitterPosition')))
    if(size(CurrentVar,2) > 1) DimId = netcdf.defDim(ncid,[CurrentVarName 'DimId'],size(CurrentVar,2)); end
  end

  % ------------------------ variables ---------------------------
  if(~isnumeric(CurrentVar)) % --- define string variables ---
  % string variable, [1 1], [1 x], [M 1], [M x]
    if(size(CurrentVar,1) == 1 & size(CurrentVar,2) == 1) % [1 1]
      VarId = netcdf.defVar(ncid,CurrentVarName,2,[ScalarDimId DimId]);
    end
    if(size(CurrentVar,1) == 1 & size(CurrentVar,2) > 1) % [1 x]
      VarId = netcdf.defVar(ncid,CurrentVarName,2,[ScalarDimId xDimId DimId]);
    end
    if(size(CurrentVar,1) == M & size(CurrentVar,2) == 1) % [M 1]
      VarId = netcdf.defVar(ncid,CurrentVarName,2,[MDimId ScalarDimId DimId]);
    end
    if(size(CurrentVar,1) == M & size(CurrentVar,2) > 1) % [M x]
      VarId = netcdf.defVar(ncid,CurrentVarName,2,[MDimId xDimId DimId]);
    end
  
  else % --- define numeric variables ---
    if(strcmp(CurrentVarName,'Data')) % -- Data, float, [N M R]
    VarId = netcdf.defVar(ncid,CurrentVarName,'double',[NDimId MDimId RDimId]);
    
    elseif(strcmp(CurrentVarName,'ListenerPosition') | strcmp(CurrentVarName,'ListenerView') | ...
           strcmp(CurrentVarName,'ListenerUp') | strcmp(CurrentVarName,'ListenerRotation') | ...
           strcmp(CurrentVarName,'SourcePosition') | strcmp(CurrentVarName,'SourceView') | ...
           strcmp(CurrentVarName,'SourceUp') | strcmp(CurrentVarName,'SourceRotation')) % TODO find more elegant solution...
      % -- positions and vectors, float, [1 3] or [M 3]
      if(size(CurrentVar,1) > 1) VarId = netcdf.defVar(ncid,CurrentVarName,float,[MDimId CoordDimId]);
      else VarId = netcdf.defVar(ncid,CurrentVarName,float,[ScalarDimId CoordDimId]); end
      
    elseif(strcmp(CurrentVarName,'ReceiverPosition') | strcmp(CurrentVarName,'TransmitterPosition'))
       % receiver/transmitter position, float, [1 3 1], [1 3 R], [M 3 1] or [M 3 R]
      if((size(CurrentVar,1) == 1) && (size(CurrentVar,3) == 1))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[ScalarDimId CoordDimId ScalarDimId]); end
      if((size(CurrentVar,1) == 1) && (size(CurrentVar,3) > 1))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[ScalarDimId CoordDimId RDimId]); end
      if((size(CurrentVar,1) > 1) && (size(CurrentVar,3) == 1))
        VarId = netcdf.defVar(ncid,CurrentVarName,float,[MDimId CoordDimId ScalarDimId]); end
      if((size(CurrentVar,1) > 1) && (size(CurrentVar,3) > 1))
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
  end
  % ------------------- write values to variables -----------------
  if(~isnumeric(CurrentVar)) % string variables
    if(size(CurrentVar,1) == 1 & size(CurrentVar,2) == 1) % [1 1]
      netcdf.putVar(ncid,VarId,char(CurrentVar));
    end
    if(size(CurrentVar,1) == 1 & size(CurrentVar,2) > 1) % [1 x]
      for n=1:size(CurrentVar,2) % write elements of cell to variable one-by-one
        netcdf.putVar(ncid,VarId,[0 n-1 0],[1 1 length(CurrentVar{n})],CurrentVar{n});
      end
    end
    if(size(CurrentVar,1) == M & size(CurrentVar,2) == 1) % [M 1]
      for n=1:M % write elements of cell to variable one-by-one
        netcdf.putVar(ncid,VarId,[n-1 0 0],[1 1 length(CurrentVar{n})],CurrentVar{n});
      end
    end
    if(size(CurrentVar,1) == M & size(CurrentVar,2) > 1) % [M x]
      for n=1:M % write elements of cell to variable one-by-one
        for m=1:size(CurrentVar,2)
          netcdf.putVar(ncid,VarId,[n-1 m-1 0],[1 1 length(CurrentVar{n,m})],CurrentVar{n,m});
        end
      end
    end
    
  else % numeric variables
    netcdf.putVar(ncid,VarId,CurrentVar);
  end
end

netcdf.close(ncid);
end % of function