%% ARI to Sofa format conversion
% load ARI .mat file
load('NH30 HRTFs.mat')

oFilename = 'NH30';

% convert audio channel to corresponding twist angle (setup at ARI lab)
angles = [-30 -20 -10 0 10 20 30 40 50 60 70 80 -25 -15 -5 5 15 25 35 45 55 65];
angles = pi*(angles/180); % convert to rad
for n=1:size(hM,2)
  Twist(n,1) = angles(meta.pos(n,3));
end

oData = {'Data',hM};
oDataType = {'DataType','FIR'};
oOrientationType = {'OrientationType','Cartesian'};
oPositionType = {'PositionType','Cartesian'};
oSamplingRate = {'SamplingRate',stimPar.SamplingRate};
oSubjectID = {'SubjectID',cellstr(stimPar.SubjectID)};
oApplicationName = {'ApplicationName','test application name'};
oApplicationVersion = {'ApplicationVersion',stimPar.Version};
oSourcePosition = {'SourcePosition',[0 0 0]};
oSourceView = {'SourceView',[1 0 0]};
oSourceUp = {'SourceUp',[0 0 1]};
oSourceRotation = {'SourceRotation',[0 0 0]};
oTransmitterPosition = {'TransmitterPosition',[0 0 0]};
oListenerPosition = {'ListenerPosition',[5 0 0]};
oListenerView = {'ListenerView',[-1 0 0]};
oListenerUp = {'ListenerUp',[0 0 1]};
oListenerRotation = {'ListenerRotation',[meta.pos(:,1) meta.pos(:,2) Twist]};
oReceiverPosition = {'ReceiverPosition',zeros(1,3,2)};
oReceiverPosition{2}(:,:,1) = [0 1 0];
oReceiverPosition{2}(:,:,2) = [0 -1 0];
oMeasurementID = {'MeasurementID',stimPar.ID};
oMeasurementParameterSourceAudioChannel = {'MeasurementParameterSourceAudioChannel',meta.pos(:,3)};
oMeasurementParameterItemIndex = {'MeasurementParameterItemIndex',meta.itemidx};
oMeasurementParameterAudioLatency = {'MeasurementParameterAudioLatency',meta.lat};
oMeasurementParameterSourceAmplitude = {'MeasurementParameterSourceAmplitude',meta.amp};

varargin = {oFilename,oData,oDataType,oOrientationType,oPositionType,oSamplingRate, ...
    oSubjectID,oApplicationName,oApplicationVersion,oSourcePosition,oSourceView, ...
    oSourceUp,oSourceRotation,oTransmitterPosition,oListenerPosition,oListenerView, ...
    oListenerUp,oListenerRotation,oReceiverPosition,oMeasurementID, ...
    oMeasurementParameterSourceAudioChannel,oMeasurementParameterItemIndex, ...
    oMeasurementParameterAudioLatency,oMeasurementParameterSourceAmplitude};

%% -- N E T C D F save

% TODO here:
% assure standard format of varargin (cells and "single variables")
% check if mandatory variables exist etc.
% check matrix dimensions

ncid = netcdf.create([varargin{1} '.nc'],'NETCDF4');

% ----------------------- dimensions ---------------------------
% M ... number of measurements (always encoded in rows, except for data)
% N ... number of samples
% R ... number of receivers
% T ... number of transmitters
% DimId ... vector which contains dimension Ids for every string
% float ... netcdf.getConstant('NC_FLOAT')
% 'ID' might be part of meta data name; 'Id' is a netcdf ID

for ii=2:size(varargin,2) % loop through all input variables
  if(strcmp(varargin{ii}{1},'Data'))
    [N M R] = size(varargin{ii}{2}); % retrieve size of data array
  end
  if(strcmp(varargin{ii}{1},'TransmitterPosition'))
   T = size(varargin{ii}{2},3); % retrieve number of transmitters
  end
end

MDimId = netcdf.defDim(ncid,'MDim',M);
NDimId = netcdf.defDim(ncid,'NDim',N);
RDimId = netcdf.defDim(ncid,'RDim',R);
float = netcdf.getConstant('NC_FLOAT');

ScalarDimId = netcdf.defDim(ncid,'ScalarDim',1);
CoordDimId = netcdf.defDim(ncid,'CoordDim',3);
UnlimitedDimId= netcdf.defDim(ncid,'UnlimitedDim',netcdf.getConstant('NC_UNLIMITED'));

netcdf.endDef(ncid);

% ------- L O O P ---------
for ii=2:size(varargin,2) % loop through all input variables
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
  if(~isnumeric(CurrentVar)) % if CurrentVar is a string
    if(size(CurrentVar,1) == 1) DimId = netcdf.defDim(ncid,[CurrentVarName 'DimId'],length(CurrentVar{1})); end
  % dimensions of length x for normal, numeric variables, [1 x] or [M x]
  elseif(~(strcmp(CurrentVarName,'Data') | ...
     strcmp(CurrentVarName,'ListenerPosition') | strcmp(CurrentVarName,'ListenerView') | ...
     strcmp(CurrentVarName,'ListenerUp') | strcmp(CurrentVarName,'ListenerRotation') | ...
     strcmp(CurrentVarName,'SourcePosition') | strcmp(CurrentVarName,'SourceView') | ...
     strcmp(CurrentVarName,'SourceUp') | strcmp(CurrentVarName,'SourceRotation') | ...
     strcmp(CurrentVarName,'ReceiverPosition') | strcmp(CurrentVarName,'TransmitterPosition')))
    if(size(CurrentVar,2) > 1) DimId = netcdf.defDim(ncid,[CurrentVarName 'DimId'],size(CurrentVar,2)); end
  end

  % ------------------------ variables ---------------------------
  if(~isnumeric(CurrentVar)) % --- define string variables ---
  % string variable, single [1 length-of-string] or [M unlimited]
    if(size(CurrentVar,1) > 1) VarId = netcdf.defVar(ncid,CurrentVarName,2,[MDimId UnlimitedDimId]);
    else VarId = netcdf.defVar(ncid,CurrentVarName,2,[ScalarDimId DimId]); end
  
  else % --- define numeric variables ---
    if(strcmp(CurrentVarName,'Data')) % -- Data, float, [N M R]
    VarId = netcdf.defVar(ncid,CurrentVarName,'double',[NDimId MDimId RDimId]);
    
    elseif(strcmp(CurrentVarName,'ListenerPosition') | strcmp(CurrentVarName,'ListenerView') | ...
           strcmp(CurrentVarName,'ListenerUp') | strcmp(CurrentVarName,'ListenerRotation') | ...
           strcmp(CurrentVarName,'SourcePosition') | strcmp(CurrentVarName,'SourceView') | ...
           strcmp(CurrentVarName,'SourceUp') | strcmp(CurrentVarName,'SourceRotation'))
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
    if(size(CurrentVar,1) > 1)
      for n=1:M % write elements of cell to variable one-by-one
        netcdf.putVar(ncid,VarId,[n-1 0],[1 length(CurrentVar{n})],CurrentVar{n});
      end
    elseif(size(CurrentVar,1) == 1) netcdf.putVar(ncid,VarId,char(CurrentVar));
    end
    
  else % numeric variables
    netcdf.putVar(ncid,VarId,CurrentVar);
  end
end

netcdf.close(ncid);

%% -- N E T C D F load

ncid = netcdf.open([varargin{1} '.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); % get number of variables in file

for ii=0:nvars-1 % loop through all variables in file
  results{ii+1} = {netcdf.inqVar(ncid,ii),netcdf.getVar(ncid,ii)};
end

netcdf.close(ncid)