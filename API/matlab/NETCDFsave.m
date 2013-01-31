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

% source/listener variables
SourceListenerVars = {'ListenerPosition','ListenerView','ListenerUp','ListenerRotation', ...
                      'SourcePosition','SourceView','SourceUp','SourceRotation'};

% transmitter/receiver variables
TransmitterReceiverVars = {'ReceiverPosition','TransmitterPosition'};


%% ----- N E T C D F save ------------------------------------------------
% FIXME: discuss what mode of NETCDF files we should create
try
	ncid = netcdf.create(Filename,'NETCDF4');

% define some constants and fixed dimensions
float = netcdf.getConstant('NC_FLOAT');
ScalarDimId = netcdf.defDim(ncid,DimNames.Scalar,1);
CoordDimId = netcdf.defDim(ncid,DimNames.Coordinates,3);
UnlimitedDimId= netcdf.defDim(ncid,DimNames.Unlimited,netcdf.getConstant('NC_UNLIMITED'));

netcdf.endDef(ncid);

% ------- L O O P ---------
for ii=1:numVar % loop through all input variables
  VarId = 0; % reset VarId (otherwise it becomes a vector!?)
  DimId = 0;
  CurrentVarName = VarNames{ii}
  CurrentVar = getfield(Dataset,VarNames{ii});
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
      switch Dataset.DataType
				case 'FIR'
					M = size(CurrentVar.FIR,1);
					R = size(CurrentVar.FIR,2);
					N = size(CurrentVar.FIR,3);
					MDimId = netcdf.defDim(ncid,DimNames.Measurements,M);
					RDimId = netcdf.defDim(ncid,DimNames.Receivers,R);
					NDimId = netcdf.defDim(ncid,DimNames.Samples,N);
					VarId = netcdf.defVar(ncid,'Data.FIR','double',[MDimId RDimId NDimId]);
				case 'SpectralMagnitudePhase'
					M = size(CurrentVar.Mag,1);
					R = size(CurrentVar.Mag,2);
					N = size(CurrentVar.Mag,3);
					MDimId = netcdf.defDim(ncid,DimNames.Measurements,M);
					RDimId = netcdf.defDim(ncid,DimNames.Receivers,R);
					NDimId = netcdf.defDim(ncid,DimNames.Samples,N);
					VarIdMag = netcdf.defVar(ncid,'Data.Mag','double',[MDimId RDimId NDimId]);
					VarIdPhase = netcdf.defVar(ncid,'Data.Phase','double',[MDimId RDimId NDimId]);
				otherwise
					error('Unknown data type');
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
    % Compression of data
    netcdf.defVarDeflate(ncid,VarId,true,true,Compression);
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
    switch Dataset.DataType
			case 'FIR'
				Temp = zeros(M,R,N); % preallocating memory
				for m=1:M
					for r=1:R
						Temp(m,r,:) = CurrentVar.FIR(m,r);
					end
				end
				netcdf.putVar(ncid,VarId,Temp);
				clear Temp;
			case 'SpectralMagnitudePhase'
				Temp1 = zeros(M,R,N); % preallocating memory
				Temp2 = zeros(M,R,N); % preallocating memory
				for m=1:M
					for r=1:R
						Temp1(m,r,:) = CurrentVar.Mag(m,r);
						Temp2(m,r,:) = CurrentVar.Phase(m,r);
					end
				end
				netcdf.putVar(ncid,VarIdMag,Temp1);
				netcdf.putVar(ncid,VarIdPhase,Temp2);
				clear Temp1;
				clear Temp2;
		
    % >->-><> write values of data for future data types here
    end
    
  else % numeric variables
    netcdf.putVar(ncid,VarId,CurrentVar);
  end
end

netcdf.close(ncid);

catch me
	netcdf.abort(ncid);
	throw(me)
end


% vim:sw=2:ts=2
