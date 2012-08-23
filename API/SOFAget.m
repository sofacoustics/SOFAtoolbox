% SOFAgetData: Read several sets of data from a SOFA file for given IDs
% results = SOFAgetData(Filename,MId)
% Filename specifies the SOFA file from which the data is read.
% MId is a vector of IDs specifying the measurement(s)
% of which the data is read.
% The function returns a cell array with the following structure:
% results{x}{y}{id}
% x ... number of variable
% y = 1: variable name; y = 2: value
% id ... adressing the different measurements

% SOFA API - function SOFAgetData
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

function results = SOFAgetData(Filename,MId)
%% --------------- check and prepare variables ------------------
if(isnumeric(Filename))
  fprintf(2,'Error: Filename must be a string.\n');
  results = 0;
  return;
end
MId = MId(:);
if(~isnumeric(MId) | size(MId,2)>1)
  fprintf(2,'Error: MId must be a numeric scalar or vector.\n');
  results = 0;
  return;
end
if(~iscell(Filename)) Filename = cellstr(Filename); end % assure that Filename is a cell
%% ------------------------- get data ---------------------------
ncid = netcdf.open([Filename{1} '.sofa'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); % get number of variables in file

for id=1:size(MId)
  for ii=0:nvars-1 % loop through all variables in file
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,ii);
    % go through all dimensions of current variable and save their lengths to vectors
    for n=1:size(dimids,2)
      [DimName DimLength(n)] = netcdf.inqDim(ncid,dimids(n));
    end
    % TODO data should be handled by second elseif later (with [M R].FIR(n)]
    if(strcmp(varname,'Data'))
      results{ii+1}{1} = varname;
      results{ii+1}{2}{id} = netcdf.getVar(ncid,ii,[0 MId(id)-1 0],[DimLength(1) 1 DimLength(3)]);
    elseif(size(dimids,2)==1) % scalar
      results{ii+1}{1} = varname;
      results{ii+1}{2}{id} = netcdf.getVar(ncid,ii,0,1);
    elseif(size(dimids,2)==2) % 2-D matrix
      if(DimLength(1)==1) % [1 x]
        results{ii+1}{1} = varname;
        results{ii+1}{2}{id} = netcdf.getVar(ncid,ii,[0 0],[1 DimLength(2)]);
      elseif(DimLength(1)>1) % [M x]
        results{ii+1}{1} = varname;
        results{ii+1}{2}{id} = netcdf.getVar(ncid,ii,[MId(id)-1 0],[1 DimLength(2)]);
      end
    elseif(size(dimids,2)==3) % 3-D matrix
      if(DimLength(1)==1) % [1 x y]
        results{ii+1}{1} = varname;
        results{ii+1}{2}{id} = netcdf.getVar(ncid,ii,[0 0 0],[1 DimLength(2) DimLength(3)]);
      elseif(DimLength(1)>1) % [M x y]
        results{ii+1}{1} = varname;
        results{ii+1}{2}{id} = netcdf.getVar(ncid,ii,[MId(id)-1 0 0],[1 DimLength(2) DimLength(3)]);
      end
    end
  end
end
netcdf.close(ncid);
end