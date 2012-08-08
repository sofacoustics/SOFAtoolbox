% SOFAgetData: Read one set of data from a SOFA file for a given ID
% results = SOFAgetData(Filename,MId)
% Filename specifies the SOFA file from which the data is read.
% MId is the number of the measurement of which the data is read.
% The function returns a cell array with the following structure:
% results{x}{y}
% x ... number of variable
% y = 1: variable name; y = 2: value

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
if(~isnumeric(MId) | sum(size(MId))>2)
  fprintf(2,'Error: MId must be a numeric scalar.\n');
  results = 0;
  return;
end

%% ------------------------- get data ---------------------------
ncid = netcdf.open([Filename '.sofa'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); % get number of variables in file

for ii=0:nvars-1 % loop through all variables in file
  [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,ii);
  for n=1:size(dimids,2)
    [DimName DimLength(n)] = netcdf.inqDim(ncid,dimids(n));
  end
  % TODO data should be handled by second elseif later (with [M R].FIR(n)]
  if(strcmp(varname,'Data'))
    results{ii+1} = {varname,netcdf.getVar(ncid,ii,[0 MId-1 0],[DimLength(1) 1 DimLength(3)])};
  elseif(size(dimids,2)==1) % scalar
    results{ii+1} = {varname,netcdf.getVar(ncid,ii,0,1)};
  elseif(size(dimids,2)==2) % 2-D matrix
    if(DimLength(1)==1) % [1 x]
      results{ii+1} = {varname,netcdf.getVar(ncid,ii,[0 0],[1 DimLength(2)])};
    elseif(DimLength(1)>1) % [M x]
      results{ii+1} = {varname,netcdf.getVar(ncid,ii,[MId-1 0],[1 DimLength(2)])};
    end
  elseif(size(dimids,2)==3) % 3-D matrix
    if(DimLength(1)==1) % [1 x y]
      results{ii+1} = {varname,netcdf.getVar(ncid,ii,[0 0 0],[1 DimLength(2) DimLength(3)])};
    elseif(DimLength(1)>1) % [M x y]
      results{ii+1} = {varname,netcdf.getVar(ncid,ii,[MId-1 0 0],[1 DimLength(2) DimLength(3)])};
    end
  end
end
netcdf.close(ncid);
end