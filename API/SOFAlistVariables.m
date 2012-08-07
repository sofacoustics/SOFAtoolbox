% SOFAlistVariables: Reads specified variables from a SOFA file.
% results = SOFAlistVariables(Filename,'var_name_1','var_name_2',...)
% The function opens the SOFA file, reads the data of
% the specified variables and returns the results as a cell with the
% following structure:
% results{x}{y}
% x ... number of variable
% y = 1: variable name; y = 2: value

% SOFA API - function SOFAlistVariables
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

function results = SOFAlistVariables(Filename,varargin)
if(isnumeric(Filename))
  fprintf(2,'Error: Filename must be a string.\n');
  results = 0;
  return;
end

ncid = netcdf.open([Filename '.sofa'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); % get number of variables in file
V = size(varargin,2); % get number of requested variables

try % issue error if MDim doesn't exist
  MDimId = netcdf.inqDimID(ncid,'MDim'); % get number of measurements
catch
  fprintf(2,'Error: Invalid SOFA file.\n');
end

[DimName M] = netcdf.inqDim(ncid,MDimId);
count = 1;

for ii=0:nvars-1 % loop through all variables in file
  for n=0:V-1 % loop through all requested variables
    if(isnumeric(varargin{n+1})) % check if input arguments are strings
      fprintf(2,'Error: Input arguments must be strings.\n');
      netcdf.close(ncid);
      results = 0;
      return;
    end
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,ii);
    if(strcmp(varargin{n+1},varname))
      results{count} = {varname,netcdf.getVar(ncid,ii)}; % return cell
      count = count + 1;
    end
  end
end

netcdf.close(ncid)

end % of function