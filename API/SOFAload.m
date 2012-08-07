% SOFAload: Read all data from a SOFA file.
% results = SOFAgetData(Filename)
% Filename specifies the SOFA file from which the data is read.
% The function returns a cell array with the following structure:
% results{x}{y}
% x ... number of variable
% y = 1: variable name; y = 2: value

% SOFA API - function SOFAload
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

function results = SOFAload(Filename)
%% -- N E T C D F load
if(isnumeric(Filename))
  fprintf(2,'Error: Filename must be a string.\n');
  results = 0;
  return;
end

ncid = netcdf.open([Filename '.sofa'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); % get number of variables in file

for ii=0:nvars-1 % loop through all variables in file
  results{ii+1} = {netcdf.inqVar(ncid,ii),netcdf.getVar(ncid,ii)};
end
netcdf.close(ncid)
end % of function