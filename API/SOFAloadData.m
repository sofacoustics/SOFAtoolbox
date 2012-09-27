function results = SOFAloadData(Filename,varargin)
%SOFALOADDATA
%   results = SOFAload(Filename,ReturnType) reads all data (no metadata) from
%   a SOFA file.
%
%   Filename specifies the SOFA file from which the data is read.
%   ReturnType is optional and specifies whether the function returns the
%   lodaded values as a struct or as a cell array. Default value is 'struct'.
%   If ReturnType is 'struct', the function returns a structure array which
%   contains the data from the SOFA file.
%   If ReturnType is 'cell', the function returns a cell array with the same
%   contents as in the case of 'struct' but with the following structure:
%   results{y}
%   y = 1: variable name; y = 2: value

% SOFA API - function SOFAloadData
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% --------------- check and prepare variables ------------------
if(isnumeric(Filename))
  error('Filename must be a string.');
end
ReturnType = 'struct'; % set default value for ReturnType
if(size(varargin,2)==1)
  varargin = cellstr(varargin);
  ReturnType = varargin{1};
end
if(isnumeric(ReturnType))
  error('ReturnType must be a string.');
end

if(strcmp(ReturnType,'struct'))
  results = struct; % initialize struct variable
elseif(strcmp(ReturnType,'cell'))
  % no need to initialize
else
  error('ReturnType must be either ''struct'' or ''cell''.');
end

%% ---------------------- N E T C D F load ----------------------
ncid = netcdf.open([char(Filename) '.sofa'],'NC_NOWRITE');
try
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); % get number of variables in file

% ---------- LOOP through all variables in file ----------
for ii=0:nvars-1
  result{1} = netcdf.inqVar(ncid,ii); % get current variable name
  % ----------------- read DATA variables -----------------
  if(strncmp(result{1},'Data.',5)) % if current variables is a Data variable
    result{2} = netcdf.getVar(ncid,ii); % get values from current variable
    CurrentFieldName = result{1}(6:end); % field name is without 'Data.'
    % ------ return type: 'struct' ------
    if(strcmp(ReturnType,'struct'))
      results.(CurrentFieldName) = result{2};
      % ------ return type: 'cell' ------
    elseif(strcmp(ReturnType,'cell'))
      results{2}.(CurrentFieldName) = result{2};
      results{1} = 'Data'; % variable name is 'Data' only (if returning cell)
    end
  end
end
catch
  if(exist('ncid','var') && ~isempty(ncid)) netcdf.close(ncid); end
  error(['An error occured during reading the SOFA file: ' lasterr()]);
  % TODO lasterr() should not be used any more...
end
netcdf.close(ncid)

end % of function