function results = SOFAgetData(Filename,MId,varargin)
%SOFAGETDATA
%   results = SOFAload(Filename,MId,ReturnType) reads several sets of data 
%   from a SOFA file for given IDs
%
%   Filename specifies the SOFA file from which the data is read.
%   MId is a vector of IDs specifying the measurement(s)
%   of which the data is read.
%   ReturnType is optional and specifies whether the function returns the 
%   lodaded values as a struct or as a cell array. Default value is 'struct'.
%   If ReturnType is 'struct', the function returns a structure array which 
%   contains the data from the SOFA file.
%   If ReturnType is 'cell', the function returns a cell array with the same 
%   contents as in the case of 'struct' but with the following structure: 
%   results{y}{id}
%   y = 1: variable name; y = 2: value
%   id ... adressing the different measurements

% SOFA API - function SOFAgetData
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

MId = MId(:);
if(~isnumeric(MId) | size(MId,2)>1)
  error('MId must be a numeric scalar or vector.');
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
ncid = netcdf.open([char(Filename) '.sofa'],'NC_NOWRITE');
try
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); % get number of variables in file

%% ------------------------- get data ---------------------------
for id=1:size(MId) % --- LOOP through all requested Ids
  for ii=0:nvars-1 % --- LOOP through all variables in file
    % get current variable name and dimension Ids
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,ii);
    % go through all dimensions of current variable and save their lengths to vectors
    for n=1:size(dimids,2)
      [DimName DimLength(n)] = netcdf.inqDim(ncid,dimids(n));
    end
    M = DimLength(1); R = DimLength(2); N = DimLength(3);
    if(strncmp(varname,'Data.',5)) % if current variable is a Data variable
      CurrentFieldName = varname(6:end); % field name is without 'Data.'
      result = netcdf.getVar(ncid,ii,[MId(id)-1 0 0],[1 R N]); % read a single data entry
      % ------ return type: 'struct' ------
      if(strcmp(ReturnType,'struct'))
        results.(CurrentFieldName){id} = result;
%         for r=1:R
%           results(r).(CurrentFieldName){id} = reshape(result(1,r,:),1,N);
%         end
        % ------ return type: 'cell' ------
      elseif(strcmp(ReturnType,'cell'))
        results{1} = 'Data'; % variable name is 'Data' only (if returning cell)
        results{2}.(CurrentFieldname){id} = result;
%         for r=1:R
%           results{2}(r).(CurrentFieldName){id} = reshape(result(1,r,:),1,N);
%         end
      end
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