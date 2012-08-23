% SOFAloadMetadata: Read all metadata from a SOFA file.
% results = SOFAloadMetadata(Filename,ReturnType,'VarName1','VarName2',...)
% Filename specifies the SOFA file from which the data is read.
% ReturnType is optional and specifies whether the function returns the
% lodaded values as a struct or as a cell array. Default value is 'struct'.
% If ReturnType is 'struct', the function returns a struct which contains
% one field for each metadata value. The name of these fields are identical
% to the names of the metadata.
% If ReturnType is 'cell', the function returns a cell array with
% the following structure:
% results{x}{y}
% x ... number of variable
% y = 1: variable name; y = 2: value
%
% Additionally, an arbitary number of Metadata variable names may be passed
% to the function. In this case, the function only returns the values of
% the specified variables. Otherwise, all Metadata variables will be returned.

% SOFA API - function SOFAloadMetadata
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

function results = SOFAloadMetadata(Filename,varargin)
%% -- N E T C D F load
if(isnumeric(Filename))
  error('Filename must be a string.');
end
for ii=1:size(varargin,2)
  if(~all(ischar(varargin{ii})))
    error('Invalid input argument type (must be strings).');
  end
end
ReturnType = 'struct'; % set default value for ReturnType
if(size(varargin,2)==1)
  varargin = cellstr(varargin);
end
if(strcmp(varargin{1},'struct') || strcmp(varargin{1},'cell'))
   ReturnType = varargin{1};
   varargin(1) = []; % delete ReturnType entry from varargin
end

if(strcmp(ReturnType,'struct'))
  results = struct; % initialize struct variable
elseif(strcmp(ReturnType,'cell'))
  % no need to initialize
else % should not happen anyway, but who knows
  error('ReturnType must be either ''struct'' or ''cell''.');
end
ncid = netcdf.open([char(Filename) '.sofa'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); % get number of variables in file

count = 0;
for ii=0:nvars-1 % loop through all variables in file
  result = {netcdf.inqVar(ncid,ii),netcdf.getVar(ncid,ii)};
  % -- conversion of string variable to cells
  if(~isnumeric(result{2}) && ~isstruct(result{2})) % string variables
    if(size(result{2},1) == 1 & size(result{2},2) >= 1) % [1 len_of_str]
      result{2} = cellstr(strtrim(result{2}));
    else % multidimensional...
      for x=1:size(result{2},2)
        for m=1:size(result{2},1)
          temp{m}{x} = cellstr(strtrim(reshape(result{2}(m,x,:),1,length(result{2}(m,x,:)))));
        end
      end
      result{2} = temp;
    end
  end % -- end of conversion of strings to cell
  % -------------- read METADATA variables ---------------
  % omit data variables & only read specified variables or all, if none given)
  if((~strncmp(result{1},'Data.',5) && (sum(strcmp(result{1},varargin)) || isempty(varargin))))
    if(strcmp(ReturnType,'struct'))
      results.(result{1}) = result{2}; % create variables
    elseif(strcmp(ReturnType,'cell'))
      results{count + 1} = result;
      count = count + 1;
    end
  end
end
netcdf.close(ncid)
end % of function