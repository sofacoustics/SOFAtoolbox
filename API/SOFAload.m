% SOFAload: Read all data from a SOFA file.
% results = SOFAload(Filename,ReturnType)
% Filename specifies the SOFA file from which the data is read.
% ReturnType is optional and specifies whether the function returns the
% lodaded values as a struct or as a cell array. Default value is 'struct'.
% If ReturnType is 'struct', the function returns a struct which contains
% one field called 'Data' for the data and additional fields for each
% metadata value. The name of these fields are identical to the names of the metadata.
% If ReturnType is 'cell', the function returns a cell array with
% the following structure:
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

function results = SOFAload(Filename,varargin)
%% -- N E T C D F load
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

if(strcmp(ReturnType,'struct'))
  results = SOFAloadMetadata(Filename,'struct');
  results.Data = SOFAloadData(Filename,'struct'); % append data
elseif(strcmp(ReturnType,'cell'))
  results = SOFAloadMetadata(Filename,'cell');
  results{size(results,2)+1} = SOFAloadData(Filename,'cell');
end
end % of function