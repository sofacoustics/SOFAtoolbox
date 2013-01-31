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

function results = NETCDFloadMetadata(Filename,varargin)

%% ----- Checking of input parameters
isargfile(Filename)

ReturnType = 'struct'; % set default value for ReturnType
for ii=1:size(varargin,2)
  if(~all(ischar(varargin{ii})))
    error('Invalid input argument type (must be strings).');
  end
end
if length(varargin)>0 && ...
    (strcmp(varargin{1},'struct') || strcmp(varargin{1},'cell'))
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


%% ----- Main ------------------------------------------------------------
% open file
nc = netcdf(Filename,'r');
% get the variables stored in the file
sofa = ncvar(nc);
nvars = length(sofa);

% Read all metadata variables
count = 0;
for ii=1:nvars % loop through all variables in file

  CurrentVarName = ncname(sofa{ii});

  % omit data variables & only read specified variables or all
  if ~strncmp(CurrentVarName,'Data.',5) && (isempty(varargin) ||...
    sum(strcmp(CurrentVarName,varargin)))
    
    % get variable
    CurrentVar = nc{CurrentVarName}(:);
    if ischar(CurrentVar) CurrentVar = CurrentVar'; end

    % --- return a struct ---
    if(strcmp(ReturnType,'struct'))
      results.(CurrentVarName) = CurrentVar;
    % --- return a cell ---
    elseif(strcmp(ReturnType,'cell'))
      result{1} = CurrentVarName;
      result{2} = CurrentVar;
      results{count + 1} = result;
      count = count + 1;
    end
  end
end

close(nc)

end % of function

% vim:sw=2:ts=2
