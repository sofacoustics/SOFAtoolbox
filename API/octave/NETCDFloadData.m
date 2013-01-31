% SOFAloadData: Read all data (no metadata) from a SOFA file.
% results = SOFAload(Filename,ReturnType)
% Filename specifies the SOFA file from which the data is read.
% ReturnType is optional and specifies whether the function returns the
% lodaded values as a struct or as a cell array. Default value is 'struct'.
% If ReturnType is 'struct', the function returns a structure array which
% contains the data from the SOFA file.
% If ReturnType is 'cell', the function returns a cell array with the same
% contents as in the case of 'struct' but with the following structure:
% results{y}
% y = 1: variable name; y = 2: value

% SOFA API - function SOFAloadData
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

function results = NETCDFloadData(Filename,varargin)
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


%% ----- Main -------------------------------------------------------------
% open file
nc = netcdf(Filename,'r');
% get the variables stored in the file
sofa = ncvar(nc);
nvars = length(sofa);

%% ---------- LOOP through all variables in file ----------
for ii=1:nvars

  CurrentVarName = ncname(sofa{ii});

  % ----------------- read DATA variables -----------------
  if strncmp(CurrentVarName,'Data.',5) % if current variables is a Data variable
    CurrentVar = nc{CurrentVarName}(:);
    M = size(CurrentVar,1); % get size of data array
    R = size(CurrentVar,2);
    N = size(CurrentVar,3);
    CurrentFieldName = CurrentVarName(6:end); % field name is without 'Data.'
    % ------ return type: 'struct' ------
    if(strcmp(ReturnType,'struct'))
      for m=1:M
        for r=1:R
          results(m,r).(CurrentFieldName) = reshape(CurrentVar(m,r,:),1,N);
        end
      end
      % ------ return type: 'cell' ------
    elseif(strcmp(ReturnType,'cell'))
      for m=1:M
        for r=1:R
          results{2}(m,r).(CurrentFieldName) = reshape(CurrentVar(m,r,:),1,N);
        end
      end
      results{1} = 'Data'; % variable name is 'Data' only (if returning cell)
    end
  end
end

close(nc)
end % of function

% vim:sw=2:ts=2
