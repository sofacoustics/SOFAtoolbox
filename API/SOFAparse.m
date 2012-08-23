% SOFAgetID: Searches measurements where a certain variable has a certain
% value and returns the Id's of these measurements as a vector.
% results = SOFAgetID(Filename,TargetVarName,TargetValue,Equality, ...
% TargetValueRange)
% Filename specifies the SOFA file which is read.
% TargetVarName is the name of the variable of which certain values are
% looked for.
% TargetValue is the value of TargetVarName that will be looked for.
% Equality is optional and can be '=', '==', '<', '>', '<=' or '>='.
% Values that are equal to/greater/less than TargetValue will be returned.
% The default value for Equality is '='. '=' and '==' are equivalent.
% TargetValueRange specifies a +/- offset to TargetValue. Values within
% that range will also be returned by the function. TargetValueRange can
% either be a scalar or a vector matching the dimension of TargetValue.
% results is a column vector that contains all Id's that have been found.

% SOFA API - function SOFAgetID
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

function results = SOFAgetID(Filename,TargetVarName,TargetValue,varargin)
%% --------------- check and prepare variables ------------------
results = 0; % default value
TargetValueRange = 0; % default value for TargetValueRange
TargetCoordinate = 0; % default value for TargetCoordinate
 % TODO all sorts of checks: exist, filetypes, dimensions

% checking varargin
if(~all(size(varargin)== [0 0])) % if varargin is NOT empty -> check
%   varargin = varargin{:} % make "column cell"
  % V ... number of input variables
  V = size(varargin,2);
  if(1) % TODO all sorts of checks: exist, filetypes, dimensions
  end
  if(~iscell(varargin)) varargin = cellstr(varargin); end % make sure it is a cell
  if(~isnumeric(varargin{1}) && (strcmp(varargin{1},'=') | strcmp(varargin{1},'<') | ...
     strcmp(varargin{1},'>') | strcmp(varargin{1},'<=') | (strcmp(varargin{1},'==') | ...
     strcmp(varargin{1},'>='))))
    Equality = varargin{1};
    varargin(1) = []; % delete first entry if it was equality
  else
    Equality = '='; % default value for equality, if not passed by user
  end
  V = size(varargin,2);
  for ii=1:V
    if(strcmp(varargin{ii}{1},'TargetValueRange')) TargetValueRange = varargin{ii}{2};
    end
    if(strcmp(varargin{ii}{1},'TargetCoordinate')) TargetCoordinate = varargin{ii}{2};
    end
  end
else % if varargin is empty
  Equality = '='; % default value for equality
  TargetValueRange = 0; % default value for TargetValueRange
  TargetCoordinate = 0; % default value for TargetCoordinate
end
% only search for given coordinate
if(~(TargetCoordinate==0)) TargetValue = TargetValue(TargetCoordinate); end
count = 1; % counts number of entries that match TargetValue

% TODO   coordinate type conversion?
%% ----------------------- loop and search ----------------------
global ncid;

dataset = SOFAload(Filename);
V = size(dataset,2); % get number of variables
for n=1:V % -- go through all variables and search TargetVarName
  if(strcmp(dataset{n}{1},TargetVarName))      
    % -- go through all measurements and check equality
    for m=1:size(dataset{n}{2},1)
      if(TargetCoordinate==0) result = dataset{n}{2}(m,:);
      else
        result = dataset{n}{2}(m,TargetCoordinate);
      end
      result = round(result*10000);
      result = cast(result,'int64');
      result = cast(result,'double');
      result = result/10000;
      
      %TargetValue
      if( ((strcmp(Equality,'=') | strcmp(Equality,'==')) && (all(result<=TargetValue+TargetValueRange) ...
         && all(result>=TargetValue-TargetValueRange))) | ...
      (strcmp(Equality,'<') && all(result<(TargetValue+TargetValueRange))) | ...
      (strcmp(Equality,'>') && all(result>(TargetValue-TargetValueRange))) | ...
      (strcmp(Equality,'<=') && all(result<=(TargetValue+TargetValueRange))) | ...          
      (strcmp(Equality,'>=') && all(result>=(TargetValue-TargetValueRange))) )
        results(count,:) = m;
        count = count + 1;
      end
    end
  end
end
end