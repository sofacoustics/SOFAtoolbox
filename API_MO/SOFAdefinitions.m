function output = SOFAdefinitions(varargin)
% SOFAdefinitions
%
%   SOFAdefinitions returns a struct containing definitions like the time
%   format used in the API.
%
%   SOFAdefinitions('dateFormat') returns the date format
%
%   SOFAdefinitions('APIName') returns the APIName
%
%   SOFAdefinitions('dimensions') returns the dimensions used in the API
%
%   SOFAdefinitions('EOL') returns the end-of-line separator used in the API
%
%   SOFAdefinitions('dateReference') returns the string with the reference
%   for the date when stored as numeric (number of seconds elapsed)
%
%   SOFAdefinitions('units') returns the units and their corresponding aliases
%

% SOFA API - function SOFAdefinitions
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.


definput.flags.type={'all','dateFormat','APIName','dimensions','EOL','dateReference','units'};
[flags,kv]=SOFAarghelper({},definput,varargin);

%% Return all definitions in a structure
if flags.do_all,
  output.APIName = SOFAdefinitions('APIName');
  output.dateFormat = SOFAdefinitions('dateFormat');
  output.dimensions = SOFAdefinitions('dimensions');
  output.EOL = SOFAdefinitions('EOL');
  output.dateReference = SOFAdefinitions('dateReference');
  output.units=SOFAdefinitions('units');
end

%% name of the API
if flags.do_APIName, 
  output = 'ARI SOFA API for Matlab/Octave';
end

%% date string to use (see help datestr)
if flags.do_dateFormat, 
  output = 'yyyy-mm-dd HH:MM:SS';
end

%% EOL to use
if flags.do_EOL, 
  output = char(10);
end

%% dimensions to use
if flags.do_dimensions,
  output.M = 'M'; % Number of Measurements
  output.R = 'R'; % Number of Receivers
  output.N = 'N'; % Number of Samples or the way you represent your data
  output.E = 'E'; % Number of Emitters
  output.C = 'C'; % Coordinates
  output.I = 'I'; % Singleton
  output.S = 'S'; % size of the largest string
end

%% reference for date when used as numeric (number of seconds elapsed)
if flags.do_dateReference, 
  output = '1970-01-01 00:00:00';
end

%% return units with defined aliases
if flags.do_units,
  output.metre={'metres','meter','meters'};
  output.degree={'degrees'};
  output.second={'seconds'};
end
end