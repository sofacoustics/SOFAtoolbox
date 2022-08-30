function SOFAstart(flags)
% SOFAstart
%
%   SOFAstart adds all needed pathes and checks if we need the Matlab or Octave
%   version of the API
%
%   SOFAstart(0) or SOFAstart('silent') will suppress any message during the start.
%   SOFAstart ('short') will show a short header only during the start.
%   SOFAstart ('full') will show all information, including all compiled
%   conventions & versions.
%
%   SOFAstart checks if SOFA has been started within the MATLAB session. If
%   it is the case, SOFAstart skips all the initialization. If the initialization
%   is required, SOFAstart('restart') performs the initialization in any case.
%

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: 'full' flag added, changed order of display output messages (11.11.2021)
% #Author: Michael Mihocic: bug fixed when adding paths (29.11.2021)
%
% SOFA Toolbox - function SOFAstart
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% Input parameters
verbose = 2;
restart = 0;
if nargin>0
	if strcmpi(flags,'silent'), verbose=0; end
	if strcmpi(flags,'short'), verbose=1; end
    if strcmpi(flags,'full'), verbose=3; end
	if isnumeric(flags), if flags==0, verbose=0; end; end
    if strcmpi(flags,'restart'), restart=1; end
end

%% do not start when already started but not forced to restart
persistent started
if ~isempty(started) && ~restart
    return;
end
started=1;
%% Check required support
if exist('OCTAVE_VERSION','builtin')
    % We're in Octave
  if compare_versions(OCTAVE_VERSION,'3.6.0','<=')   % check if the octave version is high enough
    error('You need Octave >=3.6.0 to work with SOFA.');
  end
  pkg load netcdf
  if ~which('test_netcdf') % check if octcdf is installed
    error('You have to install the netcdf package in Octave to work with SOFA.');
  end
else
    % We're in Matlab
  if verLessThan('matlab','8')
    warning('SOFA:start','SOFA for Matlab version <= 8 (2012b) not tested. Use on your risk.');
  end
end


%% Add Paths
% Get the basepath as the directory this function resides in.
% The 'which' solution below is more portable than 'mfilename'
% becase old versions of Matlab does not have "mfilename('fullpath')"
basepath=which('SOFAstart');
basepath=basepath(1:end-12); % Kill the function name from the path.
f=filesep;

% Add the base path and the needed sub-directories
% (basepath is added in case user navigated to this folder and changes dir)
if exist('addpath','file') || exist('addpath','builtin') % in Matlab it is a 'file'; in Octave it is a 'built-in' function
  addpath(basepath,[basepath f 'helpers'],[basepath f 'coordinates'],[basepath f 'converters'],[basepath f 'demos'],[basepath f 'netcdf']);
else % in case "addpath" command is not available - can this ever be the case???
  path([basepath f 'helpers'],path);
  path([basepath f 'coordinates'],path);
  path([basepath f 'converters'],path);
  path([basepath f 'demos'],path);
  path([basepath f 'netcdf'],path);
  path(path,basepath);
end


%% Provide SOFA conventions
dispOutput = SOFAcompileConventions;
convs = SOFAgetConventions;

%% Display general informations

if verbose
    disp(['SOFA Toolbox for Matlab/Octave ' SOFAgetVersion '. Copyright: Acoustics Research Institute.']);
    if verbose >= 3
        disp(dispOutput);
    end
	if verbose >= 2
		disp(['This API implements SOFA version ' SOFAgetVersion('SOFA') '.']);
		text=['Available SOFA Conventions: ' convs{1}];
		for ii=2:length(convs)
				text=[text ', ' convs{ii}];
		end
		disp(text);
		disp(['SOFAdbPath (local HRTF database): ' SOFAdbPath('reset') ]);
		disp(['SOFAdbURL (internet repository): ' SOFAdbURL('reset')]);
	end
end



% FIXME: I would check only if the URL is available in the function where it is
% needed. At the start it takes to long. Octaves urlread didn't know the TimeOut
% parameter.
%[~,stat]=urlread(SOFAdbURL);
%if ~stat, disp('  --> could not connect'); end

