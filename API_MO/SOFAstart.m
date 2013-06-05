function SOFAstart
% SOFAstart 
%
%   SOFAstart adds all needed pathes and checks if we need the Matlab or Octave
%   version of the API

% SOFA API - function SOFAstart
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% Check required support
if exist('OCTAVE_VERSION','builtin')
    % We're in Octave
  if compare_versions(OCTAVE_VERSION,'3.6.0','<=')   % check if the octave version is high enough
    error('You need Octave >=3.6.0 to work with SOFA.');
  end  
  if ~which('netcdf') % check if octcdf is installed
    error('You have to install the octcdf package in Octave to work with SOFA.');
  end
else
    % We're in Matlab
  if verLessThan('matlab','8')
    warning('SOFA for Matlab version >= 8 (2012b) not tested. Use on your risk.');
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
if exist('addpath','builtin')
  addpath(basepath);
  addpath([basepath f 'helper']);
  addpath([basepath f 'coordinates']);
  addpath([basepath f 'converters']);
  addpath([basepath f 'demos']);
  if exist('OCTAVE_VERSION','builtin')
    addpath([basepath f 'octave']);
  else
    addpath([basepath f 'matlab']);
  end
else
  path(path,basepath);
  path(path,[basepath f 'helper']);
  path(path,[basepath f 'coordinates']);
  path(path,[basepath f 'converters']);
  path(path,[basepath f 'demos']);
  if exist('OCTAVE_VERSION','builtin')
    path(path,[basepath f 'octave']);
  else
    path(path,[basepath f 'matlab']);
  end
end


%% Display general informations
disp(['SOFA Matlab/Octave API version ' SOFAgetVersion '. Copyright 2013 Acoustics Research Institute (piotr@majdak.com).']);
disp(['This API implements SOFA version ' SOFAgetVersion('SOFA') '.']);
SOFAcompileConventions;
convs=SOFAgetConventions;
text=['Available SOFA Conventions: ' convs{1}];
for ii=2:length(convs)
	text=[text ', ' convs{ii}];
end
disp(text);
disp(['SOFAdbPath (local HRTF database): ' SOFAdbPath ]);
disp(['SOFAdbURL (internet repository): ' SOFAdbURL]);
% FIXME: I would check only if the URL is available in the function where it is
% needed. At the start it takes to long. Octaves urlread didn't know the TimeOut
% parameter.
%[~,stat]=urlread(SOFAdbURL);
%if ~stat, disp('  --> could not connect'); end

