function SOFAstart
% SOFAstart 
%
%   SOFAstart adds all needed pathes and checks if we need the Matlab or Octave
%   version of the API

% SOFA API - function SOFAstart
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence.

%% Display general informations
disp(['SOFA Matlab/Octave API version ' SOFAgetVersion '. Copyright 2013 Acoustics Research Institute (piotr@majdak.com).']);
disp(['This API implements SOFA version ' SOFAgetVersion('SOFA') '.']);
convs=SOFAgetConventions;
text=['Available conventions: ' convs{1}];
for ii=2:length(convs)
	text=[text ', ' convs{ii}];
end
disp(text);
disp(['Location of the HRTF database: ' SOFAdbPath ]);

%% ---------------------------- Adding Path's -----------------------------
% Get the basepath as the directory this function resides in.
% The 'which' solution below is more portable than 'mfilename'
% becase old versions of Matlab does not have "mfilename('fullpath')"
basepath=which('SOFAstart');
% Kill the function name from the path.
basepath=basepath(1:end-12);
f=filesep;
% Add the base path and the needed sub-directories
if exist('addpath','builtin')
  addpath(basepath);
  addpath([basepath f 'helper']);
  addpath([basepath f 'coordinates']);
  addpath([basepath f 'converters']);
  addpath([basepath f 'demos']);
  if isoctave
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
  if isoctave
    path(path,[basepath f 'octave']);
  else
    path(path,[basepath f 'matlab']);
  end
end