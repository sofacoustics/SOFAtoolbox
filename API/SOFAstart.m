% SOFAstart added all needed pathes and checks if we need the Matlab or Octave
% version of the API
%
%   Usage: SOFAstart

% AUTHOR: Hagen Wierstorf


%% ----- Adding Path's ---------------------------------------------------
% Get the basepath as the directory this function resides in.
% The 'which' solution below is more portable than 'mfilename'
% becase old versions of Matlab does not have "mfilename('fullpath')"
basepath=which('SOFAstart');
% Kill the function name from the path.
basepath=basepath(1:end-12);
% Add the base path and the needed sub-directories
if exist('addpath')
  addpath(basepath);
  addpath([basepath,'/helper']);
  if isoctave
    addpath([basepath,'/octave']);
  else
    addpath([basepath,'/matlab']);
  end
else
  path(path,basepath);
  path(path,[basepath,'/helper']);
  if isoctave
    path(path,[basepath,'/octave']);
  else
    path(path,[basepath,'/matlab']);
  end
end

% vim:sw=2:ts=2
