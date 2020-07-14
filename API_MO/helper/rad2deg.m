function phi = rad2deg(phi)
%RAD2DEG returns the given angle in degree
%   phi = rad2deg(phi) converts the angle phi given in radians to degree
%
%   Input options:
%       phi     - angle, can be a scalar or matrix (rad)
%
%   Output options:
%       phi     - angle (degree)
%
%   See also: deg2rad

% SOFA API - function rad2deg
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 

% Checking of input parameters
narginchk(1,1);

% Convert angle
phi = phi./pi*180;
