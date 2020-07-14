function phi = deg2rad(phi)
%DEG2RAD returns the given angle in radians
%   phi = deg2rad(phi) converts the angle phi given in degree to radians
%
%   Input options:
%       phi     - angle, can be a scalar or matrix (degree)
%
%   Output options:
%       phi     - angle (rad)
%
%   See also: RAD2DEG

% SOFA API - function deg2rad
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
phi = phi./180*pi;
