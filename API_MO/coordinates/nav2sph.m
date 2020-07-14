function [azi,ele]=nav2sph(azi,ele)
%NAV2SPH Coordinate Transform.
%	[azi,ele] = nav2sph(azi,ele) vonverts navigational coordinates to
%	spherical coordinates.
%
%	Input:
%       azi ... azimuth (-180 <= azi <= 180)
%       ele ... elevation (-90 <= ele <= 90)
%
%   Output:
%       azi ... azimuth (0 <= azi < 360)
%       ele ... elevation (-90 <= ele <= 90)
%
%   See also SPH2HOR, SPH2NAV, SPH2VERT, VERT2SPH, HOR2SPH

% SOFA API - function nav2sph
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

idx=find(azi<0); % azi between -180 and 0 deg
azi(idx) = azi(idx)+360;
