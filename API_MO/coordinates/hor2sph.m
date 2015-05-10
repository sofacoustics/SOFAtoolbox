function [azi,ele]=hor2sph(lat,pol)
%HOR2SPH  transform horizontal-polar to spherical coordinates.
%   [azi,ele]=hor2sph(lat,pol)
% 
%   Input:
%       lat ... lateral angle (-90 <= lat <= 90)
%       pol ... polar angle (-90 <= pol < 270)
% 
%   Output:
%       azi ... azimuth (0 <= azi < 360)
%       ele ... elevation (-90 <= ele <= 90)
%
%   See also SPH2HOR, SPH2NAV, SPH2VERT, VERT2SPH, NAV2SPH

% AUTHOR: Robert Baumgartner

% SOFA API - function hor2sph
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 

% interpret horizontal polar format as rotated spherical coordinates with
% negative azimuth direction
[x,nz,y] = sph2cart(-deg2rad(pol),deg2rad(lat),ones(size(lat)));

[azi,ele,r] = cart2sph(x,y,-nz);

azi = rad2deg(azi);
ele = rad2deg(ele);

% adjust azimuth range 
[azi,ele] = nav2sph(azi,ele);

end