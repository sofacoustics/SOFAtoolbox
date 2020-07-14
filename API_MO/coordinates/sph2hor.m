function [lat,pol]=sph2hor(azi,ele)
%SPH2HOR  transform spherical to horizontal-polar coordinates.
%   [lat,pol]=sph2hor(azi,ele)
% 
%   Input:
%       azi ... azimuth (in degrees)
%       ele ... elevation (in degrees)
% 
%   Output:
%       lat ... lateral angle (-90 <= lat <= 90)
%       pol ... polar angle (-90 <= pol < 270)
%
%   See also SPH2NAV, SPH2VERT, VERT2SPH, NAV2SPH, HOR2SPH

% AUTHOR: Robert Baumgartner

% SOFA API - function sph2hor
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 

[x,y,z] = sph2cart(deg2rad(azi),deg2rad(ele),ones(size(azi)));

% remove noise below eps
x(abs(x)<eps)=0;
y(abs(y)<eps)=0;
z(abs(z)<eps)=0;

% interpret horizontal polar format as rotated spherical coordinates with
% negative azimuth direction
[pol,nlat,r] = cart2sph(x,z,-y);
pol = rad2deg(pol);
lat = rad2deg(-nlat);

% adjust polar angle range
pol = mod(pol+90,360)-90;
end