function [lat,pol]=sph2hor(azi,ele)
%SPH2HOR  Coordinate Transform
%   [lat,pol]=sph2hor(azi,ele) converts coordinates in spherical format to the horizontal polar
%   one.
% 
%   Input:
%       azi ... azimuth (0 <= azi < 360)
%       ele ... elevation (-90 <= ele <= 90)
% 
%   Output:
%       lat ... lateral angle (-90 <= lat <= 90)
%       pol ... polar angle (-90 <= pol < 270)
%
%   See also SPH2NAV, SPH2VERT, VERT2SPH, NAV2SPH, HOR2SPH

% SOFA API - function sph2hor
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and
% limitations under the Licence. 

azi=azi/180*pi;
ele=ele/180*pi;
lat=zeros(size(azi));
pol=lat;
for ii=1:length(azi)
    lat(ii)=asin(-sin(azi(ii))*cos(ele(ii)));
    pol(ii)=sign(ele(ii))*acos(cos(azi(ii))*cos(ele(ii))/cos(lat(ii)));
end
lat=-lat/pi*180;
pol=pol/pi*180;
end