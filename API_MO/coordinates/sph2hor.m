function [lat,pol]=sph2hor(azi,ele)
%SPH2HOR  Coordinate Transform
%   [lat,pol]=sph2hor(azi,ele) converts coordinates in spherical format to the horizontal polar
%   one.
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

% SOFA API - function sph2hor
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and
% limitations under the Licence. 

% azi=azi/180*pi;
% ele=ele/180*pi;
% lat=zeros(size(azi));
% pol=lat;
% for ii=1:length(azi)
%     lat(ii)=asin(-sin(azi(ii))*cos(ele(ii)));
%     pol(ii)=sign(ele(ii))*acos(cos(azi(ii))*cos(ele(ii))/cos(lat(ii)));
% end
% lat=-lat/pi*180;
% pol=pol/pi*180;
azi=mod(azi+360,360);
ele=mod(ele+360,360);

razi = azi*2*pi/360;
rele = ele*2*pi/360;
rlat=asin(sin(razi).*cos(rele));
rpol=zeros(size(rlat));
idx=find(cos(rlat)~=0);
rpol(idx)=real(asin(sin(rele(idx))./cos(rlat(idx))));
pol = rpol*360/2/pi;
lat = rlat*360/2/pi;

idx = find(razi>pi/2 & razi < 3*pi/2 & (rele < pi/2 | rele > 3*pi/2));
pol(idx)=180-pol(idx);
idx = find(~(razi>pi/2 & razi < 3*pi/2) & rele > pi/2 & rele < 3*pi/2);
pol(idx)=180-pol(idx);
