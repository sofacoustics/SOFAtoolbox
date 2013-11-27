function [azi,ele]=hor2sph(lat,pol)
%HOR2SPH  Coordinate Transform.
%   [azi,ele]=hor2sph(lat,pol) converts coordinates in horizotal polar format to the spherical
%   one.
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

% SOFA API - function hor2sph
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 

if length(lat) ~= length(pol)
  disp('lat and pol must not have different dimensions!')
  return
end

azi = zeros(size(lat));
ele = azi;
for ii = 1:length(lat)  
  if abs(lat(ii))==90
      azi(ii)=lat(ii);
      ele(ii)=0;
      azi(ii)=mod(azi(ii)+360,360);
  else
      lat(ii)=deg2rad(mod(lat(ii)+360,360));
      pol(ii)=deg2rad(mod(pol(ii)+360,360));
      ele(ii)=asin(cos(lat(ii))*sin(pol(ii)));
      if cos(ele(ii))==0
          azi(ii)=0;
      else
          azi(ii)=real(rad2deg(asin(-sin(lat(ii))/cos(ele(ii)))));
      end
      ele(ii)=rad2deg(ele(ii));
      if pol(ii) > pi/2 && pol(ii)< 3*pi/2
          azi(ii)=180-azi(ii);
      end
      azi(ii)=mod(azi(ii)+360,360);
  end
end
end