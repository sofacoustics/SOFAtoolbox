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
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and
% limitations under the Licence. 

if lat==90
    azi=lat;
    ele=0;
    azi=mod(azi+360,360);
else
    lat=deg2rad(mod(lat+360,360));
    pol=deg2rad(mod(pol+360,360));
    ele=asin(cos(lat)*sin(pol));
    if cos(ele)==0
        azi=0;
    else
        azi=real(rad2deg(asin(-sin(lat)/cos(ele))));
    end
    ele=rad2deg(ele);
    if pol > pi/2 && pol< 3*pi/2
        azi=180-azi;
    end
    azi=mod(azi+360,360);
end
end