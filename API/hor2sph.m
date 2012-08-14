function [azi,ele]=hor2sph(lat,pol)
% 
% hor2sph converts coordinates in horizotal polar format to the spherical
% one.
% 
% Input:
%   lat: lateral angles (-90 <= lat <= 90)
%     
%   pol: polar angles (-90 <= pol < 270)
% 
% Output:
%   azi: azimuth (0 <= azi < 360)
%     
%   ele: elevation (-90 <= ele <= 90)
% 
% Created by Harald Ziegelwanger, OEAW Acoustical Research Institute
% Last Modification: 09.05.2012 by Harald Ziegelwanger
% Name of function changed on 08.13.2012 by Wolfgang Hrauda

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