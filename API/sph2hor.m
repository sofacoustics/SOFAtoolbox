function [lat,pol]=sph2hor(azi,ele)
% 
% sph2hor converts coordinates in spherical format to the horizontal polar
% one.
% 
% Input:
%   azi: azimuth (0 <= azi < 360)
%     
%   ele: elevation (-90 <= ele <= 90)
% 
% Output:
%   lat: lateral angles (-90 <= lat <= 90)
%     
%   pol: polar angles (-90 <= pol < 270)
% 
% Created by Harald Ziegelwanger, OEAW Acoustical Research Institute
% Last Modification: 09.05.2012 by Harald Ziegelwanger
% Name of function changed on 08.13.2012 by Wolfgang Hrauda

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