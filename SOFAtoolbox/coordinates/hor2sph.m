function [azi,ele]=hor2sph(lat,pol)
%HOR2SPH - Transform horizontal-polar to spherical coordinates.
%   Usage: [azi,ele]=hor2sph(lat,pol)
% 
%   Input:
%       lat : lateral angle (-90 <= lat <= 90)
%       pol : polar angle (-90 <= pol < 270)
% 
%   Output:
%       azi : azimuth (0 <= azi < 360)
%       ele : elevation (-90 <= ele <= 90)
%
%   See also SPH2HOR, SPH2NAV, SPH2HOR, SPH2SH, NAV2SPH

% #Author: Robert Baumgartner
% #Author: Michael Mihocic: header documentation updated (28.10.2021)

% SOFA Toolbox - function hor2sph
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% interpret horizontal polar format as rotated spherical coordinates with
% negative azimuth direction
[x,nz,y] = sph2cart(-deg2rad(pol),deg2rad(lat),ones(size(lat)));

[azi,ele,r] = cart2sph(x,y,-nz);

azi = rad2deg(azi);
ele = rad2deg(ele);

% adjust azimuth range 
[azi,ele] = nav2sph(azi,ele);

end