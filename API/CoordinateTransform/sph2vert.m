function [azi,ele] = sph2vert(azi,ele)
%SPH2VERT  Coordinate Transform
%   [azi,ele] = sph2vert(azi,ele) converts vertical-polar
%   coordinates to spherical coordinates.
%
%   Input:
%       azi ... azimuth (0 <= azi <= 360)
%       ele ... elevation (-90 <= ele <= 90)
%
%   Output:
%       azi ... azimuth (-90 <= azi <= 90)
%       ele ... elevation (-90 <= ele <= 270)
%
%   See also SPH2HOR, SPH2NAV, VERT2SPH, NAV2SPH, HOR2SPH

% SOFA API - function sph2vert
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

  if(azi>=270 && azi_in<360) % azi_in between 270 and 360 deg
    azi = -(360-azi);
    ele = ele;
  elseif(azi>=180 && azi<270) % azi_in between 180 and 270 deg
    azi = azi - 180;
    ele = 180 - ele;
  elseif(azi>=90 && azi<180) % azi_in between 90 and 180 deg
    azi = azi - 180;
    ele = 180 - ele;
  else % azi_in between 0 and 90 deg and outside 360 deg
    azi = azi;
    ele = ele;
  end
end % end of function