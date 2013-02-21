function [azi,ele]=vert2sph(azi,ele)
%VERT2SPH Coordinate Transform
%   [azi,elev] = vert2sph(azi,elev) converts vertical-polar
%   coordinates to spherical coordinates.
%
%   Input:
%       azi ... azimuth (-90 <= azi <= 90)
%       ele ... elevation (-90 <= ele <= 270)
%
%   Output:
%       azi ... azimuth (0 <= azi < 360)
%       ele ... elevation (-90 <= ele <= 90)
%
%   See also SPH2HOR, SPH2NAV, SPH2VERT, NAV2SPH, HOR2SPH

% SOFA API - function vert2sph
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

  if(azi>-90 && azi<0) % azi between -90 and 0 deg
    azi = azi + 360;
  else % azi between 0 and 180
    azi = azi;
  end
  if(elev>90 && elev<=180)
    azi = mod(azi+180,360);
    ele = 180 - ele;
  elseif(ele>180 && ele<=270)
    azi = mod(azi+180,360);
    ele = -(ele -180);
  else
	
  end
end % end of function