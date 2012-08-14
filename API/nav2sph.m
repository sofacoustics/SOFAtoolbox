% nav2sph: Converts navigational coordinates to spherical coordinates.
% [azi_out,elev] = nav2sph(azi_in,elev)
%
% Input:
%   azi_in ... azimuth angles (-180 <= azi <= 180)
%   elev ... elevation angels (-90 <= elev <= 90)
%
% Output:
%   azi_out ... azimuth angles (0 <= azi < 360)
%   elev ... elevation angels (-90 <= elev <= 90)

% SOFA API - function nav2sph
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

function [azi_out,elev_out]=nav2sph(azi_in,elev_in)
  if(azi_in>-180 && azi_in<0) % azi_in between -180 and 0 deg
    azi_out = azi_in + 360;
  else % azi_in between 0 and 180
    azi_out = azi_in;
  end
  elev_out = elev_in; % elevation remains the same
end % end of function