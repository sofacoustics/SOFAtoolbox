function APV = SOFAcalculateAPV(Obj)
%SOFAcalculateAPV
%   APV = SOFAcalculateAPV(Obj) calculates the apparent position vector
%   (APV) which represents the position of the source relative to the
%   listener's position and view. APV is in the format [azi ele radius] 
%   with units [deg deg m].
%   Note that ListenerUp is not considered and the APV can be considered as
%   the HRTF direction usually used in HRTF databases

% SOFA API - function SOFAcalculateAPV
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

LP=Obj.ListenerPosition;
LV=Obj.ListenerView;
switch Obj.SourcePosition_Type
  case 'spherical'
    [SPx,SPy,SPz]=sph2cart(deg2rad(Obj.SourcePosition(:,1)),deg2rad(Obj.SourcePosition(:,2)),Obj.SourcePosition(:,3));
    SP=[SPx SPy SPz];
  case 'cartesian'
    SP=Obj.SourcePosition;
end
dist=bsxfun(@minus, SP, LP);
[APVazi,APVele,APVr]=cart2sph(dist(:,1),dist(:,2),dist(:,3));
[LVazi,LVele,~]=cart2sph(LV(:,1),LV(:,2),LV(:,3));
APVazi=bsxfun(@minus, APVazi,LVazi);
APVele=bsxfun(@minus, APVele,LVele);
APV=[rad2deg(APVazi) rad2deg(APVele) APVr];