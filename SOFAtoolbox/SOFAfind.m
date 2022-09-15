function [idx, azinew, elenew, rnew] = SOFAfind(Obj,azi,ele,r)
%SOFAfind - Find a trajectory of a source given the SOFA object
%   Usage: idx = SOFAfind(Obj,azi,ele,r)
%
%   SOFAfind(Obj,azi,ele) returns the indecies idx to the nearest
%   directions available in Obj according to the trajectory given in 
%   the spherical coordinate system described by [azi,ele].
%   For the directions, SourcePositions is used. 
%
%   [idx, azinew, elenew] = SOFAfind(..) returns the actual trajectory
%   described by azinew and elenew. 
%
%   SOFAfind(Obj,azi,ele,r) considers the radius as well. 
%
%   [idx, azinew, elenew, rnew] = SOFAfind(..) returs the actual trajectory
%   described by azinew, elenew, and rnew.
%
%   Input parameters: 
%		Obj:      SOFA object with the positions of the source
%		azi, ele: direction (in degrees) for azimuth and elevation
%       r:        radius (optional). If not provided, radius will be ignored.
% 
%    Output parameters: 
%		idx:      index of the data being nearest to the trajectory
%		azi, ele: azimuth and elevation of the actual direction (degrees)
%       r:        radius of the actual position 
%

% #Author: Piotr Majdak (2019)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
%
% SOFA Toolbox - SOFAfind
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 



%% create a 2D-grid with nearest positions
if ~exist('r','var'),
    r=1;
    Obj.SourcePosition(:,3)=r;
end
pos=SOFAconvertCoordinates(Obj.SourcePosition,'spherical','cartesian');
idx=zeros(length(azi),1);
for ii=1:length(azi)
    [t(:,1),t(:,2),t(:,3)] = sph2cart(deg2rad(azi(ii)),deg2rad(ele(ii)),r);
    dist = sum((pos-repmat(t,size(pos,1),1)).^2,2);
    [~,idx(ii)]=min(dist);
end

%% Output
% actually used angles
azinew = Obj.SourcePosition(idx,1);
elenew = Obj.SourcePosition(idx,2);
rnew = Obj.SourcePosition(idx,3);
end
