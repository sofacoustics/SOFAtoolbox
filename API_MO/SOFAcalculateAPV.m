function ApparentPositionVector = SOFAcalculateAPV(Obj)
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


% === Apparent azimuth ===
% Apparent azimuth is the relative direction of the sound source in the
% horizontal plane from the listener viewpoint
%
% Get head orientation of the listener and source position in spherical
% coordinates
HeadOrientation = SOFAconvertCoordinates( ...
    Obj.ListenerView,Obj.ListenerView_Type,'spherical');
SourcePosition = SOFAconvertCoordinates( ...
    Obj.SourcePosition,Obj.SourcePosition_Type,'spherical');
% Calculate the relative azimuth angle between them
APVazimuth = correctAzimuth(bsxfun( ...
    @minus,SourcePosition(:,1),HeadOrientation(:,1)));

% === Apparent elevation ===
% Apparent elevation is the relative direction of the sound source in the median
% plane from the listener viewpoint
APVelevation = correctElevation(bsxfun( ...
    @minus,SourcePosition(:,2),HeadOrientation(:,2)));

% === Apparent distance ===
% Apparent distance is the relative distance between the sound source and the
% listener
%
% Get listener positon in spherical coordinates
ListenerPosition = SOFAconvertCoordinates( ...
    Obj.ListenerPosition,Obj.ListenerPosition_Type,'spherical');
% Relative distance
APVdistance = bsxfun(@minus,SourcePosition(:,3),ListenerPosition(:,3));

% Combine to matrix
ApparentPositionVector = [ ...
    APVazimuth, ...
    APVelevation, ...
    APVdistance, ...
];
    
end

function phi = correctAzimuth(phi)
    % Ensure -360 <= phi <= 360
    phi = rem(phi,360);
    % Ensure -180 <= phi < 180
    phi(phi<-180) = phi(phi<-180) + 360;
    phi(phi>=180) = phi(phi>=180) - 360;
end

% TODO: check what convetion we are using for delta!
function delta = correctElevation(delta)
    % Ensure -180 <= delta <= 180
    delta = correctAzimuth(delta);
    % Ensure -90 <= delta <= 90
    delta(delta<-90) = -delta(delta<-90) - 180;
    delta(delta>90) = -delta(delta>90) + 180;
end
