function [] = SOFAplot(filename,a)
%SOFAPLOT
%
%##########################################################################
%################### IMPORTANT: Not tested (old stuff) ####################
%##########################################################################
%
%   [] = SOFAplot(filename,a) reads one set of data from a SOFA file for a given ID and plot
%   the measurement setup.
%
%   filename specifies the SOFA file from which the data is read.
%   SOFAplot is in a very draft phase. Do not use yet.

% SOFA API - function SOFAplot
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% ----------------- read data from SOFA file --------------------
ncid = netcdf.open(filename,'NC_NOWRITE');

% ---------------------- get variable IDs -------------------------
SubjectIDId = netcdf.inqVarID(ncid,'SubjectID');
SourcePositionTypeId = netcdf.inqVarID(ncid,'SourcePositionType');
SourceViewTypeId = netcdf.inqVarID(ncid,'SourceViewType');
SourceUpTypeId = netcdf.inqVarID(ncid,'SourceUpType');
TransmitterPositionTypeId = netcdf.inqVarID(ncid,'TransmitterPositionType');
ListenerPositionTypeId = netcdf.inqVarID(ncid,'ListenerPositionType');
ListenerViewTypeId = netcdf.inqVarID(ncid,'ListenerViewType');
ListenerUpTypeId = netcdf.inqVarID(ncid,'ListenerUpType');
ReceiverPositionTypeId = netcdf.inqVarID(ncid,'ReceiverPositionType');

SourcePositionId = netcdf.inqVarID(ncid,'SourcePosition');
SourceViewId = netcdf.inqVarID(ncid,'SourceView');
SourceUpId = netcdf.inqVarID(ncid,'SourceUp');
SourceRotationId = netcdf.inqVarID(ncid,'SourceRotation');
TransmitterPositionId = netcdf.inqVarID(ncid,'TransmitterPosition');
ListenerPositionId = netcdf.inqVarID(ncid,'ListenerPosition');
ListenerViewId = netcdf.inqVarID(ncid,'ListenerView');
ListenerUpId = netcdf.inqVarID(ncid,'ListenerUp');
ListenerRotationId = netcdf.inqVarID(ncid,'ListenerRotation');
ReceiverPositionId = netcdf.inqVarID(ncid,'ReceiverPosition');
RoomTypeId = netcdf.inqVarID(ncid,'RoomType');

% ------------------- read values from variables -----------------
oSubjectID = cellstr(netcdf.getVar(ncid,SubjectIDId));
oSourcePositionType = netcdf.getVar(ncid,SourcePositionTypeId);
oSourceViewType = netcdf.getVar(ncid,SourceViewTypeId);
oSourceUpType = netcdf.getVar(ncid,SourceUpTypeId);
oTransmitterPositionType = netcdf.getVar(ncid,TransmitterPositionTypeId);
oListenerPositionType = netcdf.getVar(ncid,ListenerPositionTypeId);
oListenerViewType = netcdf.getVar(ncid,ListenerViewTypeId);
oListenerUpType = netcdf.getVar(ncid,ListenerUpTypeId);
oReceiverPositionType = netcdf.getVar(ncid,ReceiverPositionTypeId);

oSourcePosition = netcdf.getVar(ncid,SourcePositionId);
oSourceView = netcdf.getVar(ncid,SourceViewId);
oSourceUp = netcdf.getVar(ncid,SourceUpId);
oSourceRotation = netcdf.getVar(ncid,SourceRotationId);
oTransmitterPosition = netcdf.getVar(ncid,TransmitterPositionId);
oListenerPosition = netcdf.getVar(ncid,ListenerPositionId);
oListenerView = netcdf.getVar(ncid,ListenerViewId);
oListenerUp = netcdf.getVar(ncid,ListenerUpId);
oListenerRotation = netcdf.getVar(ncid,ListenerRotationId);
oReceiverPosition = netcdf.getVar(ncid,ReceiverPositionId);
oRoomType = netcdf.getVar(ncid,RoomTypeId);

if(strcmp(oRoomType,'shoebox'))
  RoomCornerAId = netcdf.inqVarID(ncid,'RoomCornerA');
  RoomCornerBId = netcdf.inqVarID(ncid,'RoomCornerB');
elseif(strcmp(oRoomType,'collada'))
  DAEFilePath = netcdf.inqVarID(ncid,'DAEFilePath');
end

netcdf.close(ncid)

%% ---------- create data vectors ---------
% ROOM
if(strcmp(oRoomType,'shoebox'))
  RoomCorner1X = RoomCornerAId();
  RoomCorner1Y = -20;
  RoomCorner1Z = -5;
  RoomCorner2X = 15;
  RoomCorner2Y = 20;
  RoomCorner2Z = 12;
elseif(strcmp(oRoomType,'free-field'))
  RoomCorner1X = -15;
  RoomCorner1Y = -20;
  RoomCorner1Z = -5;
  RoomCorner2X = 15;
  RoomCorner2Y = 20;
  RoomCorner2Z = 12;
end



% SOURCE & TRANSMITTERS
% -- SourcePosition
if(size(oSourcePosition,1)==1)
  SourcePositionX = double(oSourcePosition(1,1));
  SourcePositionY = double(oSourcePosition(1,2));
  SourcePositionZ = double(oSourcePosition(1,3));
elseif(size(oSourcePosition,1)>1)
  SourcePositionX = double(oSourcePosition(a,1));
  SourcePositionY = double(oSourcePosition(a,2));
  SourcePositionZ = double(oSourcePosition(a,3));
end
% -- SourceView
if(size(oSourceView,1)==1)
  SourceViewX = double(oSourceView(1,1));
  SourceViewY = double(oSourceView(1,2));
  SourceViewZ = double(oSourceView(1,3));
elseif(size(oSourceView,1)>1)
  SourceViewX = double(oSourceView(a,1));
  SourceViewY = double(oSourceView(a,2));
  SourceViewZ = double(oSourceView(a,3));
end
% -- SourceUp
if(size(oSourceUp,1)==1)
  SourceUpX = double(oSourceUp(1,1));
  SourceUpY = double(oSourceUp(1,2));
  SourceUpZ = double(oSourceUp(1,3));
elseif(size(oSourceUp,1)>1)
  SourceUpX = double(oSourceUp(a,1));
  SourceUpY = double(oSourceUp(a,2));
  SourceUpZ = double(oSourceUp(a,3));
end
% -- SourceRotation
if(size(oSourceRotation,1)==1)
  SourceAzimuth = oSourceRotation(1,1);
  SourceElevation = oSourceRotation(1,2);
  SourceTwist = oSourceRotation(1,3);
else
  SourceAzimuth = oSourceRotation(a,1);
  SourceElevation = oSourceRotation(a,2);
  SourceTwist = oSourceRotation(a,3);
end
% -- TransmitterPosition
if(size(oTransmitterPosition,1)==1)
  TransmitterPositionX(1,:) = double(oTransmitterPosition(1,1,:));
  TransmitterPositionY(1,:) = double(oTransmitterPosition(1,2,:));
  TransmitterPositionZ(1,:) = double(oTransmitterPosition(1,3,:));
elseif(size(oTransmitterPosition,1)>1)
  TransmitterPositionX(a,:) = double(oTransmitterPosition(a,1,:));
  TransmitterPositionY(a,:) = double(oTransmitterPosition(a,2,:));
  TransmitterPositionZ(a,:) = double(oTransmitterPosition(a,3,:));
end

% LISTENER & RECEIVERS
% -- ListenerPosition
if(size(oListenerPosition,1)==1)
  ListenerPositionX = double(oListenerPosition(1,1));
  ListenerPositionY = double(oListenerPosition(1,2));
  ListenerPositionZ = double(oListenerPosition(1,3));
elseif(size(oListenerPosition,1)>1)
  ListenerPositionX = double(oListenerPosition(a,1));
  ListenerPositionY = double(oListenerPosition(a,2));
  ListenerPositionZ = double(oListenerPosition(a,3));
end
% -- ListenerView
if(size(oListenerView,1)==1)
  ListenerViewX = double(oListenerView(1,1));
  ListenerViewY = double(oListenerView(1,2));
  ListenerViewZ = double(oListenerView(1,3));
elseif(size(oListenerView,1)>1)
  ListenerViewX = double(oListenerView(a,1));
  ListenerViewY = double(oListenerView(a,2));
  ListenerViewZ = double(oListenerView(a,3));
end
% -- ListenerUp
if(size(oListenerUp,1)==1)
  ListenerUpX = double(oListenerUp(1,1));
  ListenerUpY = double(oListenerUp(1,2));
  ListenerUpZ = double(oListenerUp(1,3));
elseif(size(oListenerUp,1)>1)
  ListenerUpX = double(oListenerUp(a,1));
  ListenerUpY = double(oListenerUp(a,2));
  ListenerUpZ = double(oListenerUp(a,3));
end
% -- ListenerRotation
if(size(oListenerRotation,1)==1)
  ListenerAzimuth = oListenerRotation(1,1);
  ListenerElevation = oListenerRotation(1,2);
  ListenerTwist = oListenerRotation(1,3);
else
  ListenerAzimuth = oListenerRotation(a,1);
  ListenerElevation = oListenerRotation(a,2);
  ListenerTwist = oListenerRotation(a,3);
end
% -- ReceiverPosition
if(size(oReceiverPosition,1)==1)
  ReceiverPositionX(1,:) = double(oReceiverPosition(1,1,:));
  ReceiverPositionY(1,:) = double(oReceiverPosition(1,2,:));
  ReceiverPositionZ(1,:) = double(oReceiverPosition(1,3,:));
elseif(size(oReceiverPosition,1)>1)
  ReceiverPositionX(a,:) = double(oReceiverPosition(a,1,:));
  ReceiverPositionY(a,:) = double(oReceiverPosition(a,2,:));
  ReceiverPositionZ(a,:) = double(oReceiverPosition(a,3,:));
end

%% ----- coordinate type conversion -----
if(strcmp(oSourcePositionType,'spherical')) % TODO add other coordinate types
  [SourcePositionX SourcePositionY SourcePositionZ] = sph2cart(SourcePositionX,SourcePositionY,SourcePositionZ);
elseif(1)
  
end
if(strcmp(oSourceViewType,'spherical'))
  [SourceViewX SourceViewY SourceViewZ] = sph2cart(degtorad(SourceViewX),degtorad(SourceViewY),SourceViewZ);
elseif(1)
  
end
if(strcmp(oSourceUpType,'spherical'))
  [SourceUpX SourceUpY SourceUpZ] = sph2cart(degtorad(SourceUpX),degtorad(SourceUpY),SourceUpZ);
elseif(1)
  
end
if(strcmp(oTransmitterPositionType,'spherical'))
  [TransmitterPositionX TransmitterPositionY TransmitterPositionZ] = sph2cart(degtorad(TransmitterPositionX),degtorad(TransmitterPositionY),TransmitterPositionZ);
elseif(1)
  
end
if(strcmp(oListenerPositionType,'spherical'))
  [ListenerPositionX ListenerPositionY ListenerPositionZ]
  [ListenerPositionX ListenerPositionY ListenerPositionZ] = sph2cart(degtorad(ListenerPositionX),degtorad(ListenerPositionY),ListenerPositionZ);
elseif(1)
  
end
if(strcmp(oListenerViewType,'spherical'))
  [ListenerViewX ListenerViewY ListenerViewZ] = sph2cart(degtorad(ListenerViewX),degtorad(ListenerViewY),ListenerViewZ);
elseif(1)
  
end
if(strcmp(oListenerUpType,'spherical'))
  [ListenerUpX ListenerUpY ListenerUpZ] = sph2cart(degtorad(ListenerUpX),degtorad(ListenerUpY),ListenerUpZ);
elseif(1)
  
end
if(strcmp(oReceiverPositionType,'spherical'))
  [ReceiverPositionX ReceiverPositionY ReceiverPositionZ] = sph2cart(degtorad(ReceiverPositionX),degtorad(ReceiverPositionY),ReceiverPositionZ);
elseif(1)
  
end

%% ----- Transformations & Rotations -----

% -- transformation of source according to view and up vector
% write values to vectors
SourceView = [SourceViewX SourceViewY SourceViewZ]';
SourceUp = [SourceUpX SourceUpY SourceUpZ]';

SourcePosition = [SourcePositionX SourcePositionY SourcePositionZ]';
TransmitterPosition = [TransmitterPositionX; TransmitterPositionY; TransmitterPositionZ];

SourceOrientation0 = [[1 0 0]' [0 -1 0]' [0 0 1]']; % default orientation is looking in x-direction, with up in z-direction
SourceOrientation1 = [SourceView cross(SourceView,SourceUp) SourceUp];

T = SourceOrientation1 * inv(SourceOrientation0); % get transformation matrix T
% rotation matrices
RotationXII= [1 0 0; 0 cos(SourceTwist) sin(SourceTwist); 0 -sin(SourceTwist) cos(SourceTwist)];
RotationYI = [cos(SourceElevation) 0 -sin(SourceElevation); 0 1 0; sin(SourceElevation) 0 cos(SourceElevation)];
RotationZ = [cos(SourceAzimuth) sin(SourceAzimuth) 0; -sin(SourceAzimuth) cos(SourceAzimuth) 0; 0 0 1];
% TODO fix rotation
for i=1:length(TransmitterPosition(1,:)) % for every transmitter...
  TransmitterPositionOrient(:,i) = T * TransmitterPosition(:,i); % transformation
  TransmitterPositionOrientRot(:,i) = RotationZ*RotationYI*RotationXII*TransmitterPositionOrient(:,i); % rotation
end

% shift transmitters to source position ------------------
TransmitterPositionGlob = TransmitterPosition; %ensure same dimensions...
TransmitterPositionGlob = repmat(SourcePosition,1,length(TransmitterPosition(1,:))) + TransmitterPositionOrientRot;


% -- transformation of listener according to view and up vector
% write values to vectors
ListenerView = [ListenerViewX ListenerViewY ListenerViewZ]';
ListenerUp = [ListenerUpX ListenerUpY ListenerUpZ]';

ListenerPosition = [ListenerPositionX ListenerPositionY ListenerPositionZ]';
ReceiverPosition = [ReceiverPositionX; ReceiverPositionY; ReceiverPositionZ];

ListenerOrientation0 = [[1 0 0]' [0 -1 0]' [0 0 1]']; % default orientation is looking in x-direction, with up in z-direction
ListenerOrientation1 = [ListenerView cross(ListenerView,ListenerUp) ListenerUp];

T = ListenerOrientation1 * inv(ListenerOrientation0); % get transformation matrix T
% rotation matrices
RotationXII= [1 0 0; 0 cos(ListenerTwist) sin(ListenerTwist); 0 -sin(ListenerTwist) cos(ListenerTwist)];
RotationYI = [cos(ListenerElevation) 0 -sin(ListenerElevation); 0 1 0; sin(ListenerElevation) 0 cos(ListenerElevation)];
RotationZ = [cos(ListenerAzimuth) sin(ListenerAzimuth) 0; -sin(ListenerAzimuth) cos(ListenerAzimuth) 0; 0 0 1];
test = [0 1 0]';
for i=1:length(ReceiverPosition(1,:)) % for every Receiver...
  ReceiverPositionOrient(:,i) = T * ReceiverPosition(:,i); % transformation
  ReceiverPositionOrientRot(:,i) = RotationZ*RotationYI*RotationXII*ReceiverPositionOrient(:,i); % rotation
end

% shift Receivers to Listener position ------------------
ReceiverPositionGlob = ReceiverPosition; %ensure same dimensions...
ReceiverPositionGlob = repmat(ListenerPosition,1,length(ReceiverPosition(1,:))) + ReceiverPositionOrientRot;



%% ----- P L O T T I N G -----
fullscreen = get(0,'ScreenSize');
figure('Position',[0 -4 fullscreen(3) fullscreen(4)])
hold on

% ROOM
plot3(RoomCorner1X,RoomCorner1Y,RoomCorner1Z,'marker','o') % corner points
plot3(RoomCorner2X,RoomCorner2Y,RoomCorner2Z,'marker','o')
text(RoomCorner1X+0.8,RoomCorner1Y,RoomCorner1Z,'RoomCorner1') % labels
text(RoomCorner2X+0.8,RoomCorner2Y,RoomCorner2Z,'RoomCorner2')

% -- SOURCE & TRANSMITTERS
plot3(SourcePositionX,SourcePositionY,SourcePositionZ,'marker','o')
text(SourcePositionX+0.5,SourcePositionY,SourcePositionZ,'S')
for n=1:length(TransmitterPosition(1,:))
  plot3(TransmitterPositionGlob(1,n),TransmitterPositionGlob(2,n),TransmitterPositionGlob(3,n),'marker','o')
  text(double(TransmitterPositionGlob(1,n)),double(TransmitterPositionGlob(2,n)+0.5),double(TransmitterPositionGlob(3,n)),['T' num2str(n)])
end
% draw view and up vectors
quiver3(SourcePosition(1),SourcePosition(2),SourcePosition(3),SourceView(1),SourceView(2),SourceView(3))
quiver3(SourcePosition(1),SourcePosition(2),SourcePosition(3),SourceUp(1),SourceUp(2),SourceUp(3))
%text(SourcePosition(1)+SourceView(1)/2,SourcePosition(2)+SourceView(2)/2,SourcePosition(3)+SourceView(3)/2,'View')
%text(SourcePosition(1)+SourceUp(1)/2,SourcePosition(2)+SourceUp(2)/2,SourcePosition(3)+SourceUp(3)/2,'Up')

% -- LISTENER & RECEIVERS
plot3(ListenerPositionX,ListenerPositionY,ListenerPositionZ,'marker','o')
text(ListenerPositionX+0.5,ListenerPositionY,ListenerPositionZ,'L')
for n=1:length(ReceiverPositionX)
  plot3(ReceiverPositionGlob(1,n),ReceiverPositionGlob(2,n),ReceiverPositionGlob(3,n),'marker','o')
  text(double(ReceiverPositionGlob(1,n)),double(ReceiverPositionGlob(2,n)+0.5),double(ReceiverPositionGlob(3,n)),['R' num2str(n)])
end
% draw view and up vectors
quiver3(ListenerPosition(1),ListenerPosition(2),ListenerPosition(3),ListenerView(1),ListenerView(2),ListenerView(3))
quiver3(ListenerPosition(1),ListenerPosition(2),ListenerPosition(3),ListenerUp(1),ListenerUp(2),ListenerUp(3))
%text(ListenerPosition(1)+ListenerView(1)/2,ListenerPosition(2)+ListenerView(2)/2,ListenerPosition(3)+ListenerView(3)/2,'View')
%text(ListenerPosition(1)+ListenerUp(1)/2,ListenerPosition(2)+ListenerUp(2)/2,ListenerPosition(3)+ListenerUp(3)/2,'Up')

grid on
view(37.5,30)
xlabel('X') % axis labels
ylabel('Y')
zlabel('Z')
axis([RoomCorner1X RoomCorner2X RoomCorner1Y RoomCorner2Y RoomCorner1Z RoomCorner2Z])

end % end of function