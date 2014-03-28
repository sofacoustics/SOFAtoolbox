function SOFAplotGeometry(Obj, index)
% SOFAplotGeometry(Obj) plots the geometry found in the Obj.
% 
% SOFAplotGeometry(Obj, index) plots the geometry for the measurements
% given in the index. 

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


if ~exist('index','var')
  index=1:Obj.API.M;
end

switch Obj.GLOBAL_SOFAConventions
%%
  case 'SimpleFreeFieldHRIR'
    figure; hold on;
    h=[];
    OE=SOFAexpand(Obj);
      % Listener: determine unique listener geometries
    LP=OE.ListenerPosition(index,:); 
    LV=OE.ListenerView(index,:); 
    LU=OE.ListenerUp(index,:);
    x=unique([LP LV LU],'rows');
      % Listener: plot first appearance
    h(end+1)=plot3(x(1,1), x(1,2), x(1,3),'ro','MarkerFaceColor',[1 0 0]);  % ListenerPosition
    h(end+1)=plot3(x(1,4), x(1,5), x(1,6),'ro','MarkerFaceColor',[1 1 1]); % ListenerView
    line([x(1,1), x(1,4)], [x(1,2) x(1,5)], [x(1,3) x(1,6)],'Color',[1 0 0]); % Line from LP to LV
      % Listener: plot other appearances
    for ii=2:size(x,1)
      plot3(x(ii,1), x(ii,2), x(ii,3),'ro','MarkerFaceColor',[1 0 0]);  % ListenerPosition
      plot3(x(ii,4), x(ii,5), x(ii,6),'ro','MarkerFaceColor',[1 1 1]); % ListenerView
      line([x(ii,1), x(ii,4)], [x(ii,2) x(ii,5)], [x(ii,3) x(ii,6)],'Color',[1 0 0]); % Line from LP to LV
    end
      % Plot receivers (at the moment we consider the first ListenerPosition only)
    if ndims(Obj.ReceiverPosition)>2
      colors='rbgymc';
      RP=shiftdim(Obj.ReceiverPosition(:,:,index),2);
      h(end+1)=plot3(LP(1,1)+RP(:,1,1), LP(1,2)+RP(:,1,2), LP(1,3)+RP(:,1,3),'rx');
      for ii=2:size(RP,2)
        plot3(LP(1,1)+RP(:,ii,1), LP(1,2)+RP(:,ii,2), LP(1,3)+RP(:,ii,3),[colors(ii) 'x']);
      end
    else % single receivers
      RP=Obj.ReceiverPosition;
      h(end+1)=plot3(LP(1,1)+RP(:,1), LP(1,2)+RP(:,2), LP(1,3)+RP(:,3),'rx');
    end
      % Plot source
    S=Obj.SourcePosition(index,:);
    if strcmp(Obj.SourcePosition_Type,'spherical');
      [X,Y,Z]=sph2cart(deg2rad(S(:,1)),deg2rad(S(:,2)),S(:,3));
    else
      X=S(:,1); Y=S(:,2); Z=S(:,3);
    end
    h(end+1)=plot3(X,Y,Z,'k.');
    legend(h,{'ListenerPosition','ListenerView','Receivers','SourcePosition'});
    xlabel(['X (in ' Obj.ListenerPosition_Units ')']);
    ylabel(['Y (in ' Obj.ListenerPosition_Units ')']);
    zlabel(['Z (in ' Obj.ListenerPosition_Units ')']);
    title('SimpleFreeFieldHRIR');
    
%%    
  case 'SingleRoomDRIR'    
    Obj=SOFAexpand(Obj);
    figure('Position',[1 1 (Obj.RoomCornerB(1)-Obj.RoomCornerA(1))*1.2 Obj.RoomCornerB(2)-Obj.RoomCornerA(2)]*100);
    axis([Obj.RoomCornerA(1)-0.5 Obj.RoomCornerB(1)+0.5 Obj.RoomCornerA(2)-0.5 Obj.RoomCornerB(2)+0.5]);
    box on; hold on; h=[];
      % plot the room
    rectangle('Position',[Obj.RoomCornerA(1) ...  
                          Obj.RoomCornerA(2) ...
                          Obj.RoomCornerB(1)-Obj.RoomCornerA(1) ...
                          Obj.RoomCornerB(2)-Obj.RoomCornerA(2)]);
      % plot listener and receivers
    LP=Obj.ListenerPosition(index,:); LV=Obj.ListenerView(index,:);
    RP=Obj.ReceiverPosition(index,:);
    for ii=1:size(LP,1)
      line([LP(ii,1), LV(ii,1)+LP(ii,1)], [LP(ii,2) LV(ii,2)+LP(ii,2)], 'Color',[1 0 0]);
    end
    h(end+1)=plot3(LP(:,1), LP(:,2), LP(:,3),'ro','MarkerFaceColor',[1 0 0]);
    h(end+1)=plot3(LV(:,1)+LP(:,1), LV(:,2)+LP(:,2), LV(:,3)+LP(:,3),'ro','MarkerFaceColor',[1 1 1]);
%     h(end+1)=plot3(LP(:,1)+RP(:,1), LP(:,2)+RP(:,2), ....
%                     LP(:,3)+RP(:,3),'rx');
    [lva,lve,lvr]=cart2sph(LV(ii,1),LV(ii,2),LV(ii,3));
    [rpa,rpe,rpr]=cart2sph(RP(ii,1),RP(ii,2),RP(ii,3));
    [x,y,z]=sph2cart(lva+rpa,lve+rpe,rpr);
    plot3(LP(ii,1)+x, LP(ii,2)+y, LP(ii,3)+z,'rx')
      % Plot source
    S=Obj.SourcePosition(index,:);    
    X=S(:,1); Y=S(:,2); Z=S(:,3);
    h(end+1)=plot3(X,Y,Z,'b.');
%     legend(h,{'ListenerPosition','ListenerView','Receivers','SourcePosition'});
    legend(h,{'ListenerPosition','ListenerView','SourcePosition'});
    title('SingleRoomDRIR');
  otherwise
    error('SOFAConventions not supported');
end
