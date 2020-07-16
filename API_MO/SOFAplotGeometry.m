function SOFAplotGeometry(Obj, index)
% SOFAplotGeometry(Obj) plots the geometry found in the Obj.
% 
% SOFAplotGeometry(Obj, index) plots the geometry for the measurements
% given in the index. 
%
% Supported conventions: 
%    SimpleFreeFieldHRIR'
%    SingleRoomDRIR
%    SimpleFreeFieldTF
%    FreeFieldDirectivityTF
%    some special cases of GeneralFIR.
%
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


if ~exist('index','var')
  index=1:Obj.API.M;  
elseif any(index > Obj.API.M)
	error("Index out of range. Only " + Obj.API.M + ...
        " measurement(s) performed.")
elseif any(index < 1)
    error("Choose index to be >= 1.")
end

switch Obj.GLOBAL_SOFAConventions
%%
  case {'SimpleFreeFieldHRIR','SingleRoomDRIR','SimpleFreeFieldTF','FreeFieldDirectivityTF'}
    % Expand entries to the same number of measurement points
    Obj = SOFAexpand(Obj);
    % See if the room geometry is specified
%     if strcmp(Obj.GLOBAL_RoomType,'shoebox')
%         figure('Position',[1 1 (Obj.RoomCornerB(1)-Obj.RoomCornerA(1))*1.2 Obj.RoomCornerB(2)-Obj.RoomCornerA(2)]*100);
%         box on; hold on;
%         % plot the room
%         rectangle('Position',[Obj.RoomCornerA(1) ...  
%                               Obj.RoomCornerA(2) ...
%                               Obj.RoomCornerB(1)-Obj.RoomCornerA(1) ...
%                               Obj.RoomCornerB(2)-Obj.RoomCornerA(2)]);
%     else
%         figure; hold on;
%     end
    
     if strcmp(Obj.GLOBAL_RoomType,'shoebox')
        x = min(Obj.RoomCornerA(1), Obj.RoomCornerB(1));
        xd = max(Obj.RoomCornerA(1), Obj.RoomCornerB(1));
        y = min(Obj.RoomCornerA(2), Obj.RoomCornerB(2));
        yd = max(Obj.RoomCornerA(2), Obj.RoomCornerB(2));
        w = xd - x;
        h = yd - y;
        figure('Position',[1 1 w*1.2 h]*100);
        box on; hold on;
        % plot the room
        rectangle('Position',[x y w h]);
    else
        figure; hold on;
    end
    
    legendEntries = [];
    title(sprintf('%s, %s',Obj.GLOBAL_SOFAConventions,Obj.GLOBAL_RoomType));
    % Get ListenerPosition, ListenerView, ReceiverPosition, and SourcePosition
    % NOTE: ListenerPosition is set to [0 0 0] for SimpleFreeFieldHRIR
    LP = SOFAconvertCoordinates(Obj.ListenerPosition(index,:),Obj.ListenerPosition_Type,'cartesian');
    if isfield(Obj,'ListenerView')
        LV = SOFAconvertCoordinates(Obj.ListenerView(index,:),Obj.ListenerView_Type,'cartesian');
    else
        LV = repmat([1 0 0],size(LP,1),1);
    end
    
    RP = SOFAconvertCoordinates(Obj.ReceiverPosition(:,:,index),Obj.ReceiverPosition_Type,'cartesian');

    
    if strcmp(Obj.SourcePosition_Type,'cartesian')
        SP=Obj.SourcePosition(index,:);
    else
        SP  = SOFAconvertCoordinates(Obj.SourcePosition(index,:),Obj.SourcePosition_Type,'cartesian');
    end
    if isfield(Obj,'SourceView')
        if strcmp(Obj.SourceView_Type,'cartesian')
            SV=Obj.SourceView(index,:);
        else
            SV  = SOFAconvertCoordinates(Obj.SourceView(index,:),Obj.SourceView_Type,'cartesian');
        end
    else
        SV = zeros(size(LP));   % if SourceView not given, source is omnidirectional
    end
    % Use only unique listener and source positons
    uniquePoints = unique([LP LV SP SV],'rows');
    LP = uniquePoints(:,1:3);
    LV = uniquePoints(:,4:6);
    SP = uniquePoints(:,7:9);
    SV = uniquePoints(:,10:12);
    % Plot ListenerPosition
    legendEntries(end+1) = plot3(LP(:,1),LP(:,2),LP(:,3),'ro','MarkerFaceColor',[1 0 0]);
    % Plot ListenerView
    for ii=1:size(LV,1)
      % Scale size of ListenerView vector smaller
      LV(ii,:) = 0.2*LV(ii,:)./norm(LV(ii,:));
      % Plot line for ListenerView vector
%       line([LP(ii,1), LV(ii,1)+LP(ii,1)], [LP(ii,2) LV(ii,2)+LP(ii,2)], 'Color',[1 0 0]);
      plot3([LP(ii,1) LV(ii,1)+LP(ii,1)],[LP(ii,2) LV(ii,2)+LP(ii,2)],[LP(ii,3) LV(ii,3)+LP(ii,3)], 'Color',[1 0 0]);

    end
    legendEntries(end+1) = plot3(LV(:,1)+LP(:,1),LV(:,2)+LP(:,2),LV(:,3)+LP(:,3),'ro','MarkerFaceColor',[1 1 1]);
    % Plot ReceiverPositon (this is plotted only for the first ListenerPosition)
    if ndims(RP)>2
        % If ReceiverPosition has more than two dimensions reduce it to the first
        % ListenerPosition
        RP = shiftdim(RP,2);
        RP = squeeze(RP(1,:,:));
        RP = reshape(RP,[size(Obj.ReceiverPosition,1), Obj.API.C]);
    end
    legendEntries(end+1) = plot3(LP(1,1)+RP(1,1), LP(1,2)+RP(1,2), LP(1,3)+RP(1,3),'rx');
    for ii=2:size(RP,1)
      plot3(LP(1,1)+RP(ii,1), LP(1,2)+RP(ii,2), LP(1,3)+RP(ii,3),'rx');
    end
    % Plot SourcePosition
    legendEntries(end+1)=plot3(SP(:,1),SP(:,2),SP(:,3),'ks','MarkerFaceColor',[0 0 0]);
    % Plot SourceView
    for ii=1:size(SV,1)
      % Scale size of SourceView vector smaller
      SV(ii,:) = 0.2*SV(ii,:)./norm(SV(ii,:));
      % Plot line for SourceView vector
%       line([SP(ii,1), SV(ii,1)+SP(ii,1)], [SP(ii,2) SV(ii,2)+SP(ii,2)], 'Color',[0 0 0]);
      plot3([LP(ii,1) LV(ii,1)+LP(ii,1)],[LP(ii,2) LV(ii,2)+LP(ii,2)],[LP(ii,3) LV(ii,3)+LP(ii,3)], 'Color',[1 0 0]);

    end
    legendEntries(end+1) = plot3(SV(:,1)+SP(:,1),SV(:,2)+SP(:,2),SV(:,3)+SP(:,3),'ks','MarkerFaceColor',[1 1 1]);
    % create legend
    legend(legendEntries,{'ListenerPosition','ListenerView','Receivers','SourcePosition','SourceView'},'Location','NorthEastOutside');
    xlabel(['X / ' Obj.ListenerPosition_Units]);
    ylabel(['Y / ' Obj.ListenerPosition_Units]);
    zlabel(['Z / ' Obj.ListenerPosition_Units]);

  case {'GeneralFIR'}
    % Expand entries to the same number of measurement points
    ObjC=Obj; % save as compact object
    Obj = SOFAexpand(Obj);
    % See if the room geometry is specified
%     if strcmp(Obj.GLOBAL_RoomType,'shoebox')
%         figure('Position',[1 1 (Obj.RoomCornerB(1)-Obj.RoomCornerA(1))*1.2 Obj.RoomCornerB(2)-Obj.RoomCornerA(2)]*100);
%         box on; hold on;
%         % plot the room
%         rectangle('Position',[Obj.RoomCornerA(1) ...  
%                               Obj.RoomCornerA(2) ...
%                               Obj.RoomCornerB(1)-Obj.RoomCornerA(1) ...
%                               Obj.RoomCornerB(2)-Obj.RoomCornerA(2)]);
%     else
%         figure; hold on;
%     end
    
    if strcmp(Obj.GLOBAL_RoomType,'shoebox')
        x = min(Obj.RoomCornerA(1), Obj.RoomCornerB(1));
        xd = max(Obj.RoomCornerA(1), Obj.RoomCornerB(1));
        y = min(Obj.RoomCornerA(2), Obj.RoomCornerB(2));
        yd = max(Obj.RoomCornerA(2), Obj.RoomCornerB(2));
        w = xd - x;
        h = yd - y;
        figure('Position',[1 1 w*1.2 h]*100);
        box on; hold on;
        % plot the room
        rectangle('Position',[x y w h]);
    else
        figure; hold on;
    end
    
    legendEntries = [];
    title(sprintf('%s, %s',Obj.GLOBAL_SOFAConventions,Obj.GLOBAL_RoomType));
    % Get ListenerPosition, ListenerView, ReceiverPosition, and SourcePosition
    % NOTE: ListenerPosition is set to [0 0 0] for SimpleFreeFieldHRIR
    LP = SOFAconvertCoordinates(Obj.ListenerPosition(index,:),Obj.ListenerPosition_Type,'cartesian');
    if isfield(Obj,'ListenerView')
        LV = SOFAconvertCoordinates(Obj.ListenerView(index,:),Obj.ListenerView_Type,'cartesian');
    else
        LV = repmat([1 0 0],size(LP,1),1);
    end
    if strcmp(Obj.SourcePosition_Type,'cartesian')
        SP=Obj.SourcePosition(index,:);
    else
        SP  = SOFAconvertCoordinates(Obj.SourcePosition(index,:),Obj.SourcePosition_Type,'cartesian');
    end
    % Use only unique listener and source positons
    uniquePoints = unique([LP LV SP],'rows');
    LP = uniquePoints(:,1:3);
    LV = uniquePoints(:,4:6);
    SP  = uniquePoints(:,7:9);
    % Plot ListenerPosition
    legendEntries(end+1) = plot3(LP(:,1),LP(:,2),LP(:,3),'ro','MarkerFaceColor',[1 0 0]);
    % Plot ListenerView
    for ii=1:size(LV,1)
      % Scale size of ListenerView vector smaller
      LV(ii,:) = 0.2*LV(ii,:)./norm(LV(ii,:));
      % Plot line for ListenerView vector
%       line([LP(ii,1), LV(ii,1)+LP(ii,1)], [LP(ii,2) LV(ii,2)+LP(ii,2)], 'Color',[1 0 0]);
      plot3([LP(ii,1) LV(ii,1)+LP(ii,1)],[LP(ii,2) LV(ii,2)+LP(ii,2)],[LP(ii,3) LV(ii,3)+LP(ii,3)], 'Color',[1 0 0]);
    end
    legendEntries(end+1) = plot3(LV(:,1)+LP(:,1),LV(:,2)+LP(:,2),LV(:,3)+LP(:,3),'ro','MarkerFaceColor',[1 1 1]);
    % Plot ReceiverPositon (this is plotted only for the first ListenerPosition)
    RP = SOFAconvertCoordinates(ObjC.ReceiverPosition(:,:,:),ObjC.ReceiverPosition_Type,'cartesian');
    if ndims(RP)>2
        % If ReceiverPosition has more than two dimensions reduce it to the first
        % ListenerPosition
        RP = shiftdim(RP,2);
        RP = unique(RP);
        RP = squeeze(RP(1,:,:));
    end
    legendEntries(end+1) = plot3(LP(1,1)+RP(1,1), LP(1,2)+RP(1,2), LP(1,3)+RP(1,3),'rx');
    for ii=2:size(RP,1)
      plot3(LP(1,1)+RP(ii,1), LP(1,2)+RP(ii,2), LP(1,3)+RP(ii,3),'rx');
    end
    % Plot SourcePosition
    legendEntries(end+1)=plot3(SP(:,1),SP(:,2),SP(:,3),'k.');
    legend(legendEntries,{'ListenerPosition','ListenerView','Receivers','SourcePosition'},'Location','NorthEastOutside');
    xlabel(['X / ' Obj.ListenerPosition_Units]);
    ylabel(['Y / ' Obj.ListenerPosition_Units]);
    zlabel(['Z / ' Obj.ListenerPosition_Units]);
    
  otherwise
    error('This SOFAConventions is not supported for plotting');
end

% Set fixed aspect ratio
axis equal;
% Add a little bit extra space at the axis
axisLimits = axis();
paddingSpace = 0.2 * max(abs(axisLimits(:)));
axisLimits([1 3]) = axisLimits([1 3]) - paddingSpace;
axisLimits([2 4]) = axisLimits([2 4]) + paddingSpace;
axis(axisLimits);
