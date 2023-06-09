function SOFAplotGeometry(Obj0,varargin)
%SOFAplotGeometry - Plot the geometry found a SOFA object
%   Usage: SOFAplotGeometry(Obj)
%
%   SOFAplotGeometry(Obj) plots the geometry found in a SOFA object
%   for all measurements, i.e., along the M dimension.
%   Obj needs to be one of the following conventions: 
%     SimpleFreeFieldHRIR
%     SimpleFreeFieldHRTF
%     SingleRoomDRIR
%     FreeFieldDirectivityTF
%   some special cases of GeneralFIR.
% 
%   SOFAplotGeometry(Obj, index) plots the geometry for the measurements
%   given in the index. 
%
%   SOFAplotGeometry(Obj,key,value) specifies in more detail the plotting:
%     'index'     : consider only specific measurements along the dimension M. 
%                   Value is a vector with entries from 1 to M. 
%     'SHorder'   : If the coordinate system is Spherical Harmonics, value defines
%                   the order of the spherical harmonic to be plotted. 
%                   Default: largest order available.
%     'SHm'       : If the coordinate system is Spherical Harmonics, value defines
%                   the degree of the shperical harmonic to be plotted.
%                   Default: largest degree available.
%
%   SOFAplotGeometry(Obj,key,value,flags) specifies in more even detail the plotting:
%     'normalize'  : Normalize the size of the View and Up vectors. 
%     'original'   : Plot the View and Up vectors as stored in Obj. 
%                    Default: normalize.
%
%    Multiple key-value pairs combined with the flags can be provided. 


% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: bug fixed when extracting LU (listener up) coordinates (28.12.2021)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

definput.keyvals.index=1:Obj0.API.M;
definput.keyvals.shorder=Inf;
definput.keyvals.shm=Inf;
definput.flags.normalize={'normalize','original'};
argin=varargin;
for ii=1:length(argin)
    if ischar(argin{ii}), argin{ii}=lower(argin{ii}); end
end
[flags,kv] = SOFAarghelper({'index','shorder','shm'},definput,argin);
index = kv.index;
SHorder=kv.shorder;
SHm=kv.shm;
flags.do_normalize = flags.normalize;

if any(index > Obj0.API.M)
    error(['Index out of range. Only ', num2str(Obj0.API.M), ...
         ' measurement(s) performed.'])
elseif any(index < 1)
    error('Choose index to be >= 1.')
end

switch Obj0.GLOBAL_SOFAConventions
%%
  case {'SimpleFreeFieldHRTF','SimpleFreeFieldHRIR','SingleRoomDRIR','FreeFieldDirectivityTF','GeneralFIR','GeneralTFE','FreeFieldHRIR','FreeFieldHRTF','GeneralTF-E'}
    % Expand entries to the same number of measurement points
    Obj = SOFAexpand(Obj0);
    % See if the room geometry is specified
     if strcmpi(Obj.GLOBAL_RoomType,'shoebox')
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
    % Get ListenerPosition, ReceiverPosition, SourcePosition, and
    % EmitterPosition
    % NOTE: ListenerPosition is set to [0 0 0] for SimpleFreeFieldHRIR
    LP = SOFAconvertCoordinates(Obj.ListenerPosition(index,:),Obj.ListenerPosition_Type,'cartesian');
    if ~(strcmpi(Obj.ReceiverPosition_Type,'Spherical Harmonics'))
        if size(Obj.ReceiverPosition,3)==1, idx=1; else idx=index; end      
        RP = SOFAconvertCoordinates(Obj.ReceiverPosition(:,:,idx),Obj.ReceiverPosition_Type,'cartesian');
    end
    if size(Obj.SourcePosition,1)==1, idx=1; else idx=index; end    
    SP = SOFAconvertCoordinates(Obj.SourcePosition(idx,:),Obj.SourcePosition_Type,'cartesian');
    if ~(strcmpi(Obj.EmitterPosition_Type,'Spherical Harmonics'))
        if size(Obj.EmitterPosition,3)==1, idx=1; else idx=index; end
        EP = SOFAconvertCoordinates(Obj.EmitterPosition(:,:,idx),Obj.EmitterPosition_Type,'cartesian');
    end
    if isfield(Obj,'ListenerView')
        if size(Obj.ListenerView,1)==1, idx=1; else idx=index; end      
        LV = SOFAconvertCoordinates(Obj.ListenerView(idx,:),Obj.ListenerView_Type,'cartesian');
    end
    if isfield(Obj,'ListenerUp')
        try 
            if size(Obj.ListenerUp,1)==1, idx=1; else idx=index; end          
            LU = SOFAconvertCoordinates(Obj.ListenerUp(idx,:),Obj.ListenerUp_Type,'cartesian');
        catch
            % if listerUp_type is not defined try using listenerView_type
            % instead
            if size(Obj.ListenerUp,1)==1, idx=1; else idx=index; end        
            LU = SOFAconvertCoordinates(Obj.ListenerUp(idx,:),Obj.ListenerView_Type,'cartesian');
        end
    end   
    if isfield(Obj,'SourceView')
      if size(Obj.SourceView,1)==1, idx=1; else idx=index; end      
    	SV  = SOFAconvertCoordinates(Obj.SourceView(idx,:),Obj.SourceView_Type,'cartesian');
    end
    if isfield(Obj,'SourceUp')
        try 
            if size(Obj.SourceUp,1)==1, idx=1; else idx=index; end
            SU = SOFAconvertCoordinates(Obj.SourceUp(idx,:),Obj.SourceUp_Type,'cartesian');
        catch
            if size(Obj.SourceUp,1)==1, idx=1; else idx=index; end      
            SU = SOFAconvertCoordinates(Obj.SourceUp(idx,:),Obj.SourceView_Type,'cartesian');
        end
    end
    % Use only unique listener and source positons
    caseString = '';
    uniquePoints = [LP SP];
    if exist('LV')
        uniquePoints = [uniquePoints LV];
        caseString = strcat(caseString , 'LV');
    end
    if exist('LU')
        uniquePoints = [uniquePoints LU];
        caseString = strcat(caseString, 'LU');
    end
    if exist('SV')
        uniquePoints = [uniquePoints SV];
        caseString = strcat(caseString, 'SV');
    end
    if exist('SU')
        uniquePoints = [uniquePoints SU];
        caseString = strcat(caseString, 'SU');
    end
    
    uniquePoints = unique(uniquePoints,'rows');
    switch caseString
        case ''
        	LP = uniquePoints(:,1:3);
            SP = uniquePoints(:,4:6);
        case 'LV'
        	LP = uniquePoints(:,1:3);
            SP = uniquePoints(:,4:6);
            LV = uniquePoints(:,7:9);
        case 'LVLU'
        	LP = uniquePoints(:,1:3);
            SP = uniquePoints(:,4:6);
            LV = uniquePoints(:,7:9);
%             LU = uniquePoints(:,7:9); % I think this was a bug (miho)
            LU = uniquePoints(:,10:12);
        case 'LVLUSV'
        	LP = uniquePoints(:,1:3);
            SP = uniquePoints(:,4:6);
            LV = uniquePoints(:,7:9);
            LU = uniquePoints(:,10:12);
            SV = uniquePoints(:,13:15);
        case 'SV'
        	LP = uniquePoints(:,1:3);
            SP = uniquePoints(:,4:6);
            SV = uniquePoints(:,7:9);
        case 'SVSU'
        	LP = uniquePoints(:,1:3);
            SP = uniquePoints(:,4:6);
            SV = uniquePoints(:,7:9);
            SU = uniquePoints(:,10:12);
        case 'LVSV'
        	LP = uniquePoints(:,1:3);
            SP = uniquePoints(:,4:6);
            LV = uniquePoints(:,7:9);
            SV = uniquePoints(:,10:12);
        case 'LVSVSU'
        	LP = uniquePoints(:,1:3);
            SP = uniquePoints(:,4:6);
            LV = uniquePoints(:,7:9);
            SV = uniquePoints(:,10:12);
            SU = uniquePoints(:,13:15);
        case 'LVLUSVSU'
        	LP = uniquePoints(:,1:3);
            SP = uniquePoints(:,4:6);
            LV = uniquePoints(:,7:9);
            LU = uniquePoints(:,10:12);
            SV = uniquePoints(:,13:15);
            SU = uniquePoints(:,16:18);
        otherwise
            error('This SOFAConventions is not supported for plotting');
    end

    % Plot ListenerPosition
    legendEntries(end+1) = plot3(LP(:,1),LP(:,2),LP(:,3),'ro','MarkerFaceColor','r','MarkerSize',5);
    if strcmpi(Obj.ReceiverPosition_Type,'Spherical Harmonics')
        maxSHorder = sqrt(Obj.API.R)-1;
         % set SHorder to max if user didn't specify it
        if isinf(SHorder)
            SHorder = maxSHorder;
        end
        % check if chosen SHorder is possible
        if SHorder > maxSHorder
            error(['Chosen SHorder not possibile, only orders up to ', ...
                num2str(maxSHorder), ' possible.'])
        elseif SHorder < 0
            error('Chosen SHorder not possibile, as it must be positive.')
        end
        x0 = Obj.ListenerPosition(1,1);
        y0 = Obj.ListenerPosition(1,2);
        z0 = Obj.ListenerPosition(1,3);
        
        % check for m given by the user and if it is possible
        if isinf(SHm)
            % if not set to some value
            SHm = -floor(1/2 * SHorder);
        elseif abs(SHm) > SHorder
               error(['Chosen SHm not possibile, must be in range of abs(', ...
                num2str(SHorder), ').'])
        end
        % if possibile set SHmForPlotting
        SHmForPlotting = power(SHorder,2)+SHorder+SHm+1;
        
        [X,Y,Z] = sphere(50); 
        [azi_rad,elev_rad,~] = cart2sph(X,Y,Z);
        azi_length =size(azi_rad,1);
        elev_length=size(elev_rad,1);
        azi= azi_rad/pi*180;
        elev = elev_rad/pi*180;
        azi = azi(:);
        elev = elev(:);
        
        S = sph2SH([azi,elev], SHorder);
        S = S(:,SHmForPlotting);
        S = reshape(S,[azi_length,elev_length]);
        
        r_sphere = 0.7*max(max(S))*randi(2,size(S)); 
        r = abs(S) + r_sphere;        
        
        [D_x,D_y,D_z] = sph2cart(azi_rad,elev_rad,abs(r));
        legendEntries(end+1) = surf(D_x+x0,D_y+y0,D_z+z0,Y,'LineStyle','none','FaceAlpha',0.09);
%     elseif strcmpi(Obj.ReceiverPosition_Type,'spherical')
%         S = sqrt(Obj.API.R-1);
%         x0 = Obj.ListenerPosition(1,1);
%         y0 = Obj.ListenerPosition(1,2);
%         theta = -pi : 0.01 : pi;
%         r = 1;
%         phi = sin(S*theta);
%         phi_negativ = sin(-S*theta);
%         
%         [x,y] = pol2cart(theta,(r*(1+ abs(phi)+ abs(phi_negativ)))./3);
%         legendEntries(end+1)=plot(x+x0,y+y0,'LineStyle','--','Color',[0.741 0.747 0.741]);
% 
% %         text(x0,y0+r,['Order: ',num2str(S)],'HorizontalAlignment',...
% %            'center','VerticalAlignment','bottom')

    else
       % Plot ReceiverPositon (this is plotted only for the first ListenerPosition)
        if ndims(RP)>2
            % If ReceiverPosition has more than two dimensions reduce it to the first
            % ListenerPosition
            RP = shiftdim(RP,2);
            RP = squeeze(RP(1,:,:));
            RP = reshape(RP,[size(Obj.ReceiverPosition,1), Obj.API.C]);
        end
        legendEntries(end+1) = plot3(LP(1,1)+RP(1,1), LP(1,2)+RP(1,2), LP(1,3)+RP(1,3),'r*','MarkerSize',8);
        for ii=2:size(RP,1)
          plot3(LP(1,1)+RP(ii,1), LP(1,2)+RP(ii,2), LP(1,3)+RP(ii,3),'r*','MarkerSize',8);
        end
    end
    % Plot SourcePosition
    legendEntries(end+1)=plot3(SP(:,1),SP(:,2),SP(:,3),'b.','MarkerSize',7);
    % Plot EmitterPositions depending on Type
    if strcmpi(Obj.EmitterPosition_Type,'Spherical Harmonics')
        maxSHorder = sqrt(Obj.API.E)-1;
        % set SHorder to max if user didn't specify it
        if isinf(SHorder)
            SHorder = maxSHorder;
        end
        % check if chosen SHorder is possible
        if SHorder > maxSHorder
            error(['Chosen SHorder not possibile, only orders up to ', ...
                num2str(maxSHorder), ' possible.'])
        elseif SHorder < 0
            error('Chosen SHorder not possibile, as it must be positive.')
        end
        x0 = Obj.SourcePosition(1,1);
        y0 = Obj.SourcePosition(1,2);
        z0 = Obj.SourcePosition(1,3);
        
        % check for m given by the user
        if isinf(SHm)
            SHm = -floor(1/2 * SHorder);
        elseif abs(SHm) > SHorder
               error(['Chosen SHm not possibile, must be in range of abs(', ...
                num2str(SHorder), ').'])
        end
        % if possibile set SHmForPlotting
        SHmForPlotting = power(SHorder,2)+SHorder+SHm+1;
        
        [X,Y,Z] = sphere(50); 
        [azi_rad,elev_rad,~] = cart2sph(X,Y,Z);
        azi_length =size(azi_rad,1);
        elev_length=size(elev_rad,1);
        azi= azi_rad/pi*180;
        elev = elev_rad/pi*180;
        azi = azi(:);
        elev = elev(:);
        
        S = sph2SH([azi,elev], SHorder);
        S = S(:,SHmForPlotting);
        S = reshape(S,[azi_length,elev_length]);
        
        r_sphere = 0.7*max(max(S))*randi(2,size(S)); 
        r = abs(S) + r_sphere;        
        
        [D_x,D_y,D_z] = sph2cart(azi_rad,elev_rad,abs(r));
        legendEntries(end+1) = surf(D_x+x0,D_y+y0,D_z+z0,Y,'LineStyle','none','FaceAlpha',0.09);

%     elseif strcmpi(Obj.EmitterPosition_Type,'spherical')
%         S = sqrt(Obj.API.R-1);
%         x0 = Obj.SourcePosition(1,1);
%         y0 = Obj.SourcePosition(1,2);
%         theta = -pi : 0.01 : pi;
%         r = 1;
%         phi = sin(S*theta);
%         phi_negativ = sin(-S*theta);
%         
%         [x,y] = pol2cart(theta,(r*(1+ abs(phi)+ abs(phi_negativ)))./3);
%         legendEntries(end+1)=plot(x+x0,y+y0,'LineStyle','--','Color',[0.741 0.747 0.741]);
% 
% %         text(x0,y0+r,['Order: ',num2str(S)],'HorizontalAlignment',...
% %            'center','VerticalAlignment','bottom')

    else
        % Plot EmitterPosition
        if ~isequal(Obj0.EmitterPosition,[0 0 0]) % plot only if not simple emitter in the source's center
          if ndims(EP)>2
              % If EmitterPosition has more than two dimensions reduce it to the first
              % ListenerPosition
              EP = shiftdim(EP,2);
              EP = squeeze(EP(1,:,:));
              EP = reshape(EP,[size(Obj.EmitterPosition,1), Obj.API.C]);
          end
          % plot Emitters for first Source
          legendEntries(end+1) = plot3(SP(1,1)+EP(1,1), SP(1,2)+EP(1,2), SP(1,3)+EP(1,3),'b+','MarkerSize',8);
          for ii=2:size(EP,1)
              plot3(SP(1,1)+EP(ii,1), SP(1,2)+EP(ii,2), SP(1,3)+EP(ii,3),'b+','MarkerSize',8);
          end
          % plot all Emitters for each Source
          for jj=2:size(SP,1)
              for ii=1:size(EP,1)
                plot3(SP(jj,1)+EP(ii,1), SP(jj,2)+EP(ii,2), SP(jj,3)+EP(ii,3),'b+');
              end
          end
        end
    end
    if exist('LV','var')
        % Plot ListenerView
        LV=unique(LV,'rows');
        for ii = 2:size(LV,1)
            % Scale size of ListenerView vector smaller
            if flags.do_normalize
                LV(ii,:) = LV(ii,:)./norm(LV(ii,:));
            end
            % Plot line for ListenerView vector
            quiver3(LP(ii,1),LP(ii,2),LP(ii,3),LV(ii,1),LV(ii,2),LV(ii,3),'Color',[1 0 0],'MarkerFaceColor',[1 0 0]);
        end
        if flags.do_normalize
            LV(1,:) = LV(1,:)./norm(LV(1,:));
        end
        legendEntries(end+1) = quiver3(LP(1,1),LP(1,2),LP(1,3),LV(1,1),LV(1,2),LV(1,3),'Color',[1 0 0],'MarkerFaceColor',[1 0 0]);
    end
    if exist('LU','var')
        LU=unique(LU,'rows');
        for ii = 2:size(LU,1)
            if flags.do_normalize
                LU(ii,:) = LU(ii,:)./norm(LU(ii,:));
            end
            quiver3(LP(ii,1),LP(ii,2),LP(ii,3),LU(ii,1),LU(ii,2),LU(ii,3),0,'AutoScale','off','Color',[0 0 0],'MarkerFaceColor',[0 0 0]);
%               quiver3(LP(ii,1),LP(ii,2),LP(ii,3),LU(ii,1),LU(ii,2),LU(ii,3),'Color',[0 0 0],'MarkerFaceColor',[0 0 0]);
%             quiver3(LP(ii,1),LP(ii,2),LP(ii,3),LV(ii,1),LV(ii,2),LV(ii,3),'Color',[1 0 0],'MarkerFaceColor',[1 0 0]);
        end
        if flags.do_normalize
            LU(1,:) = LU(1,:)./norm(LU(1,:));
        end
        legendEntries(end+1) = quiver3(LP(1,1),LP(1,2),LP(1,3),LU(1,1),LU(1,2),LU(1,3),0,'AutoScale','off','Color',[0 0 0],'MarkerFaceColor',[0 0 0]);
%         legendEntries(end+1) = quiver3(LP(1,1),LP(1,2),LP(1,3),LU(1,1),LU(1,2),LU(1,3),'Color',[0 0 0],'MarkerFaceColor',[0 0 0]);
%         legendEntries(end+1) = quiver3(LP(1,1),LP(1,2),LP(1,3),LV(1,1),LV(1,2),LV(1,3),'Color',[1 0 0],'MarkerFaceColor',[1 0 0]);
    end
    if exist('SV','var')
        SV=unique(SV,'rows');
        % Plot ListenerView
        for ii = 2:size(SV,1)
            % Scale size of ListenerView vector smaller
            if flags.do_normalize
                SV(ii,:) = SV(ii,:)./norm(SV(ii,:));
            end
            % Plot line for ListenerView vector
            quiver3(SP(ii,1),SP(ii,2),SP(ii,3),SV(ii,1),SV(ii,2),SV(ii,3),0,...
                'AutoScale','off',...
                'Color',[0 0 1],'MarkerFaceColor',[0 0 1]);
        end
        if flags.do_normalize
            SV(1,:) = SV(1,:)./norm(SV(1,:));
        end
        legendEntries(end+1) = quiver3(SP(1,1),SP(1,2),SP(1,3),SV(1,1),SV(1,2),SV(1,3),0,...
                'AutoScale','off',...
                'Color',[0 0 1],'MarkerFaceColor',[0 0 1]);
    end
    if exist('SU','var')
        SU=unique(SU,'rows');
        for ii = 2:size(SU,1)
            if flags.do_normalize
                SU(ii,:) = SU(ii,:)./norm(SU(ii,:));
            end
            quiver3(SP(ii,1),SP(ii,2),SP(ii,3),SU(ii,1),SU(ii,2),SU(ii,3),0,...
                'AutoScale','off',...
                'Color',[0 0 0],'MarkerFaceColor',[0 0 0]);
        end
        if flags.do_normalize
            SU(1,:) = SU(1,:)./norm(SU(1,:));
        end
        legendEntries(end+1) = quiver3(SP(1,1),SP(1,2),SP(1,3),SU(1,1),SU(1,2),SU(1,3),'Color',[0 0 0],'MarkerFaceColor',[0 0 0]);
    end
    % create legend
    legendDescription = {'ListenerPosition'};
    if (strcmpi(Obj.ReceiverPosition_Type,'Spherical Harmonics'))
        legendDescription{end+1} = ['Receiver (order: ', num2str(S_R) ,')'];
    else
        legendDescription{end+1} = 'ReceiverPosition';
    end
    legendDescription{end+1} ='SourcePosition';
    if ~isequal(Obj0.EmitterPosition,[0 0 0])
      if (strcmpi(Obj.EmitterPosition_Type,'Spherical Harmonics'))
          legendDescription{end+1} = ['Emitter (order: ', num2str(SHorder),', m: ', num2str(SHm),')'];
      else
          legendDescription{end+1} = 'EmitterPosition';
      end
    end
    if exist('LV','var')
        legendDescription{end+1} = 'ListenerView';
    end
    if exist('LU','var')
        legendDescription{end+1} = 'ListenerUp';
    end
    if exist('SV','var')
        legendDescription{end+1} = 'SourceView';
    end
	if exist('SU','var')
        legendDescription{end+1} = 'SourceUp';
	end
    legend(legendEntries,legendDescription,'Location','NorthEastOutside');
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
