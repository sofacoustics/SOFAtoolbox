function [Obj,modified] = SOFAupgradeConventions(Obj)
%SOFAupgradeConventions - Upgrade conventions to a higher version
%   Usage: Obj = SOFAupgradeConventions(Obj)
%
%   Obj = SOFAupgradeConventions(Obj) searches for a new version of the 
%   convention stored in Obj and upgrades to a next higher
%   version if found. 
%
%   [Obj,MODIFIED] = SOFAupgradeConventions(..) returns MODIFIED = 1 when 
%   an upgrade was performed. In order to obtain the most recent version, 
%   SOFAupgradeConventions can be processed recursively until MODIFIED is 0. 

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
%
% SOFA Toolbox - function SOFAupgradeConventions
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

modified=0;

%% Upgrade specific to a SOFA version

switch Obj.GLOBAL_Version,
  case '0.3'
    modified=1;
    % in SOFA 0.3, only SimpleFreeFieldHRIR 0.1 was supported.
    % Updating SimpleFreeFieldHRIR 0.1 to 0.2
    Obj.GLOBAL_Version='0.4';
    Obj.GLOBAL_SOFAConventionsVersion='0.2';
    Obj.GLOBAL_TimeCreated=Obj.GLOBAL_DatabaseTimeCreated;
    Obj.GLOBAL_TimeModified=Obj.GLOBAL_DatabaseTimeModified;
    Obj.GLOBAL_History=SOFAappendText(Obj,'GLOBAL_History','Upgraded from SOFA 0.3');
      % remove dimensional variables and not used variables/attributes
    dims={'I','R','E','N','M','C','Q','SourceView',...
      'SourceUp','GLOBAL_DatabaseTimeCreated','GLOBAL_DatabaseTimeModified'}; 
    f=fieldnames(Obj);
    for ii=1:length(dims)
      for jj=1:length(f)
        if strcmp(f{jj},dims{ii}), 
          Obj=rmfield(Obj,f{jj});   % remove variable or attribute
          if isempty(strfind(f{jj},'_')),
            Obj.API.Dimensions=rmfield(Obj.API.Dimensions,f{jj}); % remove dimension
          end
        elseif strcmp(f{jj}(1:min(length(dims{ii})+1,length(f{jj}))),[dims{ii} '_']) 
          Obj=rmfield(Obj,f{jj});  % remove attributes of that variable
        end
      end
    end
    warning('SOFA:upgrade','SOFA 0.3 upgraded to 0.4.   Use warning(''off'',''SOFA:upgrade''); to switch off this warning. ');
  case '0.4'
    % in SOFA 0.4, only SimpleFreeFieldHRIR might need upgrade
    if strcmp(Obj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR')
      switch Obj.GLOBAL_SOFAConventionsVersion
        case '0.2'
          % Upgrade from SimpleFreeFieldHRIR 0.2 to 0.3
          Obj.GLOBAL_History=SOFAappendText(Obj,'GLOBAL_History','Upgraded from SimpleFreeFieldHRIR 0.2');
          % Create temp SourcePosition
          azi=bsxfun(@times,Obj.ListenerRotation(:,1),ones(size(Obj.ListenerPosition,1),1));
          ele=bsxfun(@times,Obj.ListenerRotation(:,2),ones(size(Obj.ListenerPosition,1),1));
          r=bsxfun(@times,Obj.ListenerPosition(:,1),ones(size(Obj.ListenerRotation,1),1));
          % Copy ListenerPosition
          Obj.ListenerPosition=Obj.SourcePosition;
          % Overwrite SourcePosition
          Obj.SourcePosition=[azi ele r];
          Obj.SourcePosition_Type='spherical';
          Obj.SourcePosition_Units='degree, degree, meter';
          % Mirror the ListenerView and correct ListenerUp
          Obj.ListenerView=-Obj.ListenerView;
          Obj.ListenerUp=[0 0 1];
          % Remove irrelevant fields
          if isfield(Obj,'SourceView'); Obj=rmfield(Obj,'SourceView'); end
          if isfield(Obj,'SourceView_Type'); Obj=rmfield(Obj,'SourceView_Type'); end
          if isfield(Obj,'SourceView_Units'); Obj=rmfield(Obj,'SourceView_Units'); end
          if isfield(Obj,'SourceUp'); Obj=rmfield(Obj,'SourceUp'); end
          if isfield(Obj,'SourceUp_Type'); Obj=rmfield(Obj,'SourceUp_Type'); end
          if isfield(Obj,'SourceUp_Units'); Obj=rmfield(Obj,'SourceUp_Units'); end
          Obj=rmfield(Obj,'ListenerRotation');
          Obj=rmfield(Obj,'ListenerRotation_Type');
          Obj=rmfield(Obj,'ListenerRotation_Units');
          Obj.API.Dimensions=rmfield(Obj.API.Dimensions,'ListenerRotation');
          Obj.GLOBAL_SOFAConventionsVersion='0.3';
      end
    end    
		modified=1;          
		Obj.GLOBAL_Version='0.5';
    warning('SOFA:upgrade','SOFA 0.4 upgraded to 0.5.   Use warning(''off'',''SOFA:upgrade''); to switch off this warning. ');
  case '0.5'
		% Upgrade from 0.5 to 0.6
    Obj.GLOBAL_DateCreated=Obj.GLOBAL_TimeCreated;
    Obj=rmfield(Obj,'GLOBAL_TimeCreated');
    Obj.GLOBAL_DateModified=Obj.GLOBAL_TimeModified;
    Obj=rmfield(Obj,'GLOBAL_TimeModified');
    Obj.GLOBAL_Origin=Obj.GLOBAL_Source;
    Obj=rmfield(Obj,'GLOBAL_Source');
    if isfield(Obj,'ListenerView') && ~isfield(Obj,'ListenerView_Type')
      Obj.ListenerView_Type = 'cartesian';
      Obj.ListenerView_Units = 'meter';
    end
    if isfield(Obj,'ReceiverView') && ~isfield(Obj,'ReceiverView_Type')
      Obj.ReceiverView_Type = 'cartesian';
      Obj.ReceiverView_Units = 'meter';
    end
    if isfield(Obj,'GLOBAL_SubjectID'),
      Obj.GLOBAL_ListenerShortName=Obj.GLOBAL_SubjectID;	% rename SubjectID to ListenerShortName
      Obj=rmfield(Obj,'GLOBAL_SubjectID');
    end
    switch Obj.GLOBAL_SOFAConventions
      case {'SimpleFreeFieldHRIR', 'SimpleFreeFieldTF'}
        Obj.GLOBAL_SOFAConventionsVersion='0.4';
      case {'GeneralFIR', 'GeneralTF', 'SingleRoomDRIR'}
        Obj.GLOBAL_SOFAConventionsVersion='0.2';
    end   
    Obj.GLOBAL_History=SOFAappendText(Obj,'GLOBAL_History','Upgraded from SOFA 0.5');
    Obj.GLOBAL_Version='0.6';
    modified=1;
    warning('SOFA:upgrade','SOFA 0.5 upgraded to 0.6.   Use warning(''off'',''SOFA:upgrade''); to switch off this warning. ');
  case '0.6'
    X=SOFAgetConventions(Obj.GLOBAL_SOFAConventions);
    if ~isempty(X),
      Obj.GLOBAL_History=SOFAappendText(Obj,'GLOBAL_History','Upgraded from SOFA 0.6');
      Obj.GLOBAL_Version='1.0';
      Obj.GLOBAL_SOFAConventionsVersion = X.GLOBAL_SOFAConventionsVersion;
        % replace aliases by correct unit names
      U=SOFAdefinitions('units');
      Uf=fieldnames(U);
      f=fieldnames(Obj);
      for jj=1:length(f)
        if length(f{jj}) > 6
          if strcmp(f{jj}(end-5:end),'_Units')
            for ii=1:length(Uf) % _Units found, check for alias
              Obj.(f{jj})=regexprep(Obj.(f{jj}), U.(Uf{ii}), Uf{ii}, 'ignorecase');
            end
          end
        end
      end
      f=fieldnames(Obj.Data);
      for jj=1:length(f)
        if length(f{jj}) > 6
          if strcmp(f{jj}(end-5:end),'_Units')
            for ii=1:length(Uf) % _Units found, check for alias
              Obj.Data.(f{jj})=regexprep(Obj.Data.(f{jj}), U.(Uf{ii}), Uf{ii}, 'ignorecase');
            end
          end
        end
      end
      modified=1;
      warning('SOFA:upgrade','SOFA 0.6 upgraded to 1.0.   Use warning(''off'',''SOFA:upgrade''); to switch off this warning. ');    
    else
      warning('SOFA:upgrade','Unknown conventions');
    end
end

%% Upgrade specific to conventions
if ~modified
    switch Obj.GLOBAL_SOFAConventions
    case 'MultiSpeakerBRIR'
        %% Get an empy conventions structure
        ObjNew=SOFAgetConventions('SingleRoomMIMOSRIR');
        
        %% Transfer global objects
        ObjNew.GLOBAL_ListenerShortName = Obj.GLOBAL_ListenerShortName;          
        ObjNew.GLOBAL_ApplicationName = Obj.GLOBAL_ApplicationName;
        ObjNew.GLOBAL_ApplicationVersion = Obj.GLOBAL_ApplicationVersion;
        ObjNew.GLOBAL_AuthorContact = Obj.GLOBAL_AuthorContact;
        ObjNew.GLOBAL_Comment = Obj.GLOBAL_Comment;
        % ObjNew.GLOBAL_DataType % always FIR-E
        ObjNew.GLOBAL_History = SOFAappendText(Obj, 'GLOBAL_History', 'Converted from SOFAconvertMultiSpeakerBRIR.'); % append information
        ObjNew.GLOBAL_License = Obj.GLOBAL_License;
        ObjNew.GLOBAL_Organization = Obj.GLOBAL_Organization;
        ObjNew.GLOBAL_References = Obj.GLOBAL_References;
        ObjNew.GLOBAL_RoomType = Obj.GLOBAL_RoomType; % shoebox or dae
        ObjNew.GLOBAL_Origin = Obj.GLOBAL_Origin;
        ObjNew.GLOBAL_Title = Obj.GLOBAL_Title;
        ObjNew.GLOBAL_DatabaseName = Obj.GLOBAL_DatabaseName;
        ObjNew.GLOBAL_RoomDescription = Obj.GLOBAL_RoomDescription;
        
        %% Transfer data
        ObjNew.Data=Obj.Data; 
        ObjNew.Data.IR=permute(ObjNew.Data.IR, [1 2 4 3]); % permute because dimension ordner has been changed
        ObjNew.Data.Delay=repmat(ObjNew.Data.Delay,size(ObjNew.Data.IR,1),1);
        
        %% Upgrade some objects
        ObjNew.ListenerPosition=repmat(Obj.ListenerPosition(1,:),size(ObjNew.Data.IR,1),1);
        ObjNew.SourcePosition=repmat(Obj.SourcePosition(1,:),size(ObjNew.Data.IR,1),1);
        
        %% Transfer other objects
        ObjNew.ListenerPosition_Type=Obj.ListenerPosition_Type;
        ObjNew.ListenerPosition_Units=Obj.ListenerPosition_Units;
        ObjNew.ReceiverPosition=Obj.ReceiverPosition;
        ObjNew.ReceiverPosition_Type=Obj.ReceiverPosition_Type;
        ObjNew.ReceiverPosition_Units=Obj.ReceiverPosition_Units;
        ObjNew.SourcePosition_Type=Obj.SourcePosition_Type;
        ObjNew.SourcePosition_Units=Obj.SourcePosition_Units;
        ObjNew.EmitterPosition=Obj.EmitterPosition;
        ObjNew.EmitterPosition_Type=Obj.EmitterPosition_Type;
        ObjNew.EmitterPosition_Units=Obj.EmitterPosition_Units;
        ObjNew.ListenerUp=Obj.ListenerUp;
        ObjNew.ListenerView=Obj.ListenerView;
        ObjNew.ListenerView_Type=Obj.ListenerView_Type;
        % ObjNew.SourceUp=Obj.SourceUp;
        % ObjNew.SourceView=Obj.SourceView;
        % ObjNew.SourceView_Type=Obj.SourceView_Type;
        if isfield(Obj,'RoomCornerA')
            if size(Obj.RoomCornerA,1) == 3
                ObjNew.RoomCornerA=Obj.RoomCornerA'; % must be transposed for the data existing so far
            else
                ObjNew.RoomCornerA=Obj.RoomCornerA;
            end
            ObjNew.RoomCornerA_Type=Obj.RoomCornerA_Type;
            ObjNew.RoomCornerA_Units=Obj.RoomCornerA_Units;
        end
        
        if isfield(Obj,'RoomCornerB')
            if size(Obj.RoomCornerB,1) == 3
                ObjNew.RoomCornerB=Obj.RoomCornerB'; % must be transposed for the data existing so far
            else
                ObjNew.RoomCornerB=Obj.RoomCornerB;
            end
            ObjNew.RoomCornerB_Type=Obj.RoomCornerB_Type;
            ObjNew.RoomCornerB_Units=Obj.RoomCornerB_Units;
        end
        
        %% Update dimensions
        Obj=SOFAupdateDimensions(ObjNew);
        modified=1;
        warning('SOFA:upgrade','Conventions MultiSpeakerBRIR upgraded to SingleRoomMIMOSRIR.  Use warning(''off'',''SOFA:upgrade''); to switch off this warning. ');

%       if strcmp(Obj.GLOBAL_SOFAConventionsVersion,'0.1');
%           % upgrade to 0.2
%         Obj.GLOBAL_DataType='FIRE';
%         Obj.GLOBAL_SOFAConventionsVersion='0.2';
%         %Obj.Data.Delay = 
%         if strcmp(Obj.API.Dimensions.Data.Delay,'IR')
%           Obj.API.Dimensions.Data.Delay='IRE'; 
%           Obj.Data.Delay=repmat(Obj.Data.Delay,[1 1 size(Obj.EmitterPosition,1)]);
%         end
%         if strcmp(Obj.API.Dimensions.Data.Delay,'MR')
%           Obj.API.Dimensions.Data.Delay='MRE'; 
%           Obj.Data.Delay=repmat(Obj.Data.Delay,[1 1 size(Obj.EmitterPosition,1)]);
%         end
%         modified=1;
%         warning('SOFA:upgrade','Conventions MultiSpeakerBRIR 0.1 upgraded to 0.2.  Use warning(''off'',''SOFA:upgrade''); to switch off this warning. ');
%       end
    case 'SingleRoomDRIR'
        %% Get an empy conventions SingleRoomSRIR structure to upgrade to
        ObjNew=SOFAgetConventions('SingleRoomSRIR');
        
        %% Transfer global objects
        ObjNew.GLOBAL_ListenerShortName = Obj.GLOBAL_ListenerShortName;          
        ObjNew.GLOBAL_ApplicationName = Obj.GLOBAL_ApplicationName;
        ObjNew.GLOBAL_ApplicationVersion = Obj.GLOBAL_ApplicationVersion;
        ObjNew.GLOBAL_AuthorContact = Obj.GLOBAL_AuthorContact;
        ObjNew.GLOBAL_Comment = Obj.GLOBAL_Comment;
        % ObjNew.GLOBAL_DataType % always FIR
        ObjNew.GLOBAL_History = SOFAappendText(Obj, 'GLOBAL_History', 'Converted from SingleRoomDRIR.'); % append information
        ObjNew.GLOBAL_License = Obj.GLOBAL_License;
        ObjNew.GLOBAL_Organization = Obj.GLOBAL_Organization;
        ObjNew.GLOBAL_References = Obj.GLOBAL_References;
        ObjNew.GLOBAL_RoomType = Obj.GLOBAL_RoomType; % shoebox or dae
        ObjNew.GLOBAL_Origin = Obj.GLOBAL_Origin;
        ObjNew.GLOBAL_Title = Obj.GLOBAL_Title;
        ObjNew.GLOBAL_DatabaseName = Obj.GLOBAL_DatabaseName;
        ObjNew.GLOBAL_RoomDescription = Obj.GLOBAL_RoomDescription;
        
        %% Transfer data
        ObjNew.Data=Obj.Data;
        
        %% Upgrade some objects
        ObjNew.ListenerPosition=repmat(Obj.ListenerPosition(1,:),size(ObjNew.Data.IR,1),1);
        % for ii=1:size(ObjNew.Data.IR,1)
        %     ObjNew.ListenerPosition(ii,:)=Obj.ListenerPosition(1,:);
        % end
        
        %% Transfer other objects
        ObjNew.ListenerPosition_Type=Obj.ListenerPosition_Type;
        ObjNew.ListenerPosition_Units=Obj.ListenerPosition_Units;
        ObjNew.ReceiverPosition=Obj.ReceiverPosition;
        ObjNew.ReceiverPosition_Type=Obj.ReceiverPosition_Type;
        ObjNew.ReceiverPosition_Units=Obj.ReceiverPosition_Units;
        ObjNew.SourcePosition=Obj.SourcePosition;
        ObjNew.SourcePosition_Type=Obj.SourcePosition_Type;
        ObjNew.SourcePosition_Units=Obj.SourcePosition_Units;
        ObjNew.EmitterPosition=Obj.EmitterPosition;
        ObjNew.EmitterPosition_Type=Obj.EmitterPosition_Type;
        ObjNew.EmitterPosition_Units=Obj.EmitterPosition_Units;
        ObjNew.ListenerUp=Obj.ListenerUp;
        ObjNew.ListenerView=Obj.ListenerView;
        ObjNew.ListenerView_Type=Obj.ListenerView_Type;
        ObjNew.SourceUp=Obj.SourceUp;
        ObjNew.SourceView=Obj.SourceView;
        ObjNew.SourceView_Type=Obj.SourceView_Type;
        if isfield(Obj,'RoomCornerA')
            if size(Obj.RoomCornerA,1) == 3
                ObjNew.RoomCornerA=Obj.RoomCornerA'; % must be transposed for the data existing so far
            else
                ObjNew.RoomCornerA=Obj.RoomCornerA;
            end
            ObjNew.RoomCornerA_Type=Obj.RoomCornerA_Type;
            ObjNew.RoomCornerA_Units=Obj.RoomCornerA_Units;
        end
        
        if isfield(Obj,'RoomCornerB')
            if size(Obj.RoomCornerB,1) == 3
                ObjNew.RoomCornerB=Obj.RoomCornerB'; % must be transposed for the data existing so far
            else
                ObjNew.RoomCornerB=Obj.RoomCornerB;
            end
            ObjNew.RoomCornerB_Type=Obj.RoomCornerB_Type;
            ObjNew.RoomCornerB_Units=Obj.RoomCornerB_Units;
        end
        
        %% Update dimensions
        Obj=SOFAupdateDimensions(ObjNew); % overwrite original Object
        modified=1;
        warning('SOFA:upgrade','Conventions SingleRoomDRIR upgraded to SingleRoomSRIR.  Use warning(''off'',''SOFA:upgrade''); to switch off this warning. ');
    end
end
