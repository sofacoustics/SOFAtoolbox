function [Obj,modified] = SOFAupgradeConventions(Obj)
%SOFAcompatibility 
%   [Obj,modified] = SOFAupgradeConventions(Obj) upgrades the Obj to the next higher
%   version if required. MODIFIED is 1 when an upgrade was required. 
%   In order to obtain the most recent version, SOFAupgradeConventions
%   should be processed recursively until MODIFIED is 0. 


% SOFA API - function SOFAcompatibility
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
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
    warning('SOFA:upgrade','SOFA 0.3 upgraded to 0.4');
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
    warning('SOFA:upgrade','SOFA 0.4 upgraded to 0.5');
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
    warning('SOFA:upgrade','SOFA 0.5 upgraded to 0.6');
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
      warning('SOFA:upgrade','SOFA 0.6 upgraded to 1.0');    
    else
      warning('SOFA:upgrade','Unknown conventions');
    end
end

%% Upgrade specific to conventions
if ~modified
  switch Obj.GLOBAL_SOFAConventions
    case 'MultiSpeakerBRIR'
      if strcmp(Obj.GLOBAL_SOFAConventionsVersion,'0.1');
          % upgrade to 0.2
        Obj.GLOBAL_DataType='FIRE';
        Obj.GLOBAL_SOFAConventionsVersion='0.2';
        %Obj.Data.Delay = 
        if strcmp(Obj.API.Dimensions.Data.Delay,'IR')
          Obj.API.Dimensions.Data.Delay='IRE'; 
          Obj.Data.Delay=repmat(Obj.Data.Delay,[1 1 size(Obj.EmitterPosition,1)]);
        end
        if strcmp(Obj.API.Dimensions.Data.Delay,'MR')
          Obj.API.Dimensions.Data.Delay='MRE'; 
          Obj.Data.Delay=repmat(Obj.Data.Delay,[1 1 size(Obj.EmitterPosition,1)]);
        end
        modified=1;
        warning('SOFA:upgrade','Conventions MultiSpeakerBRIR 0.1 upgraded to 0.2');
      end
  end
end
