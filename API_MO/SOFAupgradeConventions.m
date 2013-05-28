function Obj = SOFAupgradeConventions(Obj)
%SOFAcompatibility 
%   Obj = SOFAcompatibility(Obj) adapts the Obj such that it is represented
%   in the supported SOFA version.
%   

% SOFA API - function SOFAcompatibility
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Convert from SOFA 0.3 to SOFA 0.4
if strcmp(Obj.GLOBAL_Version,'0.3')
  disp('Upgrading from SOFA 0.3 to SOFA 0.4...');
  % in SOFA 0.3, only SimpleFreeFieldHRIR was supported. Adapt this conventions
    % update metadata
  X=SOFAgetConventions('SimpleFreeFieldHRIR');
  att={'GLOBAL_Version', 'GLOBAL_SOFAConventionsVersion'}; %,...
  for ii=1:length(att)
    Obj.(att{ii})=X.(att{ii});
  end
  Obj.GLOBAL_TimeCreated=Obj.GLOBAL_DatabaseTimeCreated;
  Obj.GLOBAL_TimeModified=Obj.GLOBAL_DatabaseTimeModified;
  if isfield(Obj,'GLOBAL_History'), tmp=Obj.GLOBAL_History; else tmp=''; end
  Obj.GLOBAL_History=[tmp '; Upgraded from SOFA 0.3'];
    % remove dimensional variables and not used variables/attributes
  dims={'I','R','E','N','M','C','Q','SourceView',...
    'SourceUp','GLOBAL_DatabaseTimeCreated','GLOBAL_DatabaseTimeModified'}; 
  f=fieldnames(Obj);
  for ii=1:length(dims)
    for jj=1:length(f)
      if strcmp(f{jj},dims{ii}), 
        Obj=rmfield(Obj,f{jj});   % remove variable or attribute
        if isempty(strfind(f{jj},'_')),
          Obj.Dimensions=rmfield(Obj.Dimensions,f{jj}); % remove dimension
        end
      elseif strcmp(f{jj}(1:min(length(dims{ii})+1,length(f{jj}))),[dims{ii} '_']) 
        Obj=rmfield(Obj,f{jj});  % remove attributes of that variable
      end
    end
  end
end