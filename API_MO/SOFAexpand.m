function [Obj, log] = SOFAexpand(Obj,VarName)
%SOFAexpand
%   Obj = SOFAexpand(Obj) expands the singleton dimensions of all variables.
%   Only variables will be expanded. Data and attributes won't. Note that
%   also Obj.API.Dimensions will be updated to the new dimensions.
% 
%		Obj = SOFAexpand(Obj,VarName) expands the singleton dimensions of 
%		the variable VarName.
%   
%   [Obj,log] = SOFAexpand(...) returns a log of expanded variables.

% SOFA API - function SOFAexpand
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Initial check 
OC = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'a');
log={''};

%% Update dimensions
Obj=SOFAupdateDimensions(Obj,'nodata');

%% If VarName given, expand a single variable only
if ~exist('VarName','var'),
  
  %% Expand all variables
    % create field names which should have dimensions
  X=rmfield(Obj,{'Data','API'});
  if isfield(X,'PRIVATE'), X=rmfield(X,'PRIVATE'); end
  Xf=fieldnames(X);

  % Update the dimensions structure w/o data
  for ii=1:length(Xf)
    if ~isempty(strfind(Xf{ii},'_')),	continue; end; % is an attribute --> not expandable
    if ~isfield(OC.API.Dimensions, Xf{ii}), continue; end; % is a used-defined variable --> not expandable
    dim=OC.API.Dimensions.(Xf{ii}); % all possible dimensions
    if ~iscell(dim), continue; end;	% is a variable with a single dimension definition --> not expandable
    if numel(dim)==1, continue; end;	% is a variable with a single dimension definition --> not expandable
    [varNew,dimNew]=expand(Obj,Xf{ii},dim);
    if ~isempty(varNew),
      Obj.(Xf{ii})=varNew;
      Obj.API.Dimensions.(Xf{ii})=dimNew;
      log{end+1}=[Xf{ii} ' expanded to ' dimNew];
    end
  end

  % Expand the dimensions of Data
  Xf=fieldnames(Obj.Data);
  for ii=1:length(Xf)
    if ~isempty(strfind(Xf{ii},'_')),	continue; end; % is an attribute --> not expandable
    if ~isfield(OC.API.Dimensions.Data, Xf{ii}), continue; end; % is a used-defined variable --> not expandable
    dim=OC.API.Dimensions.Data.(Xf{ii}); % all possible dimensions
    if ~iscell(dim), continue; end;	% is a variable with a single dimension definition --> not expandable
    if numel(dim)==1, continue; end;	% is a variable with a single dimension definition --> not expandable
    [varNew,dimNew]=expandData(Obj,Xf{ii},dim);
    if ~isempty(varNew),
      Obj.Data.(Xf{ii})=varNew;
      Obj.API.Dimensions.Data.(Xf{ii})=dimNew;
      log{end+1}=['Data.' Xf{ii} ' expanded to ' dimNew];
    end
  end
  
  
else	% Expand a single variable only
	if isempty(strfind(VarName,'_')),	% is an attribute --> not expandable
    if strncmp(VarName,'Data.',length('Data.'))         
      % variable within the Data. structure
      VarName=VarName(length('Data.')+1:end);
      if isfield(OC.API.Dimensions.Data, VarName), % is a used-defined variable --> not expandable
        dim=OC.API.Dimensions.Data.(VarName); % all possible dimensions
        if iscell(dim), % is a variable with a single dimension definition --> not expandable
            [varNew,dimNew]=expandData(Obj,VarName,dim);
            if ~isempty(varNew),
                Obj.Data.(VarName)=varNew;
                Obj.API.Dimensions.Data.(VarName)=dimNew;
                log{end+1}=['Data.' VarName ' expanded to ' dimNew];
            end
        end
      end
    else
      if isfield(OC.API.Dimensions, VarName), % is a used-defined variable --> not expandable
        dim=OC.API.Dimensions.(VarName); % all possible dimensions
        if iscell(dim), % is a variable with a single dimension definition --> not expandable
          [varNew,dimNew]=expand(Obj,VarName,dim);
          if ~isempty(varNew),
              Obj.(VarName)=varNew;
              Obj.API.Dimensions.(VarName)=dimNew;
              log{end+1}=[VarName ' expanded to ' dimNew];
          end
        end
      end
		end
	end
end	

%% log variable
if length(log)>1, log=log(2:end); else log={}; end;

%% expand a single variable
% Obj: the full SOFA object
% f: name of the variable
% dims: allowed dimensions of that variable (cell array)
% var: expanded variable, or empty if nothing happened
% dN: new dimension, or empty if nothing happened
function [var,dN]=expand(Obj,f,dims)
	d=cell2mat(strfind(dims,'I'));	% all choices for a singleton dimensions
	for jj=1:length(d)	% loop through all expandable dimensions
		len=size(Obj.(f),d(jj)); % size of the considered dimension
		if len>1, continue; end;	% the expandable dimension is already expanded
		dN=dims{cellfun('isempty',strfind(dims,'I'))==1};
		var=bsxfun(@times,Obj.(f),ones(getdim(Obj,dN)));
	end
	if ~exist('var','var'), var=[]; dN=[]; end;
%% Get the sizes of the dimension variables according the dimension variables in str
function vec=getdim(Obj,str)
	vec=arrayfun(@(f)(Obj.API.(f)),upper(str));
  
%% expand a single Data variable
% Obj: the full SOFA object
% f: name of the Data variable
% dims: allowed dimensions of that variable (cell array)
% var: expanded variable, or empty if nothing happened
% dN: new dimension, or empty if nothing happened
function [var,dN]=expandData(Obj,f,dims)
	d=cell2mat(strfind(dims,'I'));	% all choices for a singleton dimensions
	for jj=1:length(d)	% loop through all expandable dimensions
		len=size(Obj.Data.(f),d(jj)); % size of the considered dimension
		if len>1, continue; end;	% the expandable dimension is already expanded
		dN=dims{cellfun('isempty',strfind(dims,'I'))==1};
		var=bsxfun(@times,Obj.Data.(f),ones(getdim(Obj,dN)));
	end
	if ~exist('var','var'), var=[]; dN=[]; end;
