function [Obj, log] = SOFAcompact(Obj)
%SOFAcompact
%   Obj = SOFAcompact(Obj) compacts the unique value to singleton dimensions 
%   of variables where possible. 
%   Current limitation: Variables with 3 dimensions will be only compacted
%   when the third dimension is the compressible one.
%

% SOFA API - function SOFAcompact
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% Initial check 
OC = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'a');
log={''};

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);
dims=SOFAdefinitions('dimensions');

%% compact w/o data 
X=rmfield(Obj,{'Data','API'});
Xf=fieldnames(X);
for ii=1:length(Xf)
	if ~isempty(strfind(Xf{ii},'_')),	continue; end; % is an attribute --> not compressible
	if ~isfield(OC.API.Dimensions, Xf{ii}), continue; end; % is a used-defined variable --> not compressible
	dims=OC.API.Dimensions.(Xf{ii}); % get all possible dimensions
	if ~iscell(dims), continue; end;   %  variable with a single dimension definition --> not compressible
  if numel(dims)==1, continue; end;	% variable with a single dimension definition --> not compressible
  switch length(dims{1})  % how many dimensions do we have?
    case 1
      n=unique(Obj.(Xf{ii}));
      if length(n)>1, continue; end;  % entries not unique --> not compressible
      Obj.(Xf{ii})=n; % compressed!
    case 2
      d=cell2mat(strfind(dims,'I'));	% get all choices for a singleton dimensions
      if strcmp(dims{d}(1),'I'), n=unique(Obj.(Xf{ii}),'rows'); else n=unique(Obj.(Xf{ii})','rows')'; end;
      if size(n,1)>1, continue; end;  % entries not unique --> not compressible
      if strcmp(dims{d}(1),'I'), Obj.(Xf{ii})=n; else Obj.(Xf{ii})=n'; end; % compressed!
    case 3
      d=cell2mat(strfind(dims,'I'));
      switch d
        case 3
          y=unique(reshape(Obj.(Xf{ii}),[],size(Obj.(Xf{ii}),d))','rows');
          if size(y,1)>1, continue; end; % entries not unique --> not compressible
          Obj.(Xf{ii})=reshape(y',size(Obj.(Xf{ii}),1),size(Obj.(Xf{ii}),2));
        otherwise
          warning('SOFA:compact',[Xf{ii} ' not compressed (currently limited)']);
      end
  end
end

%% Compact the data
Xf=fieldnames(Obj.Data);
for ii=1:length(Xf)
	if ~isempty(strfind(Xf{ii},'_')),	continue; end; % is an attribute --> not compressible
	if ~isfield(OC.API.Dimensions.Data, Xf{ii}), continue; end; % is a used-defined variable --> not compressible
	dims=OC.API.Dimensions.Data.(Xf{ii}); % get all possible dimensions
	if ~iscell(dims), continue; end;   %  variable with a single dimension definition --> not compressible
  if numel(dims)==1, continue; end;	% variable with a single dimension definition --> not compressible
  switch length(dims{1})  % how many dimensions do we have?
    case 1
      n=unique(Obj.Data.(Xf{ii}));
      if length(n)>1, continue; end;  % entries not unique --> not compressible
      Obj.Data.(Xf{ii})=n; % compressed!
    case 2
      d=cell2mat(strfind(dims,'I'));	% all choices for a singleton dimensions
      if strcmp(dims{d}(1),'I'), n=unique(Obj.Data.(Xf{ii}),'rows'); else n=unique(Obj.Data.(Xf{ii})','rows')'; end;
      if size(n,1)>1, continue; end;  % entries not unique --> not compressible
      if strcmp(dims{d}(1),'I'), Obj.Data.(Xf{ii})=n; else Obj.Data.(Xf{ii})=n'; end; % compressed!
    case 3
      % missing
      warning('SOFA:compact',['Data.' Xf{ii} ' not compressed (functionality limited)']);
  end
end

%% clean up
Obj=SOFAupdateDimensions(Obj);
if length(log)>1, log=log(2:end); else log={}; end;

function vec=getdim(Obj,str)
vec=NaN(1,length(str));
for ii=1:length(str)
	vec(ii)=Obj.(upper(str(ii)));
end