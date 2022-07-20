function Obj = SOFAgetConventions(sofaconventions,flags,version)
%SOFAgetConventions
%
%    List = SOFAgetConventions() returns a list with supported conventions.
% 
%    Obj = SOFAgetConventions(sofaconvention) returns a SOFA object
%    with all metadata and data for the corresponding sofaconvention. Obj
%    will be empty if sofaconventions is not supported.
% 
%    Obj = SOFAgetConventions(sofaconvention, flags) returns only selected
%    metadata for the corresponding sofaconvention with the following encoding:
%        m: mandatory
%        r: readonly
%        a: all (default)
%
%    Obj = SOFAgetConventions(sofaconvention, version) returns only the selected
%    version of the convention. 
%
%    Obj = SOFAgetConventions(sofaconvention, flag, version) returns only the selected
%    version and metadata. 

% #Author: Piotr Majdak
% #Author: Michael Mihocic: doc & header documentation updated (28.10.2021)
%
% SOFA API - function SOFAgetConventions
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Cache
persistent AllObj;

%% If flags not provided, return the conventions with all metadata
if ~exist('version','var')
  if ~exist('flags','var')
    flags='a'; % no flags, no version --> flags = a, version = any.
  else
    switch flags % no version, check if flags is flags or version
      case {'a', 'm', 'r'}
        % flags is flags --> do nothing
      otherwise
        version=flags;  % flags is version --> copy to version 
        flags='a';
    end
  end
end


%% If sofaconventions not provided, return the list with supported conventions
if ~exist('sofaconventions','var')
  p=mfilename('fullpath');
  d=dir([p(1:length(p)-length(mfilename)) 'conventions' filesep '*_m_*.mat']);
  Obj={};
  for ii=1:length(d)
    dn=d(ii).name; 
    idx=strfind(dn,'_');
    Obj{ii,1}=dn(1:idx(1)-1);
  end
  Obj=unique(Obj);
  AllObj=[];
	return;
end

%% Load cached object

if isfield(AllObj,flags)
  if isfield(AllObj.(flags).Obj,'GLOBAL_SOFAConventions')
    if strcmp(AllObj.(flags).Obj.GLOBAL_SOFAConventions,sofaconventions)
      if ~exist('version','var')
        Obj=AllObj.(flags).Obj; % return cached convention object
      else
        if strcmp(AllObj.(flags).Obj.GLOBAL_SOFAConventionsVersion,version)
          Obj=AllObj.(flags).Obj; % return cached convention object
        end
      end
    end
  end
end

if ~exist('Obj','var')
        % cached object not loaded yet
  p=mfilename('fullpath');
  if exist('version','var')
      % load a specific version but do not store in the cache
    if ~isempty(dir([p(1:length(p)-length(mfilename)) 'conventions' filesep sofaconventions '_' flags '_' version '.mat']))
      load([p(1:length(p)-length(mfilename)) 'conventions' filesep sofaconventions '_' flags '_' version '.mat']);      
    else
      error(['Convention ' sofaconventions ' with the version ' version ' not found.']);
    end
  else
    allver=dir(fullfile(p(1:length(p)-length(mfilename)), 'conventions' , [sofaconventions '_' flags '_*.mat']));
    if isempty(allver)
      warning(['Convention ' sofaconventions ' not found.']');
      Obj=[];
    else
      vermax='0.0';
      idxmax=0;
      for ii=1:length(allver)
        x=load([p(1:length(p)-length(mfilename)) 'conventions' filesep allver(ii).name]);
        if compareversions(x.Obj.GLOBAL_SOFAConventionsVersion,vermax)>0
          vermax=x.Obj.GLOBAL_SOFAConventionsVersion;
          idxmax=ii;
        end
      end
        % store in cache
      AllObj.(flags)=load([p(1:length(p)-length(mfilename)) 'conventions' filesep allver(idxmax).name]);
      Obj=AllObj.(flags).Obj; % return the cached version
    end
  end
end

%% Overwrite some special fields
if isfield(Obj,'GLOBAL_DateCreated'), Obj.GLOBAL_DateCreated=datestr(now,SOFAdefinitions('dateFormat')); end
if isfield(Obj,'GLOBAL_APIVersion'), Obj.GLOBAL_APIVersion=SOFAgetVersion; end
if isfield(Obj,'GLOBAL_APIName'), Obj.GLOBAL_APIName=SOFAdefinitions('APIName'); end

end

function res=compareversions(ver1, ver2)

  [~,maj1,min1]=fileparts(ver1);
  [~,maj2,min2]=fileparts(ver2);
  maj1=str2num(maj1);
  min1=str2num(min1(2:end));
  maj2=str2num(maj2);
  min2=str2num(min2(2:end));
  if maj1>maj2 
     res=1;   % ver1 > ver2     
  elseif maj1==maj2 && min1>min2
     res=1;   % ver1 > ver2
  elseif maj2==maj2 && min1==min2
     res=0;   % ver1 == ver2
  else 
     res=-1;  % ver < ver2
  end
end