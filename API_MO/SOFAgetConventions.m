function Obj = SOFAgetConventions(sofaconventions,flags)
%SOFAgetConventions
%
%    List = SOFAgetConventions() returns a list with supported conventions.
% 
%    Obj = SOFAgetConventions(sofaconvention) returns a SOFA object
%    with all metadata and data for the corresponding sofaconvention. Obj
%    will be empty if sofaconventions is not supported.
% 
%    Obj = SOFAgetConventions(sofaconvention,flags) returns only selected
%    metadata for the corresponding sofaconvention with the following encoding:
%        m: mandatory
%        r: readonly
%        a: all (default)

% SOFA API - function SOFAgetConventions
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% If sofaconventions not provided, return the list with supported conventions
if ~exist('sofaconventions','var')
  p=mfilename('fullpath');
  d=dir([p(1:length(p)-length(mfilename)) 'conventions' filesep '*-m.mat']);
  Obj={};
  for ii=1:length(d)
    dn=d(ii).name;
    Obj{ii}=dn(1:end-6);
  end
	return;
end

%% If flags not provided, return the conventions with all metadata
if ~exist('flags','var')
    flags='a'; % flags: m: mandatory, r: readonly, a: all
end

%% Load cached object
persistent AllObj;

found=0;
if isfield(AllObj,flags)
  if isfield(AllObj.(flags).Obj,'GLOBAL_SOFAConventions')
    if strcmp(AllObj.(flags).Obj.GLOBAL_SOFAConventions,sofaconventions)
      found=1;
    end
  end
end
if found,
  Obj=AllObj.(flags).Obj; % return cached convention object
else
  p=mfilename('fullpath');
  if ~isempty(dir([p(1:length(p)-length(mfilename)) 'conventions' filesep sofaconventions '-' flags '.mat']))
    AllObj.(flags)=load([p(1:length(p)-length(mfilename)) 'conventions' filesep sofaconventions '-' flags '.mat']);
    Obj=AllObj.(flags).Obj; % load conventions to the cache and return
  else
    Obj=[]; % return empty array when conventions not found
  end
end


%% Overwrite some special fields
if isfield(Obj,'GLOBAL_DateCreated'), Obj.GLOBAL_DateCreated=datestr(now,SOFAdefinitions('dateFormat')); end
if isfield(Obj,'GLOBAL_APIVersion'), Obj.GLOBAL_APIVersion=SOFAgetVersion; end
if isfield(Obj,'GLOBAL_APIName'), Obj.GLOBAL_APIName=SOFAdefinitions('APIName'); end
