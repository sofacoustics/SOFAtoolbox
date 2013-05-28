function Obj = SOFAgetConventions(sofaconventions,flags)
%SOFAgetConventions
%
%	  List = SOFAgetConventions() returns a list with supported conventions.
% 
%   Obj = SOFAgetConventions(sofaconventions,flags) returns a SOFA object
%   with all metadata and data for the corresponding sofaconventions.
% 
%   Obj = SOFAgetConventions(sofaconventions,flags) returns only selected
%   metadata for the conventions with the following encoding:
%				m: mandatory, r: readonly, a: all (default)

% SOFA API - function SOFAgetConventions
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% If not conventions provided, return the list with supported conventions
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

%% If no flags provided, return the conventions with all metadata
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
if ~found,
  p=mfilename('fullpath');
  AllObj.(flags)=load([p(1:length(p)-length(mfilename)) 'conventions' filesep sofaconventions '-' flags '.mat']);
end
Obj=AllObj.(flags).Obj;

%% Get the supported dimension tokens
Def = SOFAdefinitions;

%% Overwrite some special fields
if isfield(Obj,'GLOBAL_TimeCreated'), Obj.GLOBAL_TimeCreated=datestr(now,Def.dateFormat); end
if isfield(Obj,'GLOBAL_APIVersion'), Obj.GLOBAL_APIVersion=SOFAgetVersion; end
if isfield(Obj,'GLOBAL_APIName'), Obj.GLOBAL_APIName=Def.APIName; end
