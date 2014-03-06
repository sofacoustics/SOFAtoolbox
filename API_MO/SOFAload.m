function Obj = SOFAload(fn,flags)
%SOFALOAD 
%   Obj = SOFAload(FN) reads the SOFA object OBJ with all data from
%   a SOFA file FN. 
%
%   FN can point to a remote file (containing '://') or to a local file: 
%     Remote file: FN will be downloaded to a temporary directory and
%       loaded.
%     Local file: If existing, will be loaded. If not existing, it will be
%       downloaded from the internet repository given by SOFAdbURL. For
%       this, FN must begin with the local HRTF directory given by SOFAdbPath.
%
%   Obj = SOFAload(FN,'nodata') ignores the "Data." variables and
%   loads metadata only (variables and attributes).
%
%   Obj = SOFAload(FN,'nochecks') loads the file but does not perform
%   any checks for correct conventions.
% 
%   Obj = SOFAload(FN,[START COUNT]) loads only COUNT number of 
%		measurements beginning with the index START. For remote files, or local
%   but not existing files, the full file will be downloaded.
%   

% SOFA API - function SOFAload
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Global definitions
% dims={'i';'r';'e';'n';'m';'c';'q'}; % dimensions
if ~exist('flags','var'), flags='all'; end;

%% check file name
fn=SOFAcheckFilename(fn);

%% check if local or remote and download if necessary
if strfind(fn,'://')
  % remote path: download as temporary
  newfn=[tempname '.sofa'];
  urlwrite(fn, newfn);
else
  newfn=fn;
  if ~exist(fn,'file') % file does not exist? 
    warning(['File not found: ' fn]);
      % local path: replace SOFAdbPath by SOFAdbURL, download to SOFAdbPath 
    if length(fn)>length(SOFAdbPath) % fn is longer than SOFAdbPath?
      if strcmp(SOFAdbPath,fn(1:length(SOFAdbPath))) % fn begins with SOFAdbPath
        webfn=fn(length(SOFAdbPath)+1:end);
        webfn(strfind(webfn,'\'))='/';
        webfn=[SOFAdbURL regexprep(webfn,' ','%20')];        
        disp(['Downloading ' fn(length(SOFAdbPath)+1:end) ' from ' SOFAdbURL]);
        [f,stat]=urlwrite(webfn,fn);
        if ~stat
          error(['Could not download file: ' webfn]);
        end
      else % fn not existing and not beginning with SOFAdbPath --> error
        error(['Unable to read file ''' fn ''': no such file']);
      end
    else % fn not existing and shorter than SOFAdbPath --> error
        error(['Unable to read file ''' fn ''': no such file']);
    end
  end
end 

%% Load the object
[Obj]=NETCDFload(newfn,flags);

%% Return if no checks should be performed
if strcmp(flags, 'nochecks')
  return;
end

%% Check if SOFA conventions
if ~isfield(Obj,'GLOBAL_Conventions'), error('File is not a valid SOFA file'); end
if ~strcmp(Obj.GLOBAL_Conventions,'SOFA'), error('File is not a valid SOFA file'); end
if ~isfield(Obj,'GLOBAL_SOFAConventions'), error('Information about SOFA conventions is missing'); end
try
	X=SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'m');
catch ME
	error('Unsupported SOFA conventions');
end

%% Check if DataType is present and correct
if ~isfield(Obj,'GLOBAL_DataType'), error('DataType is missing'); end
if ~strcmp(Obj.GLOBAL_DataType,X.GLOBAL_DataType), error('Wrong DataType'); end

%% Ensure backwards compatibility
modified=1;
while modified, [Obj,modified]=SOFAupgradeConventions(Obj); end

%% If data loaded, check if correct data format
if ~strcmp(flags,'nodata'),
	if ~isfield(Obj,'Data'), error('Data is missing'); end
	f=fieldnames(X.Data);
	for ii=1:length(f)
		if ~isfield(Obj.Data,f{ii})
      Obj.Data.(f{ii})=X.Data.(f{ii});
      Obj.API.Dimensions.Data.(f{ii})=X.API.Dimensions.Data.(f{ii});
      warning(['Data.' f{ii} ' was missing, set to default']);
    end
	end
	Obj=SOFAupdateDimensions(Obj);
end

%% Check if mandatory variables are present
f=fieldnames(X);
for ii=1:length(f)
	if ~isfield(Obj,f{ii})
    Obj.(f{ii})=X.(f{ii});
    if ~strfind(f{ii},'_'),
      Obj.API.Dimensions.(f{ii})=X.API.Dimensions.(f{ii});
    end
    warning([f{ii} ' was missing, set to default']);
  end
end

