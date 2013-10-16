function SOFAcompileConventions(conventions)
%SOFAcompileConventions
%
%   Obj = SOFAcompileConventions(conventions) compiles SOFA conventions
%   given by a CSV file to .mat files used later by SOFAgetConventions.
% 
%   The CSV file must be in the directory conventions and have the same
%   filename as conventions. SOFAcompileConventions generates 3 files, one
%   for each flag (r, m, and all). 
%
%   Before compiling, SOFAcompileConventions checks if the modification
%   date of the .mat files is older than that of the .csv file. Compiling
%   is not performed if all .mat files are newer than the .csv file. This
%   behaviour is required for operation in a read-only directory. 

% SOFA API 
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

p=mfilename('fullpath');

if ~exist('conventions','var')
  p=mfilename('fullpath');
  d=dir([p(1:length(p)-length(mfilename)) 'conventions' filesep '*.csv']);
  conventions={};
  for ii=1:length(d)
    dn=d(ii).name;
    [~,name,ext]=fileparts(dn);
    flags='ram';
    okcnt=0;
    for jj=1:length(flags)
      x=dir([p(1:length(p)-length(mfilename)) 'conventions' filesep name '-' flags(jj) '.mat']);
      if ~isempty(x)
        if x(1).datenum>d(ii).datenum, okcnt=okcnt+1; end
      end
    end
    if okcnt~=3, conventions{end+1}=name; end
  end
else
  if ~iscell(conventions), conventions={conventions}; end;
end
  
for jj=1:length(conventions)
  fid=fopen([p(1:length(p)-length(mfilename)) 'conventions' filesep conventions{jj} '.csv']);
  C=textscan(fid,'%s%s%s%s%s%s','Delimiter','\t','Headerlines',1);
  fclose(fid);
  Obj=compileConvention(C,'r');
  if strcmp(Obj.GLOBAL_SOFAConventions,conventions{jj})
    disp(['Compiling ' conventions{jj} ' ' Obj.GLOBAL_SOFAConventionsVersion ' for SOFA ' Obj.GLOBAL_Version]);
    save([p(1:length(p)-length(mfilename)) 'conventions' filesep conventions{jj} '-r.mat'],'Obj','-v7');    
    flags='am';
    for ii=1:2
      Obj=compileConvention(C,flags(ii));
      save([p(1:length(p)-length(mfilename)) 'conventions' filesep conventions{jj} '-' flags(ii) '.mat'],'Obj','-v7');
    end
  end
end

%%
function Obj=compileConvention(C,flags)
     
flagc=cellstr((flags)');

% create object structure
for ii=1:size(C{1,1},1)
  C{3}{ii}=[C{3}{ii} 'a']; 
	if  ~isempty(cell2mat(regexp(C{3}{ii},flagc)))
    var=regexprep(C{1}{ii},':','_');
    switch lower(C{5}{ii})
      case 'double'
      C{2}{ii}=str2num(C{2}{ii}); % convert default to double
      case 'string'
      C{2}{ii}={C{2}{ii}};
    end
    if isempty(strfind(var,'Data.'))
      Obj.(var)=C{2}{ii};
      if isempty(strfind(var,'_')) % && ~sum(strcmp(var,dims))
        x2=regexprep(C{4}{ii},' ',''); %  remove spaces
        y=regexprep(x2,',',['''' 10 '''']); % enclose in quotations and insert line breaks
        Obj.API.Dimensions.(var)=eval(['{''' y '''}']);
      end
    else      
      Obj.Data.(var(6:end))=C{2}{ii};
      if isempty(strfind(var(6:end),'_')) 
        x2=regexprep(C{4}{ii},' ',''); %  remove spaces
        y=regexprep(x2,',',['''' 10 '''']); % enclose in quatations and insert line breaks
        Obj.API.Dimensions.Data.(var(6:end))=eval(['{''' y '''}']);
      end      
    end
	end
end


% Overwrite some special fields
if isfield(Obj,'GLOBAL_APIVersion'), Obj.GLOBAL_APIVersion=SOFAgetVersion; end
if isfield(Obj,'GLOBAL_APIName'), Obj.GLOBAL_APIName='ARI Matlab/Octave API'; end

% Create dimension size variables - if not read-only
if flags=='r', return; end
  % fix dimension sizes
Obj.API.I=1;
Obj.API.C=3;
  % variable-dependent dimension sizes
dims='renm'; 
  % check all metadata variables
f=fieldnames(rmfield(Obj.API.Dimensions,'Data'));
for ii=1:length(dims)
	for jj=1:length(f)
		dim=strfind(Obj.API.Dimensions.(f{jj}),dims(ii));
		if iscell(dim), dim=cell2mat(dim); end;
		if ~isempty(dim)
			Obj.API.(upper(dims(ii)))=size(Obj.(f{jj}),dim(1));
			break;
		end
	end
end
  % check all data variables
fd=fieldnames(Obj.API.Dimensions.Data);
for ii=1:length(dims)
	for jj=1:length(fd)
		dim=strfind(Obj.API.Dimensions.Data.(fd{jj}),dims(ii));
		if iscell(dim), dim=cell2mat(dim); end;
		if ~isempty(dim)
			Obj.API.(upper(dims(ii)))=size(Obj.Data.(fd{jj}),dim(1));
			break;
		end
	end
end
