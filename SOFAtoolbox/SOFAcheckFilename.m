function newfn=SOFAcheckFilename(fn)
%SOFAcheckFilename - Check and prepare file name
%   Usage: newfn=SOFAcheckFilename(fn)
%   
%   SOFAcheckFilename checks the filename fn, and replaces
%   file separators by the actual file separators of the current system. 
%   If fn begins with 'db:', then it is replaced by the SOFAdbPath. 

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (20.10.2021)
% #Author: Michael Mihocic: fixes for strings instead of characters, and file names shorter than 7 characters (17.08.2022)

% SOFA Toolbox - function SOFAcheckFilename
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.


% filename = string?
if isstring(fn) 
    fn = char(fn);
end
if ~ischar(fn)
	error('Filename must be a string.');
end
if length(fn)>6 && strcmp(fn(1:7),'http://')
  newfn=fn; % remote file
else
  fn=strrep(fn,'/',filesep);
  fn=strrep(fn,'\',filesep);
  if length(fn)>4 
      if strcmpi(fn(1:4),['db:' filesep]), fn=[SOFAdbPath fn(5:end)]; end
  end
  newfn=fn;
end
