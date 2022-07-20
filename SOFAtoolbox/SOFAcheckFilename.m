function newfn=SOFAcheckFilename(fn)
%SOFACHECKFILENAME
%   newFN = SOFAcheckFilename(FN) checks the filename FN and:

% SOFA API - function SOFAcheckFilename
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (20.10.2021)
%

% filename = string?
if ~ischar(fn)
	error('Filename must be a string.');
end
if strcmp(fn(1:7),'http://')
  newfn=fn; % remote file
else
  fn=strrep(fn,'/',filesep);
  fn=strrep(fn,'\',filesep);
  if length(fn)>4 
      if strcmpi(fn(1:4),['db:' filesep]), fn=[SOFAdbPath fn(5:end)]; end
  end
  newfn=fn;
end
