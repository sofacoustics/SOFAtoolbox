function filename=SOFAcheckFilename(filename)
%SOFACHECKFILENAME
%   filename = SOFAcheckFilename(filename) checks the specified filename.

% SOFA API - function SOFAcheckFilename
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence.

if ~ischar(filename)
	error('Filename must be a string.');
end

idx=strfind(filename,'.sofa');
if isempty(idx)
    filename=[filename '.sofa'];
% elseif ~strcmp(filename(idx+1:end),'sofa')
%     error(['SOFA-API does not support *' filename(idx:end) '-files!'])
end

end %of function