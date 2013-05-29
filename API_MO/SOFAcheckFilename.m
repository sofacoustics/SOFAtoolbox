function newfn=SOFAcheckFilename(fn)
%SOFACHECKFILENAME
%   newFN = SOFAcheckFilename(FN) checks the filename FN and:

% SOFA API - function SOFAcheckFilename
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

% filename = string?
if ~ischar(fn)
	error('Filename must be a string.');
end

newfn=fn;

