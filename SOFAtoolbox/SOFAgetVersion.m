function results = SOFAgetVersion(flag)
%SOFAgetVersion - Get the version of the SOFA Toolbox
%   Usage: ver = SOFAgetVersion(flag)
%
%   ver = SOFAgetVersion() returns the version of the SOFA Toolbox
%
%   ver = SOFAgetVersion('Toolbox') and SOFAgetVersion('API') do the same (not case-sensitive).
%
%   ver = SOFAgetVersion('SOFA') returns the version of the SOFA supported 
%   by the SOFA Toolbox.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
%
% SOFA Toolbox - function SOFAgetVersion
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

if ~exist('flag','var')
	flag='TOOLBOX';
end

switch upper(flag)
	case {'TOOLBOX','API'}
		results = '2.2.0';
	case 'SOFA'
		results = '2.1';
end
