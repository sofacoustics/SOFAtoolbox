function results = SOFAgetVersion(flag)
%SOFAGETVERSION
%   version = SOFAgetVersion() returns the version of the SOFA API
%		version = SOFAgetVersion('API') does the same.
%
%   version = SOFAgetVersion('SOFA') returns the version of the SOFA
%   supported by this API.

% SOFA API - function SOFAgetVersion
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

if ~exist('flag','var')
	flag='API';
end

switch flag
	case 'API'
		results = '0.2.0';
	case 'SOFA'
		results = '0.3';
end
