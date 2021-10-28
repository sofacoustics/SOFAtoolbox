function dbURL=SOFAdbURL(dbURL)
% dbURL=SOFAdbURL 
%
%   dbURL=SOFAdbURL returns the internet URL to the directory containing
%   HRTFs. The default URL is http://www.sofacoustics.org/data.
% 
%   dbURL=SOFAdbURL(newURL) sets the internet URL to the newURL for further calls
%   of SOFAdbURL.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
%
% SOFA API - function SOFAdbURL
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

persistent CachedURL;

if exist('dbURL','var')
    if strcmp(dbURL,'reset')
        CachedURL='http://www.sofacoustics.org/data';
    else
        CachedURL=dbURL;
    end
elseif isempty(CachedURL)
  CachedURL='http://www.sofacoustics.org/data';
end
dbURL=CachedURL;
