%% Converts all CSV files with Conventions to WIKI tables

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated, license text added (28.10.2021)
%
% SOFA API - function SOFAaddVariable
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


clear all;
p=mfilename('fullpath');
d=dir('*.csv');
conventions={};
for ii=1:length(d)
  dn=d(ii).name;
  conventions{ii}=dn(1:end-4);
end
  
for jj=1:length(conventions)
  fid=fopen([conventions{jj} '.csv']);
  C=textscan(fid,'%s%s%s%s%s%s','Delimiter','\t','Headerlines',1); 
  fclose(fid);
  fid=fopen([conventions{jj} '.txt'],'w');
  fprintf(fid,'{| border="1"\n!Name\n!Default\n![[SOFA_conventions#AnchorFlags|Flags]]\n![[SOFA_conventions#AnchorDimensions|Dimensions]]\n!Type\n!Comment\n');
%   C2=regexprep(C{2},'''', '&rsquo;'); % replace single quota (') by &prime;
  for ii=1:length(C{1})
    fprintf(fid,['|-\n|' C{1}{ii} '||<nowiki>' C{2}{ii} '</nowiki>||' C{3}{ii} '||' C{4}{ii} '||' C{5}{ii} '||' C{6}{ii} '\n']);
  end
  fprintf(fid,'|}');
  fclose(fid);
end
disp ('   ** done **');