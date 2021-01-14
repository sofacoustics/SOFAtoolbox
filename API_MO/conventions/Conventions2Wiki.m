%% Converts all CSV files with Conventions to WIKI tables
% Piotr Majdak (2013), EUPL
clear all;
p=mfilename('fullpath');
d=dir(['*.csv']);
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
