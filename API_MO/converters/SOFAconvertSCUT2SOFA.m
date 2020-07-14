function Obj=SOFAconvertSCUT2SOFA(root,r)
% OBJ=SOFAconvertSCUT2SOFA(root,pinna) loads the SCUT HRTFs saved in a 
% directory ROOT for the radius R and converts to a SOFA object.
% R must be in meters.

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Get an empy conventions structure
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

%% Define elevations
eles = [-30 -15 0 15 30 45 60 75 90];
elecnt=[ 72 72 72 72 72 72 36 24  1];

%% Create empty matrix
M=sum(elecnt)*length(r);
Obj.SourcePosition=zeros(M,3);
Obj.Data.IR=zeros(M,2,512);

%% Fill with data 
Obj.Data.SamplingRate = 44100;
ii=1;
for jj=1:length(r)
  for ei = 1 : length(eles)
    ele = eles(ei);
    for ai = 0 : elecnt(ei)-1
      azi = 360/elecnt(ei)*ai;
      fn=fullfile(root, ['r' num2str(r(jj)*100)], ['ele' num2str(ele)], ['H' num2str(azi) 'c.pcm']);
      dirfn=dir(fn);
      fid = fopen(fn, 'r');
      H = fread(fid,'float');
      fclose(fid);
      Obj.Data.IR(ii,1,:) = single(H(1:2:1024));  % separate the left-ear HRIR 
      Obj.Data.IR(ii,2,:) = single(H(2:2:1024));  % separate the right-ear HRIR
      Obj.SourcePosition(ii,:) = [mod(360-azi,360) ele r(jj)];
      ii=ii+1;    
    end
  end
end

%% Fill with attributes
Obj.GLOBAL_ListenerShortName = 'KEMAR';
Obj.GLOBAL_History='Converted from the SCUT format';

Obj.GLOBAL_Author = 'Bosun Xie';
Obj.GLOBAL_AuthorContact = 'phbsxie@scut.edu.cn';
Obj.GLOBAL_License = 'CC 3.0 BY-SA-NC';
Obj.GLOBAL_Organization = 'South China University of Technology, Guangzhou, China';

Obj.GLOBAL_References = 'Bosun Xie, 2013, "Head-Related Transfer Function and Virtual Auditory Display", J Ross Publishing Inc., Plantation, FL, USA';
Obj.GLOBAL_RoomType = 'free field';
Obj.GLOBAL_Title = 'HRTF';
Obj.GLOBAL_DatabaseName='SCUT';

Obj.GLOBAL_DateCreated = datestr(datenum(dirfn.date),SOFAdefinitions('dateFormat'));


%% Fill the mandatory variables
Obj.ListenerPosition = [0 0 0];
Obj.ListenerView = [1 0 0];
Obj.ListenerUp = [0 0 1];

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);
