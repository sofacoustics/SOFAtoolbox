function Obj=SOFAconvertLISTEN2SOFA(LISTEN, subjectID)
% Obj=SOFAconvertLISTEN2SOFA(LISTEN, subjectID) converts the HRTFs described in LISTEN 
% (see LISTEN HRTF format) to a SOFA object.
%

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 


%% Get an empy conventions structure
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

%% Fill data with data
		% content_m is [M N], data.IR must be [M R N]
Obj.Data.IR = zeros(size(LISTEN.l_eq_hrir_S.content_m,1),2,size(LISTEN.l_eq_hrir_S.content_m,2));
Obj.Data.IR(:,2,:)=LISTEN.r_eq_hrir_S.content_m;
Obj.Data.IR(:,1,:)=LISTEN.l_eq_hrir_S.content_m;
Obj.Data.SamplingRate = LISTEN.l_eq_hrir_S.sampling_hz;

%% Fill with attributes
Obj.GLOBAL_SubjectID = subjectID;
Obj.GLOBAL_DatabaseName = 'LISTEN';
Obj.GLOBAL_ApplicationName = 'LISTEN2SOFA';
% Obj.GLOBAL_ApplicationVersion = '';
Obj.GLOBAL_Organization = 'Acoustics Research Institute';
Obj.GLOBAL_AuthorContact = 'piotr@majdak.com';

%% Fill the mandatory variables
Obj.ListenerPosition = [1.95 0 0];
Obj.ListenerView = [0 0 0];
Obj.ListenerUp = [1.95 0 1];
Obj.ListenerRotation = [LISTEN.l_eq_hrir_S.azim_v LISTEN.l_eq_hrir_S.elev_v zeros(size(LISTEN.l_eq_hrir_S.elev_v,1),1)];

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

%% Fill with some additional data
Obj.GLOBAL_History='Converted from the LSITEN format';


% LISTEN2AMTatARI
% script for converting LISTEN database to AMTatARI file system
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Harald Ziegelwanger, OEAW Acoustical Research Institute
% latest update: 2012-02-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     meta.pos=zeros(size(l_eq_hrir_S.elev_v,1));
%     meta.pos(:,1)=l_eq_hrir_S.azim_v;
%     meta.pos(:,2)=l_eq_hrir_S.elev_v;
%     meta.pos(:,3)=NaN(size(l_eq_hrir_S.elev_v,1),1);
%     [meta.pos(:,6),meta.pos(:,7)]=geo2hor(meta.pos(:,1),meta.pos(:,2));
%     
%     hM=zeros(size(l_eq_hrir_S.content_m,2),size(l_eq_hrir_S.content_m,1),2);
%     hM(:,:,1)=transpose(l_eq_hrir_S.content_m);
%     hM(:,:,2)=transpose(r_eq_hrir_S.content_m);
%     
% end


