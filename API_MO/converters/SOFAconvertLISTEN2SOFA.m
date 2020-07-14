function Obj=SOFAconvertLISTEN2SOFA(LISTEN, subjectID)
% Obj=SOFAconvertLISTEN2SOFA(LISTEN, subjectID) converts the HRTFs described in LISTEN 
% (see LISTEN HRTF format) to a SOFA object.
%

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Get an empty conventions structure
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

%% Fill data with data
		% content_m is [M N], data.IR must be [M R N]
Obj.Data.IR = zeros(size(LISTEN.l_eq_hrir_S.content_m,1),2,size(LISTEN.l_eq_hrir_S.content_m,2));
Obj.Data.IR(:,2,:)=LISTEN.r_eq_hrir_S.content_m;
Obj.Data.IR(:,1,:)=LISTEN.l_eq_hrir_S.content_m;
Obj.Data.SamplingRate = 48000; % Note: LISTEN.l_eq_hrir_S.sampling_hz contains 44100 which is wrong!

%% Fill with attributes
Obj.GLOBAL_ListenerShortName = subjectID;
Obj.GLOBAL_History='Converted from the LISTEN format';

%% Fill the mandatory variables
  % SimpleFreeFieldHRIR 0.2
    % Obj.ListenerPosition = [1.95 0 0];
    % Obj.ListenerView = [-1 0 0];
    % Obj.ListenerUp = [0 0 1];
    % Obj.ListenerRotation = [LISTEN.l_eq_hrir_S.azim_v LISTEN.l_eq_hrir_S.elev_v zeros(size(LISTEN.l_eq_hrir_S.elev_v,1),1)];
  % SimpleFreeFieldHRIR 0.3
Obj.ListenerPosition = [0 0 0];
Obj.ListenerView = [1 0 0];
Obj.ListenerUp = [0 0 1];
Obj.SourcePosition = [LISTEN.l_eq_hrir_S.azim_v LISTEN.l_eq_hrir_S.elev_v 1.95*ones(size(LISTEN.l_eq_hrir_S.elev_v,1),1)];


%% Update dimensions
Obj=SOFAupdateDimensions(Obj);
