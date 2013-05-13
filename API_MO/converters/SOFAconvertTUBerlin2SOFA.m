function Obj=SOFAconvertTUBerlin2SOFA(irs)
% OBJ=SOFAconvertTUBerlin2SOFA(irs) converts the HRTFs described in irs
% (see TU-Berlin HRTF format) to a SOFA object.
%

% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Get an empy conventions structure
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

%% Fill data with data
Obj.Data.IR = shiftdim(shiftdim(irs.left,-1),2); % irs.left is [N M], data.IR must be [M R N]
Obj.Data.IR(:,2,:) = shiftdim(shiftdim(irs.right,-1),2);
Obj.Data.SamplingRate = irs.fs;

%% Fill with attributes
Obj.GLOBAL_SubjectID = irs.head;
Obj.GLOBAL_History='Converted from the TU-Berlin format';
Obj.GLOBAL_Comment = irs.description;
Obj.GLOBAL_ListenerDescription = irs.head;
Obj.GLOBAL_ReceiverDescription = irs.ears;
Obj.GLOBAL_SourceDescription = irs.source;

%% Fill the mandatory variables
Obj.SourcePosition = irs.source_position';
Obj.SourceView = irs.source_reference';
Obj.SourceUp = [0 0 1];
Obj.ListenerPosition = irs.head_position';
Obj.ListenerView = irs.head_reference';
Obj.ListenerUp = [0 0 1];
Obj.ListenerRotation = [nav2sph(rad2deg(-irs.apparent_azimuth)') ...
    rad2deg(-irs.apparent_elevation)' ...
    zeros(length(irs.apparent_azimuth),1)];

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);
