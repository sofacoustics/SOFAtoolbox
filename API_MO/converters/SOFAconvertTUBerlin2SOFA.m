function Obj=SOFAconvertTUBerlin2SOFA(irs)
% OBJ=SOFAconvertTUBerlin2SOFA(irs) converts the HRTFs described in irs
% (see TU-Berlin HRTF format) to a SOFA object.
%

% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 


%% Get an empy conventions structure
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

%% Fill data with data
Obj.Data.IR = shiftdim(shiftdim(irs.left,-1),2); % irs.left is [N M], data.IR must be [M R N]
Obj.Data.IR(:,2,:) = shiftdim(shiftdim(irs.right,-1),2);
Obj.Data.SamplingRate = irs.fs;

%% Fill with attributes
Obj.GLOBAL_SubjectID = irs.head;
Obj.GLOBAL_DatabaseName = 'TU-Berlin'; % maybe setting the name by function parameter
Obj.GLOBAL_ApplicationName = 'TUBerlin2SOFA';
Obj.GLOBAL_ApplicationVersion = '1.0';
Obj.GLOBAL_Organization = 'Technische Universität Berlin';
Obj.GLOBAL_AuthorContact = 'hagen.wierstorf@tu-berlin.de';
Obj.GLOBAL_Comment = irs.description;
Obj.GLOBAL_ListenerDescription = irs.head;
Obj.GLOBAL_ReceiverDescription = irs.ears;
Obj.GLOBAL_SourceDescription = irs.source;

%% Fill the mandatory variables
Obj.SourcePosition = irs.source_position';
Obj.SourceView = irs.source_reference';
Obj.SourceUp = irs.source_position' + [0 0 1];
Obj.ListenerPosition = irs.head_position';
Obj.ListenerView = irs.head_reference';
Obj.ListenerUp = irs.head_position' + [0 0 1];
Obj.ListenerRotation = [nav2sph(radtodeg(-irs.apparent_azimuth)') radtodeg(-irs.apparent_elevation)' ...
    zeros(length(irs.apparent_azimuth),1)];

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

%% Fill with some additional data
Obj.GLOBAL_History='Converted from the TU-Berlin format';
