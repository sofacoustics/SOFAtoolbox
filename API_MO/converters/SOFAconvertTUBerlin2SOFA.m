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
Obj.GLOBAL_ListenerShortName = 'KEMAR';
Obj.GLOBAL_History='Converted from the TU-Berlin format';
Obj.GLOBAL_Comment = irs.description;
Obj.GLOBAL_License = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0';
Obj.GLOBAL_ApplicationName = 'HRTF from TU Berlin';
Obj.GLOBAL_ApplicationVersion = '1.0';
Obj.GLOBAL_AuthorContact = 'hagen.wierstorf@tu-berlin.de';
Obj.GLOBAL_References = ['H. Wierstorf, M. Geier, A. Raake, S. Spors. ', ...
    'A Free Database of Head-Related Impulse Response Measurements in ', ...
    'the Horizontal Plane with Multiple Distances. ', ...
    'In 130th Convention of the Audio Engineering Society, May 2011.'];
Obj.GLOBAL_Origin = 'https://dev.qu.tu-berlin.de/projects/measurements/repository/show/2010-11-kemar-anechoic/mat';
Obj.GLOBAL_DatabaseName = 'TU Berlin';
Obj.GLOBAL_Title = 'HRTF';
Obj.GLOBAL_ListenerDescription = irs.head;
Obj.GLOBAL_ReceiverDescription = irs.ears;
Obj.GLOBAL_SourceDescription = irs.source;

%% Fill the mandatory variables
% SimpleFreeFieldHRIR 0.6
% number of measurements
M = length(irs.apparent_elevation);
distance = sqrt(sum((irs.source_position-irs.head_position).^2));
Obj.SourcePosition = [nav2sph(rad2deg(irs.apparent_azimuth)') ...
    rad2deg(irs.apparent_elevation)' ...
    distance.*ones(M,1)];
Obj.ListenerPosition = [0 0 0];
Obj.ListenerView = [1 0 0];
Obj.ListenerUp = [0 0 1];

  
%% Update dimensions
Obj=SOFAupdateDimensions(Obj);
