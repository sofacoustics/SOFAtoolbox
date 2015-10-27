function Obj=SOFAconvertTUBerlinBRIR2SOFA(irs)
% OBJ=SOFAconvertTUBerlin2SOFA(irs) converts the HRTFs described in irs
% (see TU-Berlin HRTF format) to a SOFA object, using the MultiSpeakerBRIR
% Convention.
%

% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Get an empy conventions structure
Obj = SOFAgetConventions('MultiSpeakerBRIR');

%% Fill data with data
Obj.Data.IR = zeros(size(irs.left,2),2,1,size(irs.left,1));
Obj.Data.IR(:,1,:) = shiftdim(shiftdim(irs.left,-2),3); % irs.left is [N M], data.IR must be [M R E N]
Obj.Data.IR(:,2,:) = shiftdim(shiftdim(irs.right,-2),3);
Obj.Data.SamplingRate = irs.fs;

%% Fill with attributes
Obj.GLOBAL_ListenerShortName = 'KEMAR';
Obj.GLOBAL_History='Converted from the TU-Berlin format';
Obj.GLOBAL_Comment = irs.description;
Obj.GLOBAL_License = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0';
Obj.GLOBAL_ApplicationName = 'BRIR from TU Berlin';
Obj.GLOBAL_ApplicationVersion = '1.0';
Obj.GLOBAL_AuthorContact = 'hagen.wierstorf@tu-berlin.de';
Obj.GLOBAL_References = [''];
Obj.GLOBAL_Origin = 'TU Berlin';
Obj.GLOBAL_Organization = 'Quality and Usability Lab, Technische Universitaet Berlin';
Obj.GLOBAL_DatabaseName = 'TU Berlin';
Obj.GLOBAL_Title = 'BRIR TU Berlin';
Obj.GLOBAL_ListenerDescription = irs.head;
Obj.GLOBAL_ReceiverDescription = 'Large ears (KB0065 + KB0066) with G.R.A.S. 40AO pressure microphones';
Obj.GLOBAL_SourceDescription = 'Genelec 8030A';
Obj.GLOBAL_RoomType = 'reverberant';
Obj.GLOBAL_RoomDescription = '';


%% Fill the mandatory variables
% MultiSpeakerBRIR
% === Source ===
Obj.SourcePosition = [0 0 0]; % center of loudspeaker array
Obj.EmitterPosition = irs.source_position';
Obj.EmitterView = irs.head_position';
Obj.EmitterUp = [0 0 1];
% === Listener ===
% number of measurements
M = length(irs.apparent_elevation);
distance = sqrt(sum((irs.source_position-irs.head_position).^2));
Obj.ListenerPosition = irs.head_position';
[x,y,z] = sph2cart(fixnan(irs.head_azimuth'), ...
                   repmat(fixnan(irs.head_elevation'),size(irs.head_azimuth')), ...
                   repmat(distance,size(irs.head_azimuth')));
Obj.ListenerView = [x,y,z];
Obj.ListenerUp = [0 0 1];
% Receiver position for a dummy head (imported from SimpleFreeFieldHRIR)
Obj.ReceiverPosition = [0,-0.09,0; 0,0.09,0];
  
%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

end

function x = fixnan(x)
    if isnan(x), x=0; end
end
