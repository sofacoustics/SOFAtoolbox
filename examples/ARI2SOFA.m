function ARI2SOFA(fn,hM,meta,stimPar,compression)
% ARI2SOFA(FN,hM,meta,stimPar) saves the HRTFs describen in hM, meta, and
% stimPar (see ARI HRTF format) in a SOFA file with the filename FN.
%
% ARI2SOFA(FN,hM,meta,stimPar,C) uses C as compression (0..uncompressed,
% 9..most compressed). Default=1, results in a nice compression within a
% reasonable processing time.
%

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% created by Wolfgang Hrauda (summer 2012)
% Piotr Majdak (script->function, generalizations, optimizations)
% Last Updates: Michael Mihocic, 04.2013

%% Convention
Obj = SOFAgetConventions('SimpleFreeFieldHRTF');
% Obj = SOFAgetConventions('SimpleDRIRMicArray');

%% Default compression
if ~exist('compression','var'), compression=1; end
data.FIR=shiftdim(hM,1); % hM is [N M R], data.FIR must be [M R N]

%% Fill data to empty structure
	% Define the data, its type, and the data-type-specific metadata
Obj.Data.IR = data.FIR;
% a=1;
Obj.DataType = 'IR';
Obj.SamplingRate = stimPar.SamplingRate;

% 	% Define general information about the file
Obj.SubjectID = stimPar.SubjectID;
Obj.DatabaseName = 'ARI';
Obj.ApplicationName = 'ExpSuite';
Obj.ApplicationVersion = stimPar.Version;
Obj.Organization = 'Acoustics Research Institute';
Obj.AuthorContact = 'michael.mihocic@oeaw.ac.at';

	% Define the Source
Obj.TransmitterPositionType = 'cartesian';
Obj.TransmitterPosition = [0 0 0];	% one transmitter on the source
Obj.SourcePositionType = 'cartesian';
Obj.SourceUpType = 'cartesian';
Obj.SourceViewType = 'cartesian';
Obj.SourcePosition = [0 0 0];
Obj.SourceView = [1 0 0];
Obj.SourceUp = [0 0 1];

	% Define the Listener
%     Obj.R=2;
%     Obj.E=1;
%     Obj.M='unlimited';
%     Obj.C=3;
Obj.ReceiverPositionType = 'cartesian';
Obj.ReceiverPosition = [0 -0.09 0; 0 0.09 0];
Obj.ListenerPositionType = 'cartesian';
Obj.ListenerViewType = 'cartesian';
Obj.ListenerUpType = 'cartesian';
Obj.ListenerPosition = [1.2 0 0];
Obj.ListenerView = [0 0 0];
Obj.ListenerUp = [1.2 0 1];
Obj.ListenerRotation = [meta.pos(:,4) meta.pos(:,5) zeros(size(meta.pos(:,1)))];

% 	% Additional measurement metadata
% Obj.MeasurementID = stimPar.ID;
% Obj.MeasurementParameterSourceAudioChannel = meta.pos(:,3);
% Obj.MeasurementParameterItemIndex = meta.itemidx;
% Obj.MeasurementParameterAudioLatency = meta.lat;
% Obj.MeasurementParameterSourceAmplitude = meta.amp;

	% Room type
% Obj.RoomType = 'free field';

%% write data to sofa file
SOFAsave(fn,Obj,compression); 
