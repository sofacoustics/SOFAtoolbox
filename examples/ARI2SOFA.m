function Obj=ARI2SOFA(fn,hM,meta,stimPar,compression)
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
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

%% Default compression
if ~exist('compression','var'), compression=1; end

%% Fill data to empty structure
	% Define the data, its type, and the data-type-specific metadata
Obj.Data.IR = shiftdim(hM,1); % hM is [N M R], data.IR must be [M R N]
Obj.Data.SamplingRate = stimPar.SamplingRate;

% 	% Define general information about the file
Obj.GLOBAL_SubjectID = stimPar.SubjectID;
Obj.GLOBAL_DatabaseName = 'ARI';
Obj.GLOBAL_ApplicationName = 'ExpSuite';
Obj.GLOBAL_ApplicationVersion = stimPar.Version;
Obj.GLOBAL_Organization = 'Acoustics Research Institute';
Obj.GLOBAL_AuthorContact = 'piotr@majdak.com';

	 % Define the Source
% Obj.TransmitterPositionType = 'cartesian';
% Obj.TransmitterPosition = [0 0 0];	% one transmitter on the source
% Obj.SourcePositionType = 'cartesian';
% Obj.SourceUpType = 'cartesian';
% Obj.SourceViewType = 'cartesian';
% Obj.SourcePosition = [0 0 0];
% Obj.SourceView = [1 0 0];
% Obj.SourceUp = [0 0 1];

	% Define receivers and the listener
% Obj.ReceiverPositionType = 'cartesian';
% Obj.ReceiverPosition = [0 -0.09 0; 0 0.09 0];
% Obj.ListenerPositionType = 'cartesian';
% Obj.ListenerViewType = 'cartesian';
% Obj.ListenerUpType = 'cartesian';
Obj.ListenerPosition = [1.2 0 0];
Obj.ListenerView = [0 0 0];
Obj.ListenerUp = [1.2 0 1];
Obj.ListenerRotation = [meta.pos(1:size(hM,2),4) meta.pos(1:size(hM,2),5) zeros(size(hM,2),1)];

% 	% Additional measurement metadata
% Obj.MeasurementID = stimPar.ID;
% Obj.MeasurementParameterSourceAudioChannel = meta.pos(:,3);
% Obj.MeasurementParameterItemIndex = meta.itemidx;
% Obj.MeasurementParameterAudioLatency = meta.lat;
% Obj.MeasurementParameterSourceAmplitude = meta.amp;

	% Room type
% Obj.GLOBAL_RoomType = 'free field';

%% Update dimensions UMRNECQ
Obj.I=1;
Obj.M=size(Obj.Data.IR,1);
Obj.R=size(Obj.Data.IR,2);
Obj.N=size(Obj.Data.IR,3);
Obj.E=size(Obj.EmitterPosition,1);
Obj.C=3;

%% write data to sofa file
SOFAsave(fn,Obj,compression); 
