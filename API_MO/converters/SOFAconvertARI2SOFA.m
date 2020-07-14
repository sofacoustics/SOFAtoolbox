function Obj=SOFAconvertARI2SOFA(hM,meta,stimPar)
% OBJ=SOFAconvertARI2SOFA(hM,meta,stimPar) converts the HRTFs described in hM, meta, and
% stimPar (see ARI HRTF format) to a SOFA object.
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
Obj.Data.IR = shiftdim(hM,1); % hM is [N M R], data.IR must be [M R N]
Obj.Data.SamplingRate = stimPar.SamplingRate;

%% Fill with attributes
if isfield(stimPar, 'SubjectID'), Obj.GLOBAL_ListenerShortName = stimPar.SubjectID; end
if isfield(stimPar,'Application')
    if isfield(stimPar.Application,'Name'), Obj.GLOBAL_ApplicationName = stimPar.Application.Name; end
    if isfield(stimPar.Application,'Version'), Obj.GLOBAL_ApplicationVersion = stimPar.Application.Version; end
end

%% Fill the mandatory variables
  % SimpleFreeFieldHRIR 0.2
    % Obj.ListenerPosition = [1.2 0 0];
    % Obj.ListenerView = [-1 0 0];
    % Obj.ListenerUp = [0 0 1];
    % Obj.ListenerRotation = [meta.pos(1:size(hM,2),1) meta.pos(1:size(hM,2),2) zeros(size(hM,2),1)];
  % SimpleFreeFieldHRIR 0.3
Obj.ListenerPosition = [0 0 0];
Obj.ListenerView = [1 0 0];
Obj.ListenerUp = [0 0 1];
Obj.SourcePosition = [meta.pos(1:size(hM,2),1) meta.pos(1:size(hM,2),2) 1.2*ones(size(hM,2),1)];

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

%% Fill with some additional data
Obj.GLOBAL_History='Converted from the ARI format';
if size(meta.pos,2)>2, Obj=SOFAaddVariable(Obj,'MeasurementSourceAudioChannel','M',meta.pos(1:size(hM,2),3)); end
if isfield(meta,'lat'), Obj=SOFAaddVariable(Obj,'MeasurementAudioLatency','MR',meta.lat); end
