function [hM,meta,stimPar]=SOFAconvertSOFA2ARI(Obj)
% OBJ=SOFAconvertSOFA2ARI(hM,meta,stimPar) converts the  HRTFs described in hM, meta, and
% stimPar (see ARI HRTF format) to a SOFA object.
%

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Get an empy conventions structure
if ~strcmp(Obj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR'),
	error('ARI Format supports only SimpleFreeFieldHRIR SOFA conventions');
end

%% Fill data matrix
hM=shiftdim(Obj.Data.IR,2);	% data.IR is [M R N], hM must be [N M R]

%% Fill stimPar
stimPar.SamplingRate = Obj.Data.SamplingRate;
stimPar.TimeBase = 1e6/stimPar.SamplingRate;
stimPar.SubjectID = Obj.GLOBAL_ListenerShortName;

%% Fill meta
	% Fill in geodesic coordinate system where azi=(0;360) and ele=(-90;90);
meta.pos(:,1)=bsxfun(@times,Obj.SourcePosition(:,1),ones(Obj.API.M,1));
meta.pos(:,2)=bsxfun(@times,Obj.SourcePosition(:,2),ones(Obj.API.M,1));
  % Fill in the Channel
if isfield(Obj,'MeasurementSourceAudioChannel'),
	meta.pos(:,3) = bsxfun(@times,Obj.MeasurementSourceAudioChannel,ones(Obj.API.M,1));
else
	meta.pos(:,3) = NaN(Obj.API.M,1);
end
if isfield(Obj,'MeasurementAudioLatency'),
	meta.lat=bsxfun(@times,Obj.MeasurementAudioLatency,ones(Obj.API.M,1));
end
	% create horizontal-polar coordinates 
[meta.pos(:,6), meta.pos(:,7)]=sph2hor(meta.pos(:,1),meta.pos(:,2));
	% create continuous-elevation coordinates where azi=(-90;90) and ele=(0;360);
meta.pos(:,4)=meta.pos(:,1);
meta.pos(:,5)=meta.pos(:,2);
idx=find(meta.pos(:,1)>90 & meta.pos(:,1)<=270);
meta.pos(idx,4)=meta.pos(idx,4)-180;
meta.pos(idx,5)=180-meta.pos(idx,5);
idx=find(meta.pos(:,1)>270 & meta.pos(:,1)<=360);
meta.pos(idx,4)=meta.pos(idx,4)-360;

%% Fill with unknown but probably mandatory data
% stimPar.Channel = 0;
% stimPar.Electrode = 0;
% stimPar.PulseNr = 0;
% stimPar.Period = 0;
% stimPar.Offset = 0;
stimPar.Resolution = 24; % assume a 24-bit ADC/DAC resolution
% stimPar.FadeIn = 0;
% stimPar.FadeOut = 0;
% stimPar.Length = 0;
% stimPar.FittFileName = '';
% stimPar.StimFileName = '';
stimPar.GenMode = 1;
% stimPar.WorkDir = '';
stimPar.ID = 'hrtf';
stimPar.Version = '2.0.0';

