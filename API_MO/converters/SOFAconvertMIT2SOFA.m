function Obj=SOFAconvertMIT2SOFA(root,pinna)
% OBJ=SOFAconvertMIT2SOFA(root,pinna) loads the MIT HRTFs saved in a 
% directory ROOT for the PINNA and converts to a SOFA object.
% PINNA must be 'normal' or 'large'.

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Create correct parts of the file name
%	'L' use full data from left pinna (normal pinna)
%	'R' use full data from right pinna (large red pinna)
%	'H' use compact data -> not supported now
switch pinna
	case 'normal'
		postfix='L'; idx=[1 2];
		prefix='full';
	case 'large'
		postfix='R'; idx=[2 1];
		prefix='full';
	otherwise
		error('Pinna not supported');
end

%% Get an empy conventions structure
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

%% Define elevations
eles = [-40 -30 -20 -10  0 10 20 30 40 50 60 70 80 90];
elecnt=[ 56  60  72  72 72 72 72 60 56 45 36 24 12  1];

%% Determine data size
M=sum(elecnt);
Obj.SourcePosition=zeros(M,3);
Obj.Data.IR=zeros(M,2,length(wavread([root filesep prefix filesep 'elev0' filesep postfix '0e000a.wav'])));

%% Fill with data 
Obj.Data.SamplingRate = 44100;
ii=1;
for ei = 1 : length(eles)
	ele = eles(ei);
	for ai = 0 : elecnt(ei)-1
		azi = 360/elecnt(ei)*ai;
		flip_azi = mod(360-azi,360);
		fn=[root filesep prefix filesep 'elev' num2str(ele) filesep postfix num2str(ele) 'e' sprintf('%03d',round(flip_azi)) 'a.wav'];
		Obj.Data.IR(ii,idx(1),:) = wavread(fn)'; % data.IR must be [M R N]
    dirfn=dir(fn);
		fn=[root filesep prefix filesep 'elev' num2str(ele) filesep postfix num2str(ele) 'e' sprintf('%03d',round(azi)) 'a.wav'];
		Obj.Data.IR(ii,idx(2),:) = wavread(fn)';
      % SimpleFreeFieldHRIR 0.2
        % Obj.ListenerRotation(ii,:)=[azi ele 0];
      % SimpleFreeFieldHRIR 0.3
    Obj.SourcePosition(ii,:) = [azi ele 1.4];
    ii=ii+1;    
	end
end

%% Fill with attributes
Obj.GLOBAL_ListenerShortName = ['KEMAR, ' pinna ' pinna'];
Obj.GLOBAL_History='Converted from the MIT format';
Obj.GLOBAL_DateCreated = datestr(datenum(dirfn.date),SOFAdefinitions('dateFormat'));

%% Fill the mandatory variables
  % SimpleFreeFieldHRIR 0.2
%   Obj.ListenerPosition = [1.4 0 0];
%   Obj.ListenerView = [-1 0 0];
%   Obj.ListenerUp = [0 0 1];
  % SimpleFreeFieldHRIR 0.3
Obj.ListenerPosition = [0 0 0];
Obj.ListenerView = [1 0 0];
Obj.ListenerUp = [0 0 1];

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);
