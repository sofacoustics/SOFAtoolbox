function Obj = SOFAconvertBTdei2SOFA(BTdei)
% OBJ=SOFAconvertBTdei2SOFA(BTdei) converts the HRTFs described in BTdei
% to a SOFA object.
% 
% BTdei format is used by Michele Geronazzo, University of Padova.
%

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


% get empty object (the flag 'm' provides the mandatory fields only)
Obj=SOFAgetConventions('SimpleHeadphoneIR','m');

Obj.GLOBAL_Title = 'HPIR';
Obj.GLOBAL_DatabaseName = BTdei.specs.Database;

Obj.GLOBAL_History =             'Converted from the BTdei format';
Obj.GLOBAL_License =             'Creative Commons Attribution-NonCommercial-ShareAlike 3.0';
Obj.GLOBAL_ApplicationName =     'HpTFs from DEI - University of Padova';
Obj.GLOBAL_ApplicationVersion =  SOFAgetVersion('API');
Obj.GLOBAL_AuthorContact =       'geronazzo@dei.unipd.it';
Obj.GLOBAL_References = ['M. Geronazzo, F. Granza, S. Spagnol, F. Avanzini. ', ...
    'A Standardized Repository of Head-Related and Headphone Impulse Response Data ', ...
    'In 134th Convention of the Audio Engineering Society, May 2013.'];
Obj.GLOBAL_Organization = 'Department of Information Engineering, University of Padova';
Obj.GLOBAL_Comment = '';

% copy the data

% Emitter - Source
Obj.GLOBAL_SourceDescription    = [BTdei.hp.Id ' - ' BTdei.hp.Producer ' ' BTdei.hp.Model];
Obj.GLOBAL_SourceManufacturer   = BTdei.hp.Producer;
Obj.GLOBAL_SourceModel          = BTdei.hp.Model;
%Obj.GLOBAL_SourceURI           = BTdei.hp.Uri;
% Receiver - Listener
Obj.GLOBAL_SubjectID            = BTdei.specs.SubjectId;
Obj.GLOBAL_ListenerDescription  = BTdei.sbjType;
Obj.GLOBAL_ReceiverDescription  = BTdei.specs.MicrophonePosition; % qualitative data e.g. blocked ear canal, open ear canal, at the eardrum


Obj.ListenerPosition = [0 0 0];
Obj.ReceiverPosition = [0 0.09 0; 0 -0.09 0];
Obj.SourcePosition   = [0 0 0];
Obj.EmitterPosition  = [0 0.09 0; 0 -0.09 0];

%% Fill data with data
Obj.Data.SamplingRate           = BTdei.specs.SampleRate; % Sampling rate

% calculate the effective size of the data matrix
M = length(BTdei.data);         % number of repositionings
R = size(BTdei.data(1).HpIR,2); % number of channels (stereo)

len_vec = zeros(M,1);
for ii=1:M
    len_vec(ii)= length(BTdei.data(ii).HpIR);
end
N = max(len_vec);

Obj.API.M=M;
Obj.API.R=R;
Obj.API.N=N;

% store IR data
Obj.Data.IR = NaN(M,R,N); % data.IR must be [M R N]
for aa=1:M
	  HpIR = [BTdei.data(aa).HpIR; zeros((N-length(BTdei.data(aa).HpIR)),2)];
    Obj.Data.IR(aa,1,:)= HpIR(:,1)';
		Obj.Data.IR(aa,2,:)= HpIR(:,2)';	
end


% update dimensions
Obj = SOFAupdateDimensions(Obj);
end

