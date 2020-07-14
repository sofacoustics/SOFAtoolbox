%%
% This demo creates an artificial HRTF set.
% It shows how to use SOFAgetConventions and SOFAsave
% The HRTF set contains single pulses placed at sample index of 100
% which results in a broadband delay of 100 samples. 
% Each IR is 256 samples long (i.e., N=256)

% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Set parameters
% Latency of the created IRs
latency=100; % in samples, must be 1<latency<256
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time


%% Get an empy conventions structure
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

%% Define positions -  we use the standard CIPIC positions here
lat1=[-80 -65 -55 -45:5:45 55 65 80];    % lateral angles
pol1= -45 + 5.625*(0:49);                % polar angles
pol=repmat(pol1',length(lat1),1);
lat=lat1(round(0.5:1/length(pol1):length(lat1)+0.5-1/length(pol1)));

%% Create the impulse response
N=256;
IR=[zeros(latency,1); 1; zeros(N-latency-1,1)];

%% Fill data with data
M=length(lat1)*length(pol1);
Obj.Data.IR = NaN(M,2,N); % data.IR must be [M R N]

ii=1;
for aa=1:length(lat1)
	for ee=1:length(pol1)
		Obj.Data.IR(ii,1,:)=IR;
		Obj.Data.IR(ii,2,:)=IR;
		[azi,ele]=hor2sph(lat(ii),pol(ii));
    Obj.SourcePosition(ii,:)=[azi ele 1];
		Obj.SourcePosition(ii,:)=[azi ele 1];
		ii=ii+1;
	end
end

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

%% Fill with attributes
Obj.GLOBAL_ListenerShortName = 'KEMAR';
Obj.GLOBAL_History = 'created with a script';
Obj.GLOBAL_DatabaseName = 'none';
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA API';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
Obj.GLOBAL_Organization = 'Acoustics Research Institute';
Obj.GLOBAL_AuthorContact = 'piotr@majdak.com';
Obj.GLOBAL_Comment = 'Contains simple pulses for all directions';

%% save the SOFA file
SOFAfn=fullfile(SOFAdbPath,'sofa_api_mo_test','Pulse.sofa');
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj, compression);
