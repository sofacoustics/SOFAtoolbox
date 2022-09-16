%demo_GeneralFIR_E - Demonstrates the usage of the GeneralFIR-E conventions.

% #Author: Michael Mihocic
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Get an empy conventions structure
clear;
conventions='GeneralFIR-E';
disp(['Creating SOFA file with ' conventions 'conventions...']);
Obj = SOFAgetConventions(conventions);

%% Fill some data...
% Define positions -  we use the standard CIPIC positions here
lat1=[-80 -65 -55 -45:5:45 55 65 80];    % lateral angles
pol1= -45 + 5.625*(0:49);                % polar angles
pol=repmat(pol1',length(lat1),1);
lat=lat1(round(0.5:1/length(pol1):length(lat1)+0.5-1/length(pol1)));

% Create the impulse response
N=256;
IR=[zeros(100,1); 1; zeros(N-100-1,1)];

% Fill data with data
M=length(lat1)*length(pol1);
Obj.Data.IR = NaN(M,2,N); % data.IR must be [M R N]
Obj.Data.Delay=[0 0];

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
Obj.GLOBAL_ListenerShortName = 'dummy';
Obj.GLOBAL_History = 'created with a demo script';
Obj.GLOBAL_DatabaseName = 'none';
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA Toolbox';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
Obj.GLOBAL_Organization = 'Acoustics Research Institute';
Obj.GLOBAL_AuthorContact = 'michael.mihocic@oeaw.ac.at';
Obj.GLOBAL_Comment = 'Contains simple pulses for all directions';

%% save the SOFA file
SOFAfn=fullfile(SOFAdbPath,'sofatoolbox_test',[strrep(conventions,'-','_') '.sofa']);
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj);
