function Obj=SOFAconvertCIPIC2SOFA(CIPIC)
% Obj=CIPIC2SOFA(CIPIC) converts the HRTFs described in the structure CIPIC
% (see CIPIC HRTF format) to a SOFA object.
%

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 


%% Get an empy conventions structure
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

%% Define positions
lat1=[-80 -65 -55 -45:5:45 55 65 80];    % lateral angles
pol1= -45 + 5.625*(0:49);                % polar angles
pol=repmat(pol1',length(lat1),1);
ida=round(0.5:1/length(pol1):length(lat1)+0.5-1/length(pol1));
lat=lat1(ida);

%% Fill data with data
M=length(lat1)*length(pol1);
Obj.Data.IR = NaN(M,2,size(CIPIC.hrir_l,3)); % data.IR must be [M R N]
Obj.Data.SamplingRate = 44100;

ii=1;
for aa=1:length(lat1)
	for ee=1:length(pol1)
		Obj.Data.IR(ii,1,:)=CIPIC.hrir_l(aa,ee,:);
		Obj.Data.IR(ii,2,:)=CIPIC.hrir_r(aa,ee,:);
		[azi,ele]=hor2sph(lat(ii),pol(ii));
		Obj.ListenerRotation(ii,:)=[azi ele 0];
		ii=ii+1;
	end
end

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

%% Fill with attributes
Obj.GLOBAL_SubjectID = CIPIC.name;
Obj.GLOBAL_DatabaseName = 'CIPIC';
Obj.GLOBAL_ApplicationName = 'CIPIC2SOFA';
Obj.GLOBAL_ApplicationVersion = '';
Obj.GLOBAL_Organization = 'Acoustics Research Institute';
Obj.GLOBAL_AuthorContact = 'piotr@majdak.com';
Obj.GLOBAL_History = 'converted from the CIPIC database';

%% Fill the mandatory variables
Obj.ListenerPosition = [1 0 0];
Obj.ListenerView = [0 0 0];
Obj.ListenerUp = [1 0 1];
