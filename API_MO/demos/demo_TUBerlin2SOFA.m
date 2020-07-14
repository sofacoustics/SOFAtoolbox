% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% load HRTF in TU Berlin format and save as SOFA format

%% Define parameters
% Prefix to the files 
TUBfile = 'QU_KEMAR_anechoic_';
% Define vector with radii to be loaded. Available files: 0.5, 1, 2, and 3 m
% radius=[0.5 1 2 3];
radius=[0.5];

% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% Load, convert, and save the requested TU-Berlin files
for ii=1:length(radius)
		% load
	TUBfn=fullfile(fileparts(SOFAdbPath), 'TU-Berlin KEMAR', [TUBfile num2str(radius(ii)) 'm.mat']);
	disp(['Loading: ' TUBfn]);
	TUB=load(TUBfn);
		% convert and add application specific metadata
	Obj=SOFAconvertTUBerlin2SOFA(TUB.irs);
	Obj.GLOBAL_DatabaseName = 'TU-Berlin'; % maybe setting the name by function parameter
	Obj.GLOBAL_ApplicationName = 'Demo of the SOFA API';
	Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
	Obj.GLOBAL_Organization = 'Technische Universität Berlin';
	Obj.GLOBAL_AuthorContact = 'hagen.wierstorf@tu-berlin.de';
		% save
	SOFAfn=fullfile(SOFAdbPath, 'sofa_api_mo_test', ['TU-Berlin_' TUBfile 'radius_' sprintf('%g',radius(ii)) 'm.sofa']);
	disp(['Saving:  ' SOFAfn]);
	Obj=SOFAsave(SOFAfn, Obj, compression);
end
