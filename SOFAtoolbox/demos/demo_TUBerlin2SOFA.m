%demo_TUBerlin2SOFA - Load HRTF in TU Berlin format and save as SOFA format.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: bugs fixed (10.2021)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)

% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Define parameters
% Prefix to the files 
TUBfile = 'QU_KEMAR_anechoic_';
% Define vector with radii to be loaded. Available files: 0.5, 1, 2, and 3 m
% radius=[0.5 1 2 3];
radius=0.5;

% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% Load, convert, and save the requested TU-Berlin files
for ii=1:length(radius)
		% load

	TUBfn=fullfile(fileparts(SOFAdbPath), 'TU-Berlin KEMAR', [TUBfile num2str(radius(ii)) 'm.mat']);
    
    if isfile(TUBfn)
        disp(['Loading: ' TUBfn]);
        TUB=load(TUBfn);
    else
        if radius(ii) == 0.5 % catch if file name ends with "05m.mat" instead of "0.5m.mat"
            TUBfn2=fullfile(fileparts(SOFAdbPath), 'TU-Berlin KEMAR', [TUBfile '05m.mat']);
            if isfile(TUBfn2)
                disp(['Loading: ' TUBfn2]);
                TUB=load(TUBfn2);
            else
                warning(['File not existing: ' TUBfn '  -->  Please download it to: ' fullfile(fileparts(SOFAdbPath), 'TU-Berlin KEMAR')]);
                error(['Sorry.... ' mfilename ' cannot complete!']);
            end
        else
            warning(['File not existing: ' TUBfn '  -->  Please download it to: ' fullfile(fileparts(SOFAdbPath), 'TU-Berlin KEMAR')]);
            error(['Sorry.... ' mfilename ' cannot complete!']);
        end
    end

		% convert and add application specific metadata
	Obj=SOFAconvertTUBerlin2SOFA(TUB.irs);
	Obj.GLOBAL_DatabaseName = 'TU-Berlin'; % maybe setting the name by function parameter
	Obj.GLOBAL_ApplicationName = 'Demo of the SOFA Toolbox';
	Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
	Obj.GLOBAL_Organization = 'Technische Universität Berlin';
	Obj.GLOBAL_AuthorContact = 'hagen.wierstorf@tu-berlin.de';
		% save
	SOFAfn=fullfile(SOFAdbPath, 'sofatoolbox_test', ['TU-Berlin_' TUBfile 'radius_' sprintf('%g',radius(ii)) 'm.sofa']);
	disp(['Saving:  ' SOFAfn]);
	Obj=SOFAsave(SOFAfn, Obj, compression);
end
