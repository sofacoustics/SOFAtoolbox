% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% load HRTF in TU Berlin format and save as SOFA format

%% Define parameters
% Prefix to the files 
TUBfile = 'QU_KEMAR_anechoic_';
% Define radii to be loaded (per default 0.5, 1, 2, and 3 m are available)
radius=[0.5 1 2 3];
% radius=0.5;
% radius=3;

% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% start SOFA
SOFAstart;
f=filesep;

%% Load and convert the requested TU-Berlin files
clear Obj
for ii=1:length(radius)
	TUBfn=[SOFAdbPath f 'TU-Berlin KEMAR' f TUBfile num2str(radius(ii)) 'm.mat'];
	disp(['Loading: ' TUBfn]);
	TUB=load(TUBfn);
	Obj(ii)=SOFAconvertTUBerlin2SOFA(TUB.irs);
end

%% merge the different radii to a single SOFA object
ObjFull=Obj(1);
if length(radius)>1,
	for ii=2:length(radius)
		ObjFull=SOFAmerge(ObjFull,Obj(ii));
	end
end

%% save the object as a single SOFA file
SOFAfn=[SOFAdbPath f 'SOFA' f 'TU-Berlin_' TUBfile 'radius_' sprintf('%g_',radius) 'm.sofa'];
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, ObjFull, compression);