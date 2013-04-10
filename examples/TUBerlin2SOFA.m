function TUBerlin2SOFA(outfile,irs,compression)
% SOFA API - demo script
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% $LastChangedDate: $
% $LastChangedRevision: $
% $LastChangedBy: $

%% TU to Sofa format conversion
% The files can be downloaded at
% https://dev.qu.tu-berlin.de/projects/measurements/repository/show/2010-11-kemar-anechoic/mat

% === Configuration ===
%infile = 'QU_KEMAR_anechoic_3m.mat';
%outfile = 'QU_KEMAR_anechoic_3m';

% Load data set
%load(infile);

%% Convention
HRIR = SOFAgetConventions('SimpleFreeFieldHRIR');

% Get data matrix
data(:,:,1) = irs.left;
data(:,:,2) = irs.right;
% => dimension of M is [N M R], we need [M R N], where N are the sampling
% points, R the number of receiver channels, M the number of measurements
data = shiftdim(data,1);


%% ===== DATA ============================================================
% Store data
HRIR.Data.IR = data;
HRIR.Data.SamplingRate = irs.fs;
% Update dimensions
HRIR.M = size(data,1);
HRIR.N = size(data,3);


%% ===== METADATA ========================================================
% General information
HRIR.GLOBAL_SubjectID = irs.head;
HRIR.GLOBAL_DatabaseName = 'QU/AIPA, TU Berlin'; % maybe setting the name by function parameter
HRIR.GLOBAL_ApplicationName = 'HRIR';
HRIR.GLOBAL_ApplicationVersion = '1.0';
HRIR.GLOBAL_Organization = 'Technische Universität Berlin';
HRIR.GLOBAL_AuthorContact = 'hagen.wierstorf@tu-berlin.de';
% Receiver and listener
HRIR.SourcePosition = irs.source_position';
HRIR.SourceView = irs.source_reference';
HRIR.SourceUp = irs.source_position' + [0 0 1];
HRIR.ListenerPosition = irs.head_position';
HRIR.ListenerView = irs.head_reference';
HRIR.ListenerUp = irs.head_position + [0 0 1]';
HRIR.ListenerRotation = [-irs.apparent_azimuth' -irs.apparent_elevation' ...
    zeros(length(irs.apparent_azimuth),1)];

% write data to SOFA file
SOFAsave(outfile,HRIR,compression);

%% some very basic examples who the functions might be used:
%% load the entire data
%results1 = SOFAload(Filename);
%% get a list of SubjectIDs and corresponding ListenerRotations
%results2 = SOFAlistVariables(Filename,'SubjectID','ListenerRotation');
%% get all positions within a range of 80 to 100 degrees azimuth
%results3 = SOFAgetID({Filename},'ListenerRotation',[90 0 0],'=',{'TargetValueRange',10});
%% get data set for each position that was found by getID (previous line)
%for ii=1:size(results3)
%  results4{ii} = SOFAgetData(Filename,results3(ii));
%end
