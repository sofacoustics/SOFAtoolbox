% demo for SingleRoomDRIR: save DRIR data from Uni Oldenburg (Office II) as
% SOFA file and plot the speaker and listener positions. 
% It uses conventions SingleRoomDRIR and RoomType 'shoebox', defined by
% RoomCornerA and RoomCornerB.

% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% Based on the bachelor thesis of Wolfgang Hrauda (2013).

% Measurement situation:
%   T_60 = 300 ms
%   speaker positions:
%     entrance, desk, desk, at window (A, B, C, D)
%   conditions: 
%     for position A and D: door and window were open
%   Listener: headorientation = {1 2} -> 0°/-90° (straight/looking over right shoulder)
%   Receivers:
%     ch 1 & 2:  'in-ear'    IRs from in-ear microphoes     (not available for Office I)
%                'bte'       6-channel BTE-IRs:
%     ch 3 & 4:  'front'     BTE-IRs front microphone pair
%     ch 5 & 6:  'middle'    BTE-IRs middle microphone pair
%     ch 7 & 8:  'rear'      BTE-IRs rear microphone pair
%   ambient noise: telepohne, keyboard typing, ventilation (all for both
%     orientations); opening and closing the door (15 times)


%% Define parameters


% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% --- loading Uni Oldenburg data ----
PathO=fullfile(fileparts(SOFAdbPath), 'Oldenburg','HRIR_database_mat','hrir','office_II');
disp(['Loading: ' PathO]);
A1 = load([PathO filesep 'office_II_1_A.mat']); % has dimension: N x R
B1 = load([PathO filesep 'office_II_1_B.mat']);
C1 = load([PathO filesep 'office_II_1_C.mat']);
D1 = load([PathO filesep 'office_II_1_D.mat']);
A2 = load([PathO filesep 'office_II_2_A.mat']);
B2 = load([PathO filesep 'office_II_2_B.mat']);
C2 = load([PathO filesep 'office_II_2_C.mat']);
D2 = load([PathO filesep 'office_II_2_D.mat']);

%% Get an empy conventions structure
Obj = SOFAgetConventions('SingleRoomDRIR');

%% Listener and Receiver
Obj.ReceiverPosition = [0 0.09 0; 0 -0.09 9; ... % in-ear
                    0.02 0.09 0.02; 0.02 -0.09 0.02; ... % front bte
                    0.02 0.09 0.02; 0.02 -0.09 0.02; ... % middle bte
                    0.02 0.09 0.02; 0.02 -0.09 0.02]; % rear bte
Obj.ListenerPosition = [2.16 4.4 1.1];
Obj.ListenerUp = [0 0 1];
Obj.ListenerView = [repmat([1 0 0],4,1); repmat([0 -1 0],4,1)];

%% Source and Transmitter

Obj.EmitterPosition = [0 0 0];
A = [0.52 5.27 1.8];
B = [0.79 1.25 1.1];
C = [2.52 1.23 1.1];
D = [2.38 0 1.8];
Obj.SourcePosition = [A; B; C; D; A; B; C; D;];
Av = [-1 0 0];
Bv = [1 0 0];
Cv = [1 0 0];
Dv = [1 0 0];
Obj.SourceView = [Av; Bv; Cv; Dv; Av; Bv; Cv; Dv;];
Obj.SourceUp = [0 0 1];

%% Fill Data with data
% change data matrix dimension: N x R -> M x R x N
A1.data = shiftdim(A1.data,1);
data(1,:,:) = A1.data;
B1.data = shiftdim(B1.data,1);
data(2,:,:) = B1.data;
C1.data = shiftdim(C1.data,1);
data(3,:,:) = C1.data;
D1.data = shiftdim(D1.data,1);
data(4,:,:) = D1.data;
A2.data = shiftdim(A2.data,1);
data(5,:,:) = A2.data;
B2.data = shiftdim(B2.data,1);
data(6,:,:) = B2.data;
C2.data = shiftdim(C2.data,1);
data(7,:,:) = C2.data;
D2.data = shiftdim(D2.data,1);
data(8,:,:) = D1.data;
Obj.Data.IR = data; 
Obj.Data.Delay = zeros(1,size(Obj.Data.IR,2));
Obj.Data.SamplingRate = A1.fs;

%% Fill with attributes
Obj.GLOBAL_ListenerShortName = 'HATS';
Obj.GLOBAL_History='Converted from the Uni Oldenburg database';
Obj.GLOBAL_License='http://medi.uni-oldenburg.de/hrir/html/download.html';
Obj.GLOBAL_References='H. Kayser, S. D. Ewert, J. Anemüller, T. Rohdenburg, V. Hohmann, and B. Kollmeier, "Database of Multichannel In-Ear and Behind-the-Ear Head-Related and Binaural Room Impulse Responses," EURASIP Journal on Advances in Signal Processing, vol. 2009, doi:10.1155/2009/298605';

%% Setup the room
Obj.GLOBAL_RoomType='shoebox';
Obj.GLOBAL_RoomDescription='Office II at the University of Oldenburg, T_60 = 300 ms';

Obj.RoomCornerA = [0; 0; 0];
Obj.RoomCornerA_Type = 'cartesian';
Obj.RoomCornerA_Units = 'meter';
Obj.API.Dimensions.RoomCornerA='C';

Obj.RoomCornerB = [3.3; 6; 0];
Obj.RoomCornerB_Type = 'cartesian';
Obj.RoomCornerB_Units = 'meter';
Obj.API.Dimensions.RoomCornerB='C';
 
%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

%% Save SOFA file
SOFAfn=fullfile(SOFAdbPath,'sofa_api_mo_test','Oldenburg_OfficeII.sofa');
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj, compression);
clear Obj;

%% Re-load the file
disp(['Reloading:  ' SOFAfn]);
X=SOFAload(SOFAfn);

%% Plot the 2D-plan of the measurement setup
SOFAplotGeometry(X);
title('Office II from Kayser et al. (2009) saved as SingleRoomDRIR');
