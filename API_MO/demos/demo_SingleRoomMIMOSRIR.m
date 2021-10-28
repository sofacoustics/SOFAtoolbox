% Demonstrates the usage of SingleRoomMIMOSRIR. 
% Idea: demo_SingleRoomMIMOSRIR loads SRIRs from ???, 
% converts to SH as SingleRoomMIMOSRIR, does planewave decompositon, and stores
% the data as a ???.sofa file.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% SOFA API - demo script
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Data
% https://phaidra.kug.ac.at/view/o:104376, LS: 3 oder 5
% documentation: https://www.mdpi.com/2076-3417/10/11/3747

% Download IR dataset from https://phaidra.kug.ac.at/view/o:104376
% Copy the folder 'LigetiHall_CubeletToSt450' to your current working directory

Obj = SOFAgetConventions('SingleRoomMIMOSRIR');
Obj.GLOBAL_Title = 'LigetiHall_CubeletToSt450';
Obj.GLOBAL_Organization = 'Institut f�r elektronische musik und akustik, Graz, Austria';
Obj.GLOBAL_AuthorContact = 'julien.demuynke@eurecat.org, markus.zaunschirm@atmoky.com, zotter@iem.at';
Obj.GLOBAL_References = 'Auralization of High-Order Directional Sources from First-Order RIR Measurements (2020)';
Obj.GLOBAL_Comment = '1rst order RIR measurements with Cubelet loudspeaker array protoype (6 loudspeakers) and TSL ST450 microphone (4 channels) in the György Ligeti Saal, Graz, Austria. The source was surrounded by 4 reflecting baffles (0.9 x 1.8 m) in its rear halfspace.';
Obj.GLOBAL_ListenerShortName = 'TSL ST450';
Obj.GLOBAL_SourceShortName = 'Cubelet loudspeaker array protoype';
Obj.GLOBAL_DatabaseName = 'LigetiHall_CubeletToSt450_1rst_order_RIR';
Obj.GLOBAL_ListenerDescription = ' 4-channel Ambisonic B-format microphone array with r = 2 cm';
Obj.GLOBAL_SourceDescription = ' Spherical (r = 7.5 cm) 6-channel loudspeaker array prototype with loudspeakers arranged on surfaces of a cube';
Obj.GLOBAL_RoomDescription = 'György-Ligeti Room in building LG14 - MUMUTH in Universität für Musik und darstellende Kunst Graz, Austria. Shoebox shaped room of surface 511.15m2 and volume 5630m3 with floor made of wooden panels';

Obj.ReceiverPosition_Type = 'spherical harmonics';%B-format
Obj.EmitterPosition = [0 0 0.075;90 0 0.075;180 0 0.075;270 0 0.075;0 90 0.075;0 -90 0.075];
Obj.SourcePosition = [4.2 0 0];
Obj.SourceView = [-1 0 0];%The source and the listener are facing each other

addpath('./LigetiHall_CubeletToSt450/');
IR_list = dir('./LigetiHall_CubeletToSt450/*.wav');
if isempty(IR_list)
    error([' Folder does not exist or is empty: ' cd filesep 'LigetiHall_CubeletToSt450.' newline ' Download IR dataset from https://phaidra.kug.ac.at/view/o:104376 and copy the folder ''LigetiHall_CubeletToSt450'' to your current working directory.'])
%     error([' Folder does not exist or is empty: ' mfilename('fullpath') filesep 'LigetiHall_CubeletToSt450. Download IR dataset from https://phaidra.kug.ac.at/view/o:104376 and copy the folder ''LigetiHall_CubeletToSt450'' to your current working directory.'])
end
IR_INFO = audioinfo([IR_list(1).folder filesep IR_list(1).name]);

C = 3;
I = 1;
M = 1;
R = IR_INFO.NumChannels;
N = IR_INFO.TotalSamples;
E = length(IR_list);

fs = IR_INFO.SampleRate;
Obj.Data.IR = zeros(M,R,N,E);
for i = 1:6
    IR = audioread([IR_list(i).folder filesep IR_list(i).name]);
    Obj.Data.IR(1,:,:,i) = transpose(IR);
end
Obj.Data.SamplingRate = fs;
Obj.Data.Delay = zeros(M,R,E);

Obj.ReceiverView = zeros(R,C,I);
Obj.ReceiverUp = zeros(R,C,I);
Obj.EmitterView = zeros(E,C,I);
Obj.EmitterUp = zeros(E,C,I);

Obj = SOFAupdateDimensions(Obj);
SOFAsave('LigetiHall_CubeletToSt450_IRs.sofa',Obj,0);
disp(['SOFA file saved to: ' pwd filesep 'LigetiHall_CubeletToSt450_IRs.sofa']);

%
% Paper for the algorithms:
% Park, M. and Rafaely, B. (2005). "Sound-field analysis by plane-wave decomposition using spherical microphone array,"
% JASA, 118, 3094�3103.  https://doi.org/10.1121/1.2063108



%% Story line...
% First, the coefficients pnm up to order N are calculated using Eq. 5, 

% pnm = sum j=1 to M, a_j p(omega_j) Yn m_j.


% and then the waves directional amplitude
% density at each frequency, i.e., plane-wave decomposition, is
% computed at the desired directions using Eq. 9.
%
%  b_n (kr,ka) = 4pi i^n(j_n(kr) ? ...) Eq. 7
% 
% w (omega_l) = sum n=0 to N sum m=?n to n p_nm / b_n Y (omega_l)

% to be saved: p_nm and b_n(kr,ka) or simple a. 
% 

%% Also this might help: Khaykin, D. and Rafaely, B. (2012). "Acoustic analysis by spherical 
%  microphone array processing of room impulse responses," JASA, 132, 261�270.
