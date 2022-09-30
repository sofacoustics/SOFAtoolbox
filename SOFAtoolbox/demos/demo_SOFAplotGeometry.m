%demo_SOFAplotGeometry - Script demonstrating the usage of SOFAplotGeometry.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% load a SOFA file in SimpleFreeFieldHRIR

SOFAfile=fullfile(SOFAdbPath,'database','widespread','ICO1m_00139.sofa');
disp(['Loading: ' SOFAfile]);
Obj=SOFAload(SOFAfile);

% plot all measurements
SOFAplotGeometry(Obj);
title(['Geometry SimpleFreeFieldHRIR, ' num2str(Obj.API.M) ' position(s)'])
set(gcf, 'Name', mfilename);

% % only show every 45th measurement
index = 1:45:Obj.API.M;
SOFAplotGeometry(Obj,index);
title(['Geometry SimpleFreeFieldHRIR, reduced to ' num2str(size(index,2)) ' position(s)'])
set(gcf, 'Name', mfilename);

% % %% load a SingleRoomDRIR SOFA file (outdated)
% disp(['Loading: ' 'db://' fullfile('database','thk','DRIR_LBS_VSA_1202RS_SBL.sofa')]);
% Obj=SOFAload(['db://' ...
%   fullfile('database','thk','DRIR_LBS_VSA_1202RS_SBL.sofa')]);
% 
% % plot SOFA Object with 1202 Receivers
% SOFAplotGeometry(Obj);
% set(gcf, 'Name', mfilename)

% % remove all but one Receiver
Obj.ReceiverPosition = [0 0.09 0];
Obj.ReceiverPosition_Type = 'cartesian';
Obj.Data.IR = Obj.Data.IR(:,1,:);
Obj.Data.Delay = Obj.Data.Delay(:,1,:);
Obj = SOFAupdateDimensions(Obj);

SOFAplotGeometry(Obj);
title(['Geometry SimpleFreeFieldHRIR, ' num2str(Obj.API.R) ' receiver(s), ' num2str(Obj.API.M) ' position(s)'])
set(gcf, 'Name', mfilename);

% %% load a GeneralFIR SOFA file
SOFAfile=fullfile(SOFAdbPath,'database', 'tu-berlin','FABIAN_CTF_modeled.sofa');
Obj = SOFAload(SOFAfile);

SOFAplotGeometry(Obj);
title(['Geometry GeneralFIR, ' num2str(Obj.API.R) ' receiver(s), ' num2str(Obj.API.M) ' position(s)'])
set(gcf, 'Name', mfilename);

% % %% load example with room geometry (outdated)
% disp(['Loading: ' SOFAfile]);
% SOFAfile = fullfile(SOFAdbPath,'sofatoolbox_test', 'Oldenburg_OfficeII.sofa');
% Obj = SOFAload(SOFAfile);
% 
% SOFAplotGeometry(Obj);
% set(gcf, 'Name', mfilename)

% %% if exists try plotting SOFA file containing spherical harmonic emitter
% if exist(fullfile(SOFAdbPath,'demo_SHforHRTFs_SH.sofa'))
%    Obj =  SOFAload(fullfile(SOFAdbPath,'demo_SHforHRTFs_SH.sofa'));
%     SOFAplotGeometry(Obj);
% else
%     error('Run demoSHforHRTFS.m first.')
% end