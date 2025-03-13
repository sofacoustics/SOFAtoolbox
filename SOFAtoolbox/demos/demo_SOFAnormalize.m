%demo_SOFAnormalize - Demo script applying normalization to HRTF SOFA files, using SOFAnormalize.

% #Author: Helene Bahu: creation of test script (11.03.2025)
% #Author: Michael Mihocic: adaptions to SOFA Toolbox, header's format updated (12.03.2025)
% 
% SOFA Toolbox - demo script
% Copyright (C) Helene Bahu, helenebahu(at)gmail.com; Michael Mihocic, Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Clean up (optional)
% close all;
% clear;
% clear all % 'clear all' can break SOFA functionality
% clc

%% Adapt system parameters here (optional)
% addpath to SOFAstart
% addpath(genpath('D:\Projects\SOFA\Github\SOFAtoolbox-development\SOFAtoolbox'))
% addpath(genpath('/Users/bahu/Documents/MATLAB/SOFA Toolbox/SOFAtoolbox/'))
SofaFiles = {...
    'clubfritz/ClubFritz1.sofa', ...
    'clubfritz/ClubFritz2.sofa',...
    'clubfritz/ClubFritz3.sofa',...
    'clubfritz/ClubFritz4.sofa',...
    'clubfritz/ClubFritz5.sofa',...
    'clubfritz/ClubFritz6.sofa',...
    'clubfritz/ClubFritz7.sofa',...
    'clubfritz/ClubFritz8.sofa',...
    'clubfritz/ClubFritz9.sofa',...
    'clubfritz/ClubFritz10.sofa',...
    'clubfritz/ClubFritz11.sofa',...
    'clubfritz/ClubFritz12.sofa',...
    'bili (dtf)/IRC_1131_C_HRIR_96000.sofa',...
    % 'ari/dtf b_nh6.sofa',... % caused an error
    }; 


%% Process data
% SOFAstart;
% Define a list of file names

% Loop through all SOFA files
for i = 1:length(SofaFiles)
    SofaFile = SofaFiles{i};
    % SOFAload downloads files from SOFA Conventions repository.
    Obj = SOFAload(['db://database/' SofaFile]);
    % plot original files
    figure('Name',SofaFile);
    subplot(2,2,1);
    SOFAplotHRTF(Obj,'ETCHorizontal',1);
    title('ETC Horizontal')
    % figure('Name',SofaFile);
    subplot(2,2,2);
    SOFAplotHRTF(Obj,'MagMedian',2);
    title('Mag Median')

    % To modify the normalization parameters, use the following, and add param_S as a 2nd input parameter in HRTF_normalization_v1
    % param_S.do_gain_norm_b = 1;
    % param_S.do_resamp_b    = 0;
    % param_S.do_lp_b        = 0;
    % param_S.do_talign_b    = 0;
    % param_S.do_win_b       = 0;
    % param_S.do_zp_b        = 0;
    % param_S.do_eq_b        = 0;
    % param_S.do_LFext_b     = 0;
    % param_S.do_dist_b      = 0;

    % Apply normalization
    Objnorm_S = SOFAnormalize(Obj);
    % plot normalized data
    % figure('Name',[SofaFile ' (normalized)']);
    subplot(2,2,3);
    SOFAplotHRTF(Objnorm_S,'ETCHorizontal',1);
    title('ETC Horizontal, normalized')
    % figure('Name',[SofaFile ' (normalized)']);
    subplot(2,2,4);
    SOFAplotHRTF(Objnorm_S,'MagMedian',2);
     title('Mag Median, normalized')
end

