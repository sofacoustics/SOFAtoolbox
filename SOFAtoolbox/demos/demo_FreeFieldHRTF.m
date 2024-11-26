%demo_FreeFieldHRTF - Demonstrates the usage of spherical harmonics (SH) for HRTF interpolation. 
% 
% demo_FreeFieldHRTF loads an HRTF set, transforms to TF, then to SH, then 
% samples the horizontal and median plane in steps of 0.5 degrees. Finally, 
% the files are saved in SOFAdbPath as demo_FreeFieldHRTF_{TF, SH, TFrec}.sofa.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: save figures as optional parameter added, figures are saved with respective titles as names (10.11.2021)
% #Author: Michael Mihocic: minor bugs fixed (28.12.2021)
% #Author: Piotr Majdak: fix the propagation of global metadata, using Parents now (26.11.2024)

% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

savefigures=0; % save all created figures as fig and png (0: no, 1: yes)
if savefigures==1
    close all; % clean-up first
    warning('off', 'MATLAB:MKDIR:DirectoryExists'); % suppress warning if folder exists
    mkdir ([mfilename('fullpath') ' figures']); % output folder for figures
    folder=[mfilename('fullpath') ' figures' filesep];
end

%% Let's start, load a SimpleFreeFieldHRIR SOFA object
filename = 'database/thk/HRIR_L2354.sofa';
IR=SOFAload(['db://' filename]);
fs=IR.Data.SamplingRate;
IR.GLOBAL_APIVersion=SOFAgetVersion;
%% Figures
figure('Name',mfilename);
SOFAplotHRTF(IR,'magmedian'); 
tit='SimpleFreeFieldHRIR (FIR, mag), for reference'; % title
title(tit);
if savefigures==1
    saveas(gcf,[folder tit '.fig']); saveas(gcf,[folder tit '.png']); % save figures
%    print([folder tit '.png'], '-dpng');
end

SOFAsave(fullfile(SOFAdbPath,'sofatoolbox_test','demo_FreeFieldHRTF_1_IR.sofa'),IR);

%% Convert to TF
  % create a new SOFA structure and copy/update the global metadata
TF=SOFAgetConventions('SimpleFreeFieldHRTF');
TF.GLOBAL_AuthorContact = 'piotr.majdak@oeaw.ac.at;michael.mihocic@oeaw.ac.at';
TF.GLOBAL_Comment = 'Demonstration of the usage of the FreeFieldHRTF convention';
TF.GLOBAL_History = SOFAappendText(IR, 'GLOBAL_History', 'Converted to TF');
TF.GLOBAL_License = IR.GLOBAL_License;
TF.GLOBAL_Organization = 'ARI/ÖAW';
TF.GLOBAL_Origin = SOFAappendText(IR, 'GLOBAL_Origin', 'SOFA Toolbox');
TF.GLOBAL_Title = 'demo_FreeFieldHRTF';
TF.GLOBAL_DatabaseName = 'SOFA Toolbox Demos'; 
TF.GLOBAL_ListenerShortName = IR.GLOBAL_ListenerShortName;
TF.GLOBAL_Author = 'Piotr Majdak and Michael Mihocic';
TF.GLOBAL_ListenerDescription = IR.GLOBAL_ListenerDescription;
TF.GLOBAL_ReceiverDescription = IR.GLOBAL_ReceiverDescription;
TF.GLOBAL_SourceDescription = IR.GLOBAL_SourceDescription;
TF.GLOBAL_EmitterDescription = IR.GLOBAL_EmitterDescription;
TF.GLOBAL_RoomDescription = IR.GLOBAL_RoomDescription;
TF.GLOBAL_Parents = [SOFAdbURL '/' filename];
  % copy/update variables
TF.ListenerPosition=IR.ListenerPosition;
TF.ListenerPosition_Type=IR.ListenerPosition_Type;
TF.ListenerPosition_Units=IR.ListenerPosition_Units;
TF.ListenerView=IR.ListenerView;
TF.ListenerView_Type=IR.ListenerView_Type;
TF.ListenerView_Units=IR.ListenerView_Units;
TF.ListenerUp=IR.ListenerUp;
TF.SourcePosition=IR.SourcePosition;
TF.SourcePosition_Type=IR.SourcePosition_Type;
TF.SourcePosition_Units=IR.SourcePosition_Units;
TF.EmitterPosition=IR.EmitterPosition;
TF.EmitterPosition_Type=IR.EmitterPosition_Type;
TF.EmitterPosition_Units=IR.EmitterPosition_Units;
TF.ReceiverPosition=IR.ReceiverPosition;
TF.ReceiverPosition_Type=IR.ReceiverPosition_Type;
TF.ReceiverPosition_Units=IR.ReceiverPosition_Units;
  % allocate data
TF.Data.Real=zeros(IR.API.M,IR.API.R,IR.API.N+1);
TF.Data.Imag=zeros(IR.API.M,IR.API.R,IR.API.N+1);
  % convert to frequency domain
for ii=1:IR.API.M
  for jj=1:IR.API.R
   sp=fft(squeeze(IR.Data.IR(ii,jj,:)),2*IR.API.N); % Delay not considered!
   TF.Data.Real(ii,jj,:)=real(sp(1:IR.API.N+1,:));
   TF.Data.Imag(ii,jj,:)=imag(sp(1:IR.API.N+1,:));
  end
end
  % calculate frequencies
TF.N=(0:fs/2/IR.API.N:fs/2)';
  % updata dimensions and save
TF=SOFAupdateDimensions(TF);
SOFAsave(fullfile(SOFAdbPath,'sofatoolbox_test','demo_FreeFieldHRTF_2_TF.sofa'),TF);

%% Plot median plane and horizontal planes for reference
figure('Name',mfilename);
SOFAplotHRTF(TF,'magmedian');
tit='SimpleFreeFieldHRTF (TF, mag), for reference'; % title
title(tit);
if savefigures==1
    saveas(gcf,[folder tit '.fig']); saveas(gcf,[folder tit '.png']); % save figures
end    
    
figure('Name',mfilename);
SOFAplotHRTF(TF,'etchorizontal'); 
tit='SimpleFreeFieldHRTF (TF, etc), for reference'; % title
title(tit);
if savefigures==1
    saveas(gcf,[folder tit '.fig']); saveas(gcf,[folder tit '.png']); % save figures
end    

%% Convert to an emitter-based representation, TFE
  % Create a new SOFA structure (based on TF because only a few metadata change)
TFE=TF; 
TFE.GLOBAL_SOFAConventions = 'GeneralTF-E';
TFE.GLOBAL_DataType = 'TF-E';
TFE.GLOBAL_History = SOFAappendText(TFE, 'GLOBAL_History', 'Converted to TFE');

  % Switch to M Emitters but only 1 Measurement
TFE.API.E=TF.API.M;
TFE.API.M=1;
  % Switch from TF to TFE
TFE.Data=rmfield(TFE.Data,{'Real','Imag'});
TFE.Data.Real(1,:,:,:)=shiftdim(TF.Data.Real,1); % MRN --> 1RNM --> MRNE with M=1
TFE.API.Dimensions.Data.Real='MRNE';
TFE.Data.Imag(1,:,:,:)=shiftdim(TF.Data.Imag,1);
TFE.API.Dimensions.Data.Imag='MRNE';
  % Update variables
TFE.EmitterPosition=TF.SourcePosition;
TFE.EmitterPosition_Type=TF.SourcePosition_Type;
TFE.EmitterPosition_Units=TF.SourcePosition_Units;
TFE.API.Dimensions.EmitterPosition='ECI';
TFE.SourcePosition=[0 0 0];
TFE.API.Dimensions.SourcePosition='IC';
  % Update dimensions and save
TFE=SOFAupdateDimensions(TFE);
SOFAsave(fullfile(SOFAdbPath,'sofatoolbox_test','demo_FreeFieldHRTF_3_TFE.sofa'),TFE);

%% Convert to SH
  % Create a new SOFA structure (based on TF because only a few metadata change)
SH=TFE;
SH.GLOBAL_SOFAConventions = 'FreeFieldHRTF';
SH.GLOBAL_History = SOFAappendText(TFE, 'GLOBAL_History', 'Converted to Spherical Harmonics');
  % Update variables
SH.EmitterPosition = [nan nan mean(TFE.EmitterPosition(:,3))];
SH.EmitterPosition_Type = 'Spherical Harmonics';
SH.EmitterPosition_Units = TFE.EmitterPosition_Units;  
  % Update dimensions
Lmax=floor(sqrt(TFE.API.E)-1); % Max SH order
L=min(33,Lmax); % don't go higher than 30th order
[S, SH.API.E]=sph2SH(TFE.EmitterPosition(:,1:2), L);
  % Allocate data
SH.Data.Real=zeros(1, SH.API.R, SH.API.N, SH.API.E);
SH.Data.Imag=zeros(1, SH.API.R, SH.API.N, SH.API.E);
  % Convert from SH to TF
Sinv=pinv(S);
for ii=1:SH.API.R
  for jj=1:SH.API.N
   SH.Data.Real(1,ii,jj,:)=Sinv*squeeze(TFE.Data.Real(1,ii,jj,:));
   SH.Data.Imag(1,ii,jj,:)=Sinv*squeeze(TFE.Data.Imag(1,ii,jj,:));
  end
end
  % update dimensions and save
SH = SOFAupdateDimensions(SH);
SOFAsave(fullfile(SOFAdbPath,'sofatoolbox_test','demo_FreeFieldHRTF_4_SH.sofa'),SH);

%% plot median and horizonal planes - spatially continuous
figure('Name',mfilename);
SOFAplotHRTF(SH,'magmedian'); %title('');
tit='FreeFieldHRTF (TFE, mag) in Spherical Harmonics'; % title
title(tit);
if savefigures==1
    saveas(gcf,[folder tit '.fig']); saveas(gcf,[folder tit '.png']); % save figures
end

figure('Name',mfilename);
SOFAplotHRTF(SH,'etchorizontal'); %title ('');
tit='FreeFieldHRTF (TFE, etc) in Spherical Harmonics'; % title
title(tit);
if savefigures==1
    saveas(gcf,[folder tit '.fig']); saveas(gcf,[folder tit '.png']); % save figures
end

%% plot spatially continuous geometry
SOFAplotGeometry(SH);
set(gcf, 'Name', mfilename);
tit='FreeFieldHRTF (TFE, geometry) in Spherical Harmonics'; % title
% if ~isoctave
if exist('OCTAVE_VERSION','builtin') == 0
  fig=gcf;
  fig.Position(3:4)=[600,400]; % increase size (supported in Matlab only)
end
title(tit);
if savefigures==1
    saveas(gcf,[folder tit '.fig']); saveas(gcf,[folder tit '.png']); % save figures
end

%% plot spatial spectra
figure('Name',mfilename);
for y=1:4
  subplot(2,2,y); hold on;
  fax=0:fs/(2*SH.API.N):(fs-(fs/(2*SH.API.N)))/2;
  r=1; plot(fax,20*log10(squeeze(abs(SH.Data.Real(1,r,:,y)+1i*SH.Data.Imag(1,r,:,y)))));
  r=2; plot(fax,20*log10(squeeze(abs(SH.Data.Real(1,r,:,y)+1i*SH.Data.Imag(1,r,:,y)))),'r');
  xlabel('Frequency (Hz)');
  ylabel('Magnitude (dB)');
  tit=['SH spectra, coefficient index (ACN)= ' num2str(y)];
  title(tit);
end
subplot(2,2,1);legend('First receiver','Second receiver');
if savefigures==1
    tit='SH spectra';
    saveas(gcf,[folder tit '.fig']); saveas(gcf,[folder tit '.png']); % save figures
end

%% plot all coefficients for a given frequency
figure('Name',mfilename); hold on;
f=10000; k=round(f/fs*2*SH.API.N);
plot(squeeze(20*log10(abs(SH.Data.Real(1,1,k,:)+1i*SH.Data.Imag(1,1,k,:)))),'ob');
hold on;
plot(squeeze(20*log10(abs(SH.Data.Real(1,2,k,:)+1i*SH.Data.Imag(1,2,k,:)))),'xr');
% title();
xlabel('Coefficients (ACN index)');
ylabel('Magnitude (dB)');
legend('First receiver','Second receiver');

tit=['SH representation, frequency = ' num2str(f) ' Hz']; % title
title(tit);
if savefigures==1
    tit='SH representations'; % title
    saveas(gcf,[folder tit '.fig']); saveas(gcf,[folder tit '.png']); % save figures
end


%% interpolate for the horizontal and median planes to FreeFieldHRTF (TFE)
TFEint=SH;
elemin=-90;
ele=[elemin:1:90 89:-1:elemin zeros(1,length(1:355))]';
azi=[zeros(length(elemin:1:90),1); 180*ones(length(89:-1:elemin),1); (1:355)'];
radius=SH.EmitterPosition(:,3)*ones(size(ele));
TFEint.EmitterPosition=[azi ele radius];
TFEint.EmitterPosition_Type='spherical';
TFEint.EmitterPosition_Units=SH.EmitterPosition_Units;
Sint = sph2SH(TFEint.EmitterPosition(:,1:2), L);
TFEint.API.E=size(Sint,1);
TFEint.Data.Real=zeros(1,2,TFEint.API.N,TFEint.API.E);
TFEint.Data.Imag=zeros(1,2,TFEint.API.N,TFEint.API.E);
for ii=1:TFEint.API.R
  for jj=1:TFEint.API.N
    TFEint.Data.Real(1,ii,jj,:)=Sint*squeeze(SH.Data.Real(1,ii,jj,:));
    TFEint.Data.Imag(1,ii,jj,:)=Sint*squeeze(SH.Data.Imag(1,ii,jj,:));
  end
end

TFEint=SOFAupdateDimensions(TFEint);

SOFAsave(fullfile(SOFAdbPath,'sofatoolbox_test','demo_FreeFieldHRTF_5_TFEint.sofa'),TFEint);

%% interpolate for the horizontal and median planes to SimpleFreeFieldHRTF (TF)
  % create a new structure and copy metadata from SH
TFint=SOFAgetConventions('SimpleFreeFieldHRTF');
TFint.GLOBAL_AuthorContact = SH.GLOBAL_AuthorContact;
TFint.GLOBAL_Comment = SH.GLOBAL_Comment;
TFint.GLOBAL_History = SOFAappendText(SH, 'GLOBAL_History', 'Converted to TF');
TFint.GLOBAL_License = SH.GLOBAL_License;
TFint.GLOBAL_Organization = SH.GLOBAL_Organization;
TFint.GLOBAL_Origin = SH.GLOBAL_Origin;
TFint.GLOBAL_Title = SH.GLOBAL_Title;
TFint.GLOBAL_DatabaseName = SH.GLOBAL_DatabaseName; 
TFint.GLOBAL_ListenerShortName = SH.GLOBAL_ListenerShortName;
TFint.GLOBAL_Author = SH.GLOBAL_Author;
TFint.GLOBAL_ListenerDescription = SH.GLOBAL_ListenerDescription;
TFint.GLOBAL_ReceiverDescription = SH.GLOBAL_ReceiverDescription;
TFint.GLOBAL_SourceDescription = SH.GLOBAL_SourceDescription;
TFint.GLOBAL_EmitterDescription = SH.GLOBAL_EmitterDescription;
TFint.GLOBAL_RoomDescription = SH.GLOBAL_RoomDescription;
TFint.GLOBAL_Parents = SH.GLOBAL_Parents;
  % Copy/Update variables
TFint.ListenerPosition=SH.ListenerPosition;
TFint.ListenerPosition_Type=SH.ListenerPosition_Type;
TFint.ListenerPosition_Units=SH.ListenerPosition_Units;
TFint.ListenerView=SH.ListenerView;
TFint.ListenerView_Type=SH.ListenerView_Type;
TFint.ListenerView_Units=SH.ListenerView_Units;
TFint.ListenerUp=SH.ListenerUp;
TFint.ReceiverPosition=SH.ReceiverPosition;
TFint.ReceiverPosition_Type=SH.ReceiverPosition_Type;
TFint.ReceiverPosition_Units=SH.ReceiverPosition_Units;
TFint.N=SH.N;
TFint.N_LongName=SH.N_LongName;
TFint.N_Units=SH.N_Units;
  % Create new Source directions
ele=[-90:0.5:90 89:-.5:-90 zeros(1,length(1:0.5:355))]';
azi=[zeros(length(-90:.5:90),1); 180*ones(length(89:-.5:-90),1); (1:0.5:355)'];
radius=1.2*ones(size(ele));
TFint.SourcePosition=[azi ele radius];
  % Update dimensions
Sint = sph2SH(TFint.SourcePosition(:,1:2), sqrt(SH.API.E)-1);
TFint.API.M=size(Sint,1);
TFint.API.R=SH.API.R; 
TFint.API.N=SH.API.N;
  % Allocate data
TFint.Data.Real=zeros(TFint.API.M,2,TFint.API.N);
TFint.Data.Imag=zeros(TFint.API.M,2,TFint.API.N);
for ii=1:TFint.API.R
  for jj=1:TFint.API.N
    TFint.Data.Real(:,ii,jj)=Sint*squeeze(SH.Data.Real(1,ii,jj,:));
    TFint.Data.Imag(:,ii,jj)=Sint*squeeze(SH.Data.Imag(1,ii,jj,:));
  end
end
  % Update dimensions and save
TFint=SOFAupdateDimensions(TFint);
SOFAsave(fullfile(SOFAdbPath,'sofatoolbox_test','demo_FreeFieldHRTF_6_TFint.sofa'),TFint);

%% compare
figure('Name',mfilename);
SOFAplotHRTF(TFint,'magmedian'); %title('');
tit='SimpleFreeFieldHRTF (TF, mag), interpolated'; % title
title(tit);
if savefigures==1
    saveas(gcf,[folder tit '.fig']); saveas(gcf,[folder tit '.png']); % save figures
end

figure('Name',mfilename);
SOFAplotHRTF(TFint,'etchorizontal'); % title('');
tit='SimpleFreeFieldHRTF (TF, etc), interpolated'; % title
title(tit);
if savefigures==1
    saveas(gcf,[folder tit '.fig']); saveas(gcf,[folder tit '.png']); % save figures
end
