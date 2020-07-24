% Demonstrates the usage of spherical harmonics (SH) for HRTF interpolation. 
% demo_SHforHRTFs loads an HRTF set, transforms to TF, then to SH, then 
% sampled the horizontal and median plane in steps of 0.5 degrees. Finally, 
% the files are saved in SOFAdbPath as demo_SHforHRTFs_{TF, SH, TFrec}.sofa.

% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

IR=SOFAload('db://database/thk/HRIR_L2354.sofa');
fs=IR.Data.SamplingRate;
IR.GLOBAL_APIVersion=SOFAgetVersion;
figure;
SOFAplotHRTF(IR,'magmedian'); title('IR for reference');

%% Convert to TF
TF=IR;
TF.Data=rmfield(TF.Data,{'IR','Delay','SamplingRate','SamplingRate_Units'});
TF.API.Dimensions.Data=rmfield(TF.API.Dimensions.Data,{'IR','SamplingRate','Delay'});
TF.GLOBAL_SOFAConventions = 'SimpleFreeFieldHRTF';
TF.GLOBAL_DataType = 'TF';
TF.Data.Real=zeros(IR.API.M,IR.API.R,TF.API.N+1);
TF.Data.Imag=zeros(IR.API.M,IR.API.R,TF.API.N+1);
TF.Data.Real_LongName='frequency';
TF.Data.Imag_LongName='frequency';
TF.Data.Real_Units='hertz';
TF.Data.Imag_Units='hertz';
TF.N_LongName='frequency';
TF.N_Units='hertz';
for ii=1:IR.API.M
  for jj=1:IR.API.R
   sp=fft(squeeze(IR.Data.IR(ii,jj,:)),2*IR.API.N); % Delay not considered!
   TF.Data.Real(ii,jj,:)=real(sp(1:IR.API.N+1,:));
   TF.Data.Imag(ii,jj,:)=imag(sp(1:IR.API.N+1,:));
  end
end
TF.N=(0:fs/2/IR.API.N:fs/2)';

TF=SOFAupdateDimensions(TF);
figure;
SOFAplotHRTF(TF,'magmedian'); title('TF for reference');
figure;
SOFAplotHRTF(TF,'etchorizontal'); title ('TF for reference');

%% Convert to an emitter-based representation, TFE
TFE=TF; 
TFE.GLOBAL_SOFAConventions = 'GeneralTF-E';
TFE.GLOBAL_Version = '1.1';
TFE.GLOBAL_DataType = 'TF-E';
TFE.API.E=TF.API.M;
TFE.API.M=1;
TFE.Data=rmfield(TFE.Data,{'Real','Imag'});
TFE.Data.Real(1,:,:,:)=shiftdim(TF.Data.Real,1); % MRN --> 1RNM --> MRNE with M=1
TFE.API.Dimensions.Data.Real='MRNE';
TFE.Data.Imag(1,:,:,:)=shiftdim(TF.Data.Imag,1);
TFE.API.Dimensions.Data.Imag='MRNE';
TFE.EmitterPosition=TF.SourcePosition;
TFE.API.Dimensions.EmitterPosition='ECI';
TFE.SourcePosition=[0 0 0];
TFE.API.Dimensions.SourcePosition='IC';

TFE=SOFAupdateDimensions(TFE);

%% Convert to SH
SH=TFE;
SH.GLOBAL_Conventions = 'SHFreeFieldHRTF';

Lmax=floor(sqrt(size(SH.EmitterPosition,1))-1); % Max SH order
L=40; % actual SH order
[S, SH.API.E]=sph2SH(SH.EmitterPosition(:,1:2), L);

Sinv=pinv(S);
SH.Data.Real=zeros(1, SH.API.R, SH.API.N, SH.API.E);
SH.Data.Imag=zeros(1, SH.API.R, SH.API.N, SH.API.E);
for ii=1:TFE.API.R
  for jj=1:TFE.API.N
   SH.Data.Real(1,ii,jj,:)=Sinv*squeeze(TFE.Data.Real(1,ii,jj,:));
   SH.Data.Imag(1,ii,jj,:)=Sinv*squeeze(TFE.Data.Imag(1,ii,jj,:));
  end
end

SH.EmitterPosition=zeros(SH.API.E, SH.API.C);
SH.EmitterPosition_Type='Harmonics';
SH.EmitterPosition_Units='Spherical';

%% plot spatial spectra
figure;
for y=1:4
  subplot(2,2,y); hold on;
  fax=0:fs/(2*SH.API.N):(fs-(fs/(2*SH.API.N)))/2;
  r=1; plot(fax,20*log10(squeeze(abs(SH.Data.Real(1,r,:,y)+1i*SH.Data.Imag(1,r,:,y)))));
  r=2; plot(fax,20*log10(squeeze(abs(SH.Data.Real(1,r,:,y)+1i*SH.Data.Imag(1,r,:,y)))),'r');
  title(['SH spectra, coefficient index (ACN): ' num2str(y)]);
  xlabel('Frequency (Hz)');
  ylabel('Magnitude (dB)');
end

%% plot all coefficients for a given frequency
figure; hold on;
f=10000; k=round(f/fs*2*SH.API.N);
plot(squeeze(20*log10(abs(SH.Data.Real(1,1,k,:)+1i*SH.Data.Imag(1,1,k,:)))),'ob');
hold on;
plot(squeeze(20*log10(abs(SH.Data.Real(1,2,k,:)+1i*SH.Data.Imag(1,2,k,:)))),'xr');
title(['SH representation, frequency: ' num2str(f) ' Hz']);
xlabel('Coefficients (ACN index)');
ylabel('Magnitude (dB)');

%% reconstruct TF for the horizontal and median planes
TFrec=TF;
ele=[-90:0.5:90 89:-.5:-90 zeros(1,length(1:0.5:355))]';
azi=[zeros(length(-90:.5:90),1); 180*ones(length(89:-.5:-90),1); (1:0.5:355)'];
radius=1.2*ones(size(ele));
TFrec.SourcePosition=[azi ele radius];
Srec = sph2SH(TFrec.SourcePosition(:,1:2), sqrt(SH.API.E)-1);
TFrec.API.M=size(Srec,1);
TFrec.Data.Real=zeros(TFrec.API.M,2,TFrec.API.N);
TFrec.Data.Imag=zeros(TFrec.API.M,2,TFrec.API.N);
for ii=1:TFrec.API.R
  for jj=1:TFrec.API.N
    TFrec.Data.Real(:,ii,jj)=Srec*squeeze(SH.Data.Real(1,ii,jj,:));
    TFrec.Data.Imag(:,ii,jj)=Srec*squeeze(SH.Data.Imag(1,ii,jj,:));
  end
end
%% compare
figure;
SOFAplotHRTF(TFrec,'magmedian'); title('TF Reconstructed');
figure;
SOFAplotHRTF(TFrec,'etchorizontal'); title('TF Reconstructed');

%% save SOFA files to compare the sizes
SOFAsave(fullfile(SOFAdbPath,'demo_SHforHRTFs_TF.sofa'),TF);
SOFAsave(fullfile(SOFAdbPath,'demo_SHforHRTFs_SH.sofa'),SH);
SOFAsave(fullfile(SOFAdbPath,'demo_SHforHRTFs_TFrec.sofa'),TFrec);

