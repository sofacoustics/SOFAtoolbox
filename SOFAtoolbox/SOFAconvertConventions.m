function Obj=SOFAconvertConventions(Obj,varargin)
%SOFAconvertConventions - Convert from one conventions to an otherwise
%   Usage: Obj=SOFAconvertConventions(Obj)
%   
%   SOFAconvertConventions(Obj) converts a SOFA object Obj to 
%   the convention SimpleFreeFieldHRIR. The source conventions can be: 
%      SimpleFreeFieldSOS
%      SimpleFreeFieldTF
%      SimpleFreeFieldHRTF
%      FreeFieldHRTF
%      and some special cases of GeneralTF, GeneralTF-E.
%

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: rmfield replaced by SOFAremoveVariable for non data variables, improved stability (20.09.2022)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% get convention to be converted
if ~isempty(varargin)
    definput.flags.convention=SOFAgetConventions();
    definput.flags.convention=strrep(definput.flags.convention,'-','_');

    if contains(varargin,'-')
        varargin=strrep(varargin,'-','_');
    end
    [flags,~]=SOFAarghelper({},definput,varargin);

    newConvention = flags.convention;
    if contains(newConvention,'_')
        newConvention=strrep(newConvention,'_','-');
    end
    % check if given convention is available
    if isempty( SOFAgetConventions(newConvention))
       error([newConvention, ' not a valid convention.'])
    end
else
    newConvention = 'SimpleFreeFieldHRIR';
end
oldConvention = Obj.GLOBAL_SOFAConventions;
%% no conversion needed if convention is the same
if(strcmp(oldConvention,newConvention))
%     warning(['No conversion done, as Obj already conforms to ', newConvention,'.'])
    return
end

%% convert Object
newObj = SOFAgetConventions(newConvention);
oldObj = SOFAgetConventions(oldConvention);

Obj.GLOBAL_SOFAConventions=newObj.GLOBAL_SOFAConventions;
Obj.GLOBAL_SOFAConventionsVersion=newObj.GLOBAL_SOFAConventionsVersion;
Obj.GLOBAL_Conventions = newObj.GLOBAL_Conventions;
% convert Data
try
    switch Obj.GLOBAL_DataType
        case 'SOS'
            if strcmp(newObj.GLOBAL_DataType,'SOS')
                 % no data conversion needed
            elseif strcmp(newObj.GLOBAL_DataType,'FIR')
                Obj.GLOBAL_DataType=newObj.GLOBAL_DataType;
                N=512;
                impulse=[1; zeros(N-1,1)];
                Obj.API.Dimensions.Data.IR=Obj.API.Dimensions.Data.SOS;
                Obj.Data.IR=zeros(Obj.API.M, Obj.API.R, N);
                for ii=1:Obj.API.M
                  for jj=1:Obj.API.R
                    Obj.Data.IR(ii,jj,:)=sosfilt(reshape(squeeze(Obj.Data.SOS(ii,jj,:)),6,[])',impulse);
                  end
                end
                Obj.Data=rmfield(Obj.Data,'SOS');
                Obj.API.Dimensions.Data=rmfield(Obj.API.Dimensions.Data,'SOS');
            elseif strcmp(newObj.GLOBAL_DataType,'TF')
                Obj.GLOBAL_DataType=newObj.GLOBAL_DataType;
                N=512;
                impulse=[1; zeros(N-1,1)];
                Obj.API.Dimensions.Data.IR=Obj.API.Dimensions.Data.SOS;
                Data.SOS = Obj.Data.SOS;
                fs=Obj.Data.SamplingRate;
                Obj=rmfield(Obj,{'Data'});
                Obj.N=linspace(0,fs/2,Obj.API.N/2+1)';
                Obj.Data.Real = zeros(Obj.API.M,Obj.API.R,length(Obj.N));
                Obj.Data.Imag = zeros(Obj.API.M,Obj.API.R,length(Obj.N));
                Obj.N_LongName='frequency';
                Obj.N_Units='hertz';
                for ii=1:Obj.API.M
                  for jj=1:Obj.API.R
                    IR=sosfilt(reshape(squeeze(Data.SOS(ii,jj,:)),6,[])',impulse);
                    TF=fft(IR,Obj.API.N);
                    TF=TF(:,:,1:length(Obj.N));
                    Obj.Data.Real(ii,jj,:)=real(TF);
                    Obj.Data.Imag(ii,jj,:)=imag(TF);
                  end
                end
                Obj.API.Dimensions.Data=rmfield(Obj.API.Dimensions.Data,'SOS');
            elseif strcmp(newObj.GLOBAL_DataType,'TF-E')
                Obj.GLOBAL_DataType=newObj.GLOBAL_DataType;
                N=512;
                impulse=[1; zeros(N-1,1)];
                Obj.API.Dimensions.Data.IR=Obj.API.Dimensions.Data.SOS;
                Data.SOS = Obj.Data.SOS;
                fs=Obj.Data.SamplingRate;
                Obj=rmfield(Obj,{'Data'});
                Obj.N=linspace(0,fs/2,Obj.API.N/2+1)';
                Obj.Data.Real = zeros(Obj.API.M,Obj.API.R,length(Obj.N));
                Obj.Data.Imag = zeros(Obj.API.M,Obj.API.R,length(Obj.N));
                Obj.N_LongName='frequency';
                Obj.N_Units='hertz';
                for ii=1:Obj.API.M
                  for jj=1:Obj.API.R
                    IR=sosfilt(reshape(squeeze(Data.SOS(ii,jj,:)),6,[])',impulse);
                    TF=fft(IR,Obj.API.N);
                    TF=TF(:,:,1:length(Obj.N));
                    Obj.Data.Real(ii,jj,:)=real(TF);
                    Obj.Data.Imag(ii,jj,:)=imag(TF);
                  end
                end
                Obj.API.Dimensions.Data=rmfield(Obj.API.Dimensions.Data,'SOS');
                Obj=SOFAupdateDimensions(Obj);

                Obj.API.E=Obj.API.M;
                Obj.API.M=1;
                realBuffer=shiftdim(Obj.Data.Real,1);
                Obj.Data.Real=zeros(1,Obj.API.R,Obj.API.N,Obj.API.E);
                Obj.Data.Real(1,:,:,:)=realBuffer;% MRN --> 1RNM --> MRNE with M=1
                Obj.API.Dimensions.Data.Real='MRNE';
                imagBuffer=shiftdim(Obj.Data.Imag,1);
                Obj.Data.Imag=zeros(1,Obj.API.R,Obj.API.N,Obj.API.E);
                Obj.Data.Imag(1,:,:,:)=imagBuffer;
                Obj.API.Dimensions.Data.Imag='MRNE';
                Obj.EmitterPosition=Obj.SourcePosition;
                Obj.API.Dimensions.EmitterPosition='ECI';
                Obj.SourcePosition=[0 0 0];
                Obj.API.Dimensions.SourcePosition='IC';
                Data.Real=Obj.Data.Real;
                Data.Imag=Obj.Data.Imag;
                L=40; % actual SH order
                [S, Obj.API.E]=sph2SH(Obj.EmitterPosition(:,1:2), L);

                Sinv=pinv(S);
                Obj.Data.Real=zeros(1, Obj.API.R, Obj.API.N, Obj.API.E);
                Obj.Data.Imag=zeros(1, Obj.API.R, Obj.API.N, Obj.API.E);
                for ii=1:Obj.API.R
                  for jj=1:Obj.API.N
                    Obj.Data.Real(1,ii,jj,:)=Sinv*squeeze(Data.Real(1,ii,jj,:));
                  	Obj.Data.Imag(1,ii,jj,:)=Sinv*squeeze(Data.Imag(1,ii,jj,:));
                  end
                end

                Obj.EmitterPosition=mean(Obj.EmitterPosition(:,3));
                Obj.EmitterPosition_Type='Spherical Harmonics';
                Obj.EmitterPosition_Units='Metre';
            end

        case 'TF'
            if strcmp(newObj.GLOBAL_DataType,'SOS')
                % TODO
                error('Not supported yet')
            elseif strcmp(newObj.GLOBAL_DataType,'FIR')
                Obj.GLOBAL_DataType=newObj.GLOBAL_DataType;
                if sum(diff(diff(Obj.N)))
                  fs=max(Obj.N)*2;  % irregular grid, find the smallest frequency difference
                  N=fs/min([min(diff(Obj.N)) Obj.N(1)]);
                  N=2*(round(N/2+1)-1);
                  Nidx=Obj.N*N/fs+1;
                  Nsize=floor(N/2+1);
                else
                  N=2*(length(Obj.N)-1);  % regular grid (from an DFT probably), works for odd length only
                  fs=max(Obj.N)*2;
                  Nidx=1:length(Obj.N);
                  Nsize=length(Obj.N);
                end
                Obj.API.Dimensions.Data.IR=Obj.API.Dimensions.Data.Real;
                Obj.Data.SamplingRate=fs;
                Obj.Data.SamplingRate_Units='hertz';
                Obj.Data.IR=zeros(Obj.API.M, Obj.API.R, N);
                for ii=1:Obj.API.M
                  for jj=1:Obj.API.R
                    s=zeros(Nsize,1);
                    s(Nidx)=squeeze(Obj.Data.Real(ii,jj,:))+1i*squeeze(Obj.Data.Imag(ii,jj,:));
                    Obj.Data.IR(ii,jj,:)=myifftreal(s,N);
                  end
                  Obj.SourcePosition(ii,:)=SOFAconvertCoordinates(Obj.SourcePosition(ii,:),Obj.SourcePosition_Type,newObj.SourcePosition_Type,Obj.SourcePosition_Units,newObj.SourcePosition_Units);
                end
                Obj.Data.Delay=zeros(1,Obj.API.R);
                Obj.SourcePosition_Type=newObj.SourcePosition_Type;
                Obj.SourcePosition_Units=newObj.SourcePosition_Units;
                Obj=SOFAremoveVariable(Obj,'N');
%                 Obj=rmfield(Obj,{'N','N_LongName','N_Units'});
                Obj.Data=rmfield(Obj.Data,{'Real','Imag'});
                if isfield(Obj.API.Dimensions,'Data')
                    Obj.API.Dimensions.Data=rmfield(Obj.API.Dimensions.Data,{'Real','Imag'});
                end
            elseif strcmp(newObj.GLOBAL_DataType,'TF')
                 % no data conversion needed
            elseif strcmp(newObj.GLOBAL_DataType,'TF-E')
                Obj.GLOBAL_DataType = newObj.GLOBAL_DataType;
                Obj.API.E=Obj.API.M;
                Obj.API.M=1;
                realBuffer=shiftdim(Obj.Data.Real,1);
                Obj.Data.Real=zeros(1,Obj.API.R,Obj.API.N,Obj.API.E);
                Obj.Data.Real(1,:,:,:)=realBuffer;% MRN --> 1RNM --> MRNE with M=1
                Obj.API.Dimensions.Data.Real='MRNE';
                imagBuffer=shiftdim(Obj.Data.Imag,1);
                Obj.Data.Imag=zeros(1,Obj.API.R,Obj.API.N,Obj.API.E);
                Obj.Data.Imag(1,:,:,:)=imagBuffer;
                Obj.API.Dimensions.Data.Imag='MRNE';
                Obj.EmitterPosition=Obj.SourcePosition;
                Obj.API.Dimensions.EmitterPosition='ECI';
                Obj.SourcePosition=[0 0 0];
                Obj.API.Dimensions.SourcePosition='IC';
                Data.Real=Obj.Data.Real;
                Data.Imag=Obj.Data.Imag;
                L=40; % actual SH order
                [S, Obj.API.E]=sph2SH(Obj.EmitterPosition(:,1:2), L);

                Sinv=pinv(S);
                Obj.Data.Real=zeros(1, Obj.API.R, Obj.API.N, Obj.API.E);
                Obj.Data.Imag=zeros(1, Obj.API.R, Obj.API.N, Obj.API.E);
                for ii=1:Obj.API.R
                  for jj=1:Obj.API.N
                    Obj.Data.Real(1,ii,jj,:)=Sinv*squeeze(Data.Real(1,ii,jj,:));
                  	Obj.Data.Imag(1,ii,jj,:)=Sinv*squeeze(Data.Imag(1,ii,jj,:));
                  end
                end

                Obj.EmitterPosition=mean(Obj.EmitterPosition(:,3));
                Obj.EmitterPosition_Type='Spherical Harmonics';
                Obj.EmitterPosition_Units='Metre';
            end


        case 'TF-E'
            if strcmp(newObj.GLOBAL_DataType,'SOS')
                % TODO
                error('Not supported yet')
            elseif strcmp(newObj.GLOBAL_DataType,'FIR')
                Obj.GLOBAL_DataType=newObj.GLOBAL_DataType;
                Obj.API.Dimensions.Data.IR=Obj.API.Dimensions.Data.Real;
                Obj.Data.SamplingRate=max(Obj.N)*2;
                Obj.Data.SamplingRate_Units='hertz';
                % convert sperical harmonics
                if strcmpi(Obj.EmitterPosition_Type,'spherical harmonics')
                    [X,Y,Z]=sphere(60);
                    [azi_rad,ele_rad,radius]=cart2sph(X,Y,Z);
                    azi=azi_rad/pi*180;
                    ele=ele_rad/pi*180;
                    [azi,ele]=nav2sph(azi,ele);
                    radius=1.2*radius;
                    azi=azi(:);
                    ele=ele(:);
                    radius=radius(:);

                    S = sph2SH([azi ele], sqrt(Obj.API.E)-1);
                    Obj.API.M=size(S,1);
                    Obj.SourcePosition=[azi ele radius];
                    Obj.SourcePosition_Type='spherical';
                    Obj.SourcePosition_Units='degrees,degrees,metre';

                    Data.Real = zeros(Obj.API.M,2,Obj.API.N);
                    Data.Imag = zeros(Obj.API.M,2,Obj.API.N);
                    for ii=1:Obj.API.R
                      for jj=1:Obj.API.N
                        Data.Real(:,ii,jj)=S*squeeze(Obj.Data.Real(1,ii,jj,:));
                        Data.Imag(:,ii,jj)=S*squeeze(Obj.Data.Imag(1,ii,jj,:));
                      end
                    end
                else
%                     Obj.EmitterPosition=Obj.SourcePosition;
%                     Obj.EmitterPosition_Type=Obj.SourcePosition_Type;
%                     Obj.EmitterPosition_Units=Obj.SourcePosition_Units;
%                     Data.Real = shiftdim(squeeze(Obj.Data.Real(1,:,:,:)),2);
%                     Data.Imag = shiftdim(squeeze(Obj.Data.Imag(1,:,:,:)),2);
%                     Obj.API.M=Obj.API.E;
                    error(['Converting ' Obj.GLOBAL_SOFAConventions ' to ' newObj.GLOBAL_SOFAConventions ' not supported yet']);
                end

                if sum(diff(diff(Obj.N)))
                  fs=max(Obj.N)*2;  % irregular grid, find the smallest frequency difference
                  N=fs/min([min(diff(Obj.N)) Obj.N(1)]);
                  N=2*(round(N/2+1)-1);
                  Nidx=Obj.N*N/fs+1;
                  Nsize=floor(N/2+1);
                else
                  N=2*(length(Obj.N)-1);  % regular grid (from an DFT probably), works for odd length only
                  fs=max(Obj.N)*2;
                  Nidx=1:length(Obj.N);
                  Nsize=length(Obj.N);
                end

                Obj.Data.IR=zeros(Obj.API.M, Obj.API.R, N);
                for ii=1:Obj.API.M
                  for jj=1:Obj.API.R
                    s=zeros(Nsize,1);
                    s(Nidx)=squeeze(Data.Real(ii,jj,:))+1i*squeeze(Data.Imag(ii,jj,:));
                    Obj.Data.IR(ii,jj,:)=myifftreal(s,N);
                  end
                  Obj.SourcePosition(ii,:)=SOFAconvertCoordinates(Obj.SourcePosition(ii,:),Obj.SourcePosition_Type,newObj.SourcePosition_Type,Obj.SourcePosition_Units,newObj.SourcePosition_Units);
                end
                Obj.Data.Delay=zeros(Obj.API.M,Obj.API.R,1);
                Obj.Data.SamplingRate=fs;
                Obj=SOFAremoveVariable(Obj,'N');
%                 Obj=rmfield(Obj,{'N','N_LongName','N_Units'});
                Obj.Data=rmfield(Obj.Data,{'Real','Imag'});
                Obj.API.Dimensions.Data=rmfield(Obj.API.Dimensions.Data,{'Real','Imag'});
            elseif strcmp(newObj.GLOBAL_DataType,'TF')
                Obj.GLOBAL_DataType=newObj.GLOBAL_DataType;
                % convert sperical harmonics
                if strcmpi(Obj.EmitterPosition_Type,'harmonics')
                    [X,Y,Z]=sphere(60);
                    [azi_rad,ele_rad,radius]=cart2sph(X,Y,Z);
                    azi=azi_rad/pi*180;
                    ele=ele_rad/pi*180;
                    [azi,ele]=nav2sph(azi,ele);
                    radius=1.2*radius;
                    azi=azi(:);
                    ele=ele(:);
                    radius=radius(:);

                    S = sph2SH([azi ele], sqrt(Obj.API.E)-1);
                    Obj.API.M=size(S,1);
                    Obj.SourcePosition=radius;
                    Obj.SourcePosition_Type='spherical harmonics';
                    Obj.SourcePosition_Units='metre';

                    oldObj.Data.Real = Obj.Data.Real;
                    oldObj.Data.Imag = Obj.Data.Imag;
                    Obj.Data.Real = zeros(Obj.API.M,2,Obj.API.N);
                    Obj.Data.Imag = zeros(Obj.API.M,2,Obj.API.N);
                    for ii=1:Obj.API.R
                      for jj=1:Obj.API.N
                        Obj.Data.Real(:,ii,jj)=S*squeeze(oldObj.Data.Real(1,ii,jj,:));
                        Obj.Data.Imag(:,ii,jj)=S*squeeze(oldObj.Data.Imag(1,ii,jj,:));
                      end
                    end
                    Obj.EmitterPosition_Type=newObj.EmitterPosition_Type;
                    Obj.EmitterPosition_Units=newObj.EmitterPosition_Units;
                else
                    error('Conventions not supported');
                end
            elseif strcmp(newObj.GLOBAL_DataType,'TF-E')
                 % no data conversion needed
            end

        case 'FIR'
            if strcmp(newObj.GLOBAL_DataType,'SOS')
                % TODO
                error('Not supported yet')
            elseif strcmp(newObj.GLOBAL_DataType,'FIR')
                % no data conversion needed
            elseif strcmp(newObj.GLOBAL_DataType,'TF')
                Obj.GLOBAL_DataType=newObj.GLOBAL_DataType;
                IR=Obj.Data.IR;
                fs=Obj.Data.SamplingRate;
                Obj=rmfield(Obj,{'Data'});
                Obj.N=linspace(0,fs/2,Obj.API.N/2+1)';
                Obj.Data.Real = zeros(Obj.API.M,Obj.API.R,length(Obj.N));
                Obj.Data.Imag = zeros(Obj.API.M,Obj.API.R,length(Obj.N));
                Obj.N_LongName='frequency';
                Obj.N_Units='hertz';
                for ii=1:Obj.API.M
                    for jj=1:Obj.API.R
                        TF=fft(IR(ii,jj,:),Obj.API.N);
                        TF=TF(:,:,1:length(Obj.N));
                        Obj.Data.Real(ii,jj,:)=real(TF);
                        Obj.Data.Imag(ii,jj,:)=imag(TF);
                    end
                end
            elseif strcmp(newObj.GLOBAL_DataType,'TF-E')
                Obj.GLOBAL_DataType=newObj.GLOBAL_DataType;
                IR=Obj.Data.IR;
                fs=Obj.Data.SamplingRate;
                Obj=rmfield(Obj,{'Data'});
                Obj.N=linspace(0,fs/2,Obj.API.N/2+1)';
                Obj.Data.Real = zeros(Obj.API.M,Obj.API.R,length(Obj.N));
                Obj.Data.Imag = zeros(Obj.API.M,Obj.API.R,length(Obj.N));
                Obj.N_LongName='frequency';
                Obj.N_Units='hertz';
                for ii=1:Obj.API.M
                    for jj=1:Obj.API.R
                        TF=fft(IR(ii,jj,:),Obj.API.N);
                        TF=TF(:,:,1:length(Obj.N));
                        Obj.Data.Real(ii,jj,:)=real(TF);
                        Obj.Data.Imag(ii,jj,:)=imag(TF);
                    end
                end
                Obj=SOFAupdateDimensions(Obj);

                Obj.API.E=Obj.API.M;
                Obj.API.M=1;
                realBuffer=shiftdim(Obj.Data.Real,1);
                Obj.Data.Real=zeros(1,Obj.API.R,Obj.API.N,Obj.API.E);
                Obj.Data.Real(1,:,:,:)=realBuffer;% MRN --> 1RNM --> MRNE with M=1
                Obj.API.Dimensions.Data.Real='MRNE';
                imagBuffer=shiftdim(Obj.Data.Imag,1);
                Obj.Data.Imag=zeros(1,Obj.API.R,Obj.API.N,Obj.API.E);
                Obj.Data.Imag(1,:,:,:)=imagBuffer;
                Obj.API.Dimensions.Data.Imag='MRNE';
                Obj.EmitterPosition=Obj.SourcePosition;
                Obj.API.Dimensions.EmitterPosition='ECI';
                Obj.SourcePosition=[0 0 0];
                Obj.API.Dimensions.SourcePosition='IC';
                Data.Real=Obj.Data.Real;
                Data.Imag=Obj.Data.Imag;
                L=40; % actual SH order
                [S, Obj.API.E]=sph2SH(Obj.EmitterPosition(:,1:2), L);

                Sinv=pinv(S);
                Obj.Data.Real=zeros(1, Obj.API.R, Obj.API.N, Obj.API.E);
                Obj.Data.Imag=zeros(1, Obj.API.R, Obj.API.N, Obj.API.E);
                for ii=1:Obj.API.R
                  for jj=1:Obj.API.N
                    Obj.Data.Real(1,ii,jj,:)=Sinv*squeeze(Data.Real(1,ii,jj,:));
                  	Obj.Data.Imag(1,ii,jj,:)=Sinv*squeeze(Data.Imag(1,ii,jj,:));
                  end
                end

                Obj.EmitterPosition=mean(Obj.EmitterPosition(:,3));
                Obj.EmitterPosition_Type='Spherical Harmonics';
                Obj.EmitterPosition_Units='Metre';
            end
        otherwise
            error(['Turning ',oldConvention,' into ',...
                newObj.GLOBAL_SOFAConventions,' not supported.']);
    end
catch ME
    rethrow(ME);
    error(['Turning ',oldConvention,' into ',...
    newObj.GLOBAL_SOFAConventions,' not supported.']);
end
Obj.SourcePosition=SOFAconvertCoordinates(Obj.SourcePosition,Obj.SourcePosition_Type,newObj.SourcePosition_Type,Obj.SourcePosition_Units,newObj.SourcePosition_Units);
Obj.SourcePosition_Type=newObj.SourcePosition_Type;
Obj.SourcePosition_Units=newObj.SourcePosition_Units;
if strcmpi(Obj.EmitterPosition_Type,'Spherical Harmonics') && ~strcmpi(newObj.EmitterPosition_Type,'Spherical Harmonics')
    Obj.EmitterPosition=[0 0 0];
    Obj.EmitterPosition_Type=newObj.EmitterPosition_Type;
    Obj.EmitterPosition_Units=newObj.EmitterPosition_Units;
elseif ~strcmpi(newObj.EmitterPosition_Type,'Spherical Harmonics')
    Obj.EmitterPosition=SOFAconvertCoordinates(Obj.EmitterPosition,Obj.EmitterPosition_Type,newObj.EmitterPosition_Type,Obj.EmitterPosition_Units,newObj.EmitterPosition_Units);
    Obj.EmitterPosition_Type=newObj.EmitterPosition_Type;
    Obj.EmitterPosition_Units=newObj.EmitterPosition_Units;
end
Obj.ListenerPosition=SOFAconvertCoordinates(Obj.ListenerPosition,Obj.ListenerPosition_Type,newObj.ListenerPosition_Type,Obj.ListenerPosition_Units,newObj.ListenerPosition_Units);
Obj.ListenerPosition_Type=newObj.ListenerPosition_Type;
Obj.ListenerPosition_Units=newObj.ListenerPosition_Units;
Obj.ReceiverPosition=SOFAconvertCoordinates(Obj.ReceiverPosition,Obj.ReceiverPosition_Type,newObj.ReceiverPosition_Type,Obj.ReceiverPosition_Units,newObj.ReceiverPosition_Units);
Obj.ReceiverPosition_Type=newObj.ReceiverPosition_Type;
Obj.ReceiverPosition_Units=newObj.ReceiverPosition_Units;
Obj=SOFAupdateDimensions(Obj);

function f=myifftreal(c,N) % thanks goto the LTFAT <http://ltfat.sf.net>
if rem(N,2)==0
  f=[c; flipud(conj(c(2:end-1,:)))];
else
  f=[c; flipud(conj(c(2:end,:)))];
end
f=real(ifft(f,N,1));
