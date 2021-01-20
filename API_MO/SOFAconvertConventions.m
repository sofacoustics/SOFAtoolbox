function Obj=SOFAconvertConventions(Obj)
% SOFAconvertConventions(OBJ) converts an object to SimpleFreeFieldHRIR
%
%  Supported conventions: 
%    SimpleFreeFieldSOS
%    SimpleFreeFieldTF
%    some special cases of GeneralTF.

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Convert data to SimpleFreeFieldHRIR
switch Obj.GLOBAL_SOFAConventions
  case 'SimpleFreeFieldSOS'
    N=512;
    impulse=[1; zeros(N-1,1)];
    T=SOFAgetConventions('SimpleFreeFieldHRIR');
    Obj.GLOBAL_SOFAConventions=T.GLOBAL_SOFAConventions;
    Obj.GLOBAL_SOFAConventionsVersion=T.GLOBAL_SOFAConventionsVersion;
    Obj.GLOBAL_DataType=T.GLOBAL_DataType;
    Obj.API.Dimensions.Data.IR=Obj.API.Dimensions.Data.SOS;
    Obj.Data.IR=zeros(Obj.API.M, Obj.API.R, N);
    for ii=1:Obj.API.M
      for jj=1:Obj.API.R
        Obj.Data.IR(ii,jj,:)=sosfilt(reshape(squeeze(Obj.Data.SOS(ii,jj,:)),6,[])',impulse);
      end
    end
    Obj.Data=rmfield(Obj.Data,'SOS');
    Obj.API.Dimensions.Data=rmfield(Obj.API.Dimensions.Data,'SOS');
    Obj.API.Dimensions=rmfield(Obj.API.Dimensions,'N');
    Obj=SOFAupdateDimensions(Obj);
  case {'SimpleFreeFieldTF' 'GeneralTF'}
    fs=48000; %max(Obj.N)*2;
    N=fs/min([min(diff(Obj.N)) Obj.N(1)]);
    N=2*(round(N/2+1)-1);
    T=SOFAgetConventions('SimpleFreeFieldHRIR');
    Obj.GLOBAL_SOFAConventions=T.GLOBAL_SOFAConventions;
    Obj.GLOBAL_SOFAConventionsVersion=T.GLOBAL_SOFAConventionsVersion;
    Obj.GLOBAL_DataType=T.GLOBAL_DataType;
    Obj.API.Dimensions.Data.IR=Obj.API.Dimensions.Data.Real;
    Obj.Data.SamplingRate=fs;
    Obj.Data.SamplingRate_Units='Hertz';
    Obj.Data.IR_LongName=Obj.Data.Real_LongName;
    Obj.Data.IR_Units=Obj.Data.Real_Units;
    Obj.Data.IR=zeros(Obj.API.M, Obj.API.R, N);
    for ii=1:Obj.API.M
      for jj=1:Obj.API.R
        s=zeros(N/2+1,1);
        s(Obj.N*N/fs+1)=squeeze(Obj.Data.Real(ii,jj,:))+1i*squeeze(Obj.Data.Imag(ii,jj,:));
        Obj.Data.IR(ii,jj,:)=myifftreal(s,N);
      end
      Obj.SourcePosition(ii,:)=SOFAconvertCoordinates(Obj.SourcePosition(ii,:),Obj.SourcePosition_Type,T.SourcePosition_Type,Obj.SourcePosition_Units,T.SourcePosition_Units);
    end
    Obj.Data.Delay=zeros(1,Obj.API.R);
    Obj.SourcePosition_Type=T.SourcePosition_Type;
    Obj.SourcePosition_Units=T.SourcePosition_Units;  
    Obj=rmfield(Obj,{'N','N_LongName','N_Units'});
    Obj.Data=rmfield(Obj.Data,{'Real','Imag','Real_LongName','Imag_LongName','Real_Units','Imag_Units'});
    Obj.API.Dimensions.Data=rmfield(Obj.API.Dimensions.Data,{'Real','Imag'});
    Obj.API.Dimensions=rmfield(Obj.API.Dimensions,'N');
    Obj=SOFAupdateDimensions(Obj);    
  otherwise
    error('Conventions not supported');
end


function f=myifftreal(c,N) % thanks goto the LTFAT <http://ltfat.sf.net>
if rem(N,2)==0
  f=[c; flipud(conj(c(2:end-1,:)))];
else
  f=[c; flipud(conj(c(2:end,:)))];
end;
f=real(ifft(f,N,1));
