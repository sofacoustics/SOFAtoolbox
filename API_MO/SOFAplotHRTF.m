function [M,meta,h]=SOFAplotHRTF(Obj,type,varargin)
% SOFAplotHRTF(OBJ, TYPE, CH, DIR, COLOR) plots the CH channel of HRTFs given in OBJ.
%  The following TYPEs are supported:
%  'EtcHorizontal'  energy-time curve in the horizontal plane (+/- THR)
%  'EtcMedian'      energy-time curve in the median plane (+/- THR)
%  'MagHorizontal'  magnitude spectra in the horizontal plane (+/- THR)
%  'MagMedian'      magnitude spectra in the median plane (+/- THR)
%  'magspectrum'    single magnitude spectrum for direction(s) DIR in COLOR
%  'MagSagittal'    magnitude spectra in a sagittal plane specified by OFFSET +/- THR
%
%  More options are available by SOFAplotHRTF(Obj,type,parameter,value)
%
%   Parameter
%     'ch'     receiver channel to be plotted. Default: 1
%     'dir'    fixes the positions to be plotted:
%              [azi]: shows all direction for that azimuth
%              [azi, ele]: shows all distances for that direction 
%              [azi, ele, distance]: shows only that position 
%              default: [0,0]
%     'offset' chooses a plane to be plotted. Default: 0 deg.
%     'thr'    threshold for selecting positions around a plane. Default: 2 deg.
%     'color'  color for plotting as used by PLOT
%
%
%  Supported conventions: 
%    SimpleFreeFieldHRIR
%    SimpleFreeFieldSOS
%    SimpleFreeFieldTF
%    some special cases of GeneralTF.
%
% [M,meta,h]=SOFAplotHRTF... returns the matrix M and axes (meta) displayed in the figure.
%    h is the handle of the plot.
%

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

% for backward compatibility (type as position-dependent input parameter)
if nargin == 3 && ischar(type) && isscalar(varargin{1})
%   varargin = flipud(varargin(:));
    ch = varargin{1};
    flags.do_normalize=1;
    dir=[0,0];
    color='b';
    thr=2;
    offset=0;
else
    definput.keyvals.ch=1;
    definput.keyvals.dir=[0,0];
    definput.keyvals.thr=2;
    definput.keyvals.offset=0;
    definput.flags.color={'b','r','k','y','g','c','m'};
    definput.flags.level={'normalize','absolute'};
    argin=varargin;
    for ii=1:length(argin)
        if ischar(argin{ii}), argin{ii}=lower(argin{ii}); end
    end
    [flags,kv] = SOFAarghelper({'ch','dir','thr','offset'},definput,argin);
    ch = kv.ch;
    dir = kv.dir;
    thr=kv.thr;
    color = flags.color;
    offset = kv.offset;
end
M=[];
meta=[];

%% Convert data to FIR
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
    Obj=SOFAupdateDimensions(Obj);
  case {'SimpleFreeFieldTF' 'GeneralTF'}
    fs=max(Obj.N)*2;
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
    Obj=SOFAupdateDimensions(Obj);
  case {'SimpleFreeFieldHRIR','GeneralFIR'}
  otherwise
    error('Conventions not supported');
end

fs=Obj.Data.SamplingRate;
%% Convert to spherical if cartesian
if strcmp(Obj.SourcePosition_Type,'cartesian')
    for ii=1:Obj.API.M
        [Obj.SourcePosition(ii,1),Obj.SourcePosition(ii,2),Obj.SourcePosition(ii,3)]=cart2sph(Obj.SourcePosition(ii,1),Obj.SourcePosition(ii,2),Obj.SourcePosition(ii,3));
        Obj.SourcePosition(ii,2)=rad2deg(Obj.SourcePosition(ii,2));
        Obj.SourcePosition(ii,1)=rad2deg(Obj.SourcePosition(ii,1));
        Obj.SourcePosition(ii,1)=npi2pi(Obj.SourcePosition(ii,1),'degrees');
    end
    Obj.SourcePosition_Type='spherical';
    Obj.SourcePosition_Units='degrees';
end

%% Plot according to the type
switch lower(type)
    % Energy-time curve (ETC) in the horizontal plane
  case 'etchorizontal'
    noisefloor=-50;
    ele=0;
    Obj=SOFAexpand(Obj,'Data.Delay');
    hM=double(squeeze(Obj.Data.IR(:,ch,:)));
    pos=Obj.SourcePosition;
    pos(pos(:,1)>180,1)=pos(pos(:,1)>180,1)-360;
    idx=find(pos(:,2)<(ele+thr) & pos(:,2)>(ele-thr));
    M=(20*log10(abs(hM(idx,:))));
    pos=pos(idx,:);
    del=round(Obj.Data.Delay(idx,ch));
    M2=noisefloor*ones(size(M)+[0 max(del)]);
    for ii=1:size(M,1)
      M2(ii,del(ii)+(1:Obj.API.N))=M(ii,:);
    end
    [azi,i]=sort(pos(:,1));
    M=M2(i,:);
    if flags.do_normalize
      M=M-max(max(M));
    end
    M(M<=noisefloor)=noisefloor;
    meta.time = 0:1/fs*1000:(size(M,2)-1)/fs*1000;
    meta.azi = azi;
    h=surface(meta.time,azi,M(:,:));
    set(gca,'FontName','Arial','FontSize',10);
    set(gca, 'TickLength', [0.02 0.05]);
    set(gca,'LineWidth',1);
    cmap=colormap(hot);
    cmap=flipud(cmap);
    shading flat
    colormap(cmap);
    box on;
    colorbar;
    xlabel('Time (ms)');
    ylabel('Azimuth (deg)');
    title([Obj.GLOBAL_Title '; channel: ' num2str(ch)],'Interpreter','none');

    % Magnitude spectrum in the horizontal plane
  case 'maghorizontal'
    noisefloor=-50;
    ele=0;
    hM=double(squeeze(Obj.Data.IR(:,ch,:)));
    pos=Obj.SourcePosition;
    pos(pos(:,1)>180,1)=pos(pos(:,1)>180,1)-360;
    idx=find(pos(:,2)<(ele+thr) & pos(:,2)>(ele-thr));
%     idx=find(abs(pos(:,1))>90);
%     pos(idx,2)=180-pos(idx,2);
%     pos(idx,1)=180-pos(idx,1);
%     idx=find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
    M=(20*log10(abs(fft(hM(idx,:)')')));
    M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
    pos=pos(idx,:);
    if flags.do_normalize
      M=M-max(max(M));
    end
    M(M<noisefloor)=noisefloor;
    [azi,i]=sort(pos(:,1));
    M=M(i,:);
    meta.freq = 0:fs/size(hM,2):(size(M,2)-1)*fs/size(hM,2);
    meta.azi = azi;
    h=surface(meta.freq,azi,M(:,:));
    shading flat
    xlabel('Frequency (Hz)');
    ylabel('Azimuth (deg)');
    title([Obj.GLOBAL_Title '; channel: ' num2str(ch)],'Interpreter','none');

    % Magnitude spectrum in the median plane
  case 'magmedian'
    noisefloor=-50;
    azi=0;
    hM=double(squeeze(Obj.Data.IR(:,ch,:)));
    pos=Obj.SourcePosition;
    idx=find(abs(pos(:,1))>90);
    pos(idx,2)=180-pos(idx,2);
    pos(idx,1)=180-pos(idx,1);
    idx=find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
    M=(20*log10(abs(fft(hM(idx,:)')')));
    M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
    pos=pos(idx,:);
    if flags.do_normalize
      M=M-max(max(M));
    end
    M(M<noisefloor)=noisefloor;
    [ele,i]=sort(pos(:,2));
    M=M(i,:);
    meta.freq = 0:fs/size(hM,2):(size(M,2)-1)*fs/size(hM,2);
    meta.ele = ele;
    h=surface(meta.freq,ele,M(:,:));
    shading flat
    xlabel('Frequency (Hz)');
    ylabel('Elevation (deg)');
    title([Obj.GLOBAL_Title '; channel: ' num2str(ch)],'Interpreter','none');

    % Magnitude spectrum in the median plane
  case 'magsagittal'
    noisefloor=-50;
    hM=double(squeeze(Obj.Data.IR(:,ch,:)));
    [lat,pol]=sph2hor(Obj.SourcePosition(:,1),Obj.SourcePosition(:,2));
    pos=[lat pol];
%     idx=find(abs(pos(:,1))>90);
%     pos(idx,2)=180-pos(idx,2);
%     pos(idx,1)=180-pos(idx,1);
    idx=find(pos(:,1)<(offset+thr) & pos(:,1)>(offset-thr));
    M=(20*log10(abs(fft(hM(idx,:)')')));
    M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
    pos=pos(idx,:);
    if flags.do_normalize
      M=M-max(max(M));
    end
    M(M<noisefloor)=noisefloor;
    [ele,i]=sort(pos(:,2));
    M=M(i,:);
    meta.freq = 0:fs/size(hM,2):(size(M,2)-1)*fs/size(hM,2);
    meta.ele = ele;
    h=surface(meta.freq,ele,M(:,:));
    shading flat
    xlabel('Frequency (Hz)');
    ylabel('Polar angle (deg)');
    title([Obj.GLOBAL_Title '; channel: ' num2str(ch) '; Lateral angle: ' num2str(offset) 'deg'],'Interpreter','none');
 

    % ETC in the median plane
  case 'etcmedian'
    noisefloor=-50;
    azi=0;
    Obj=SOFAexpand(Obj,'Data.Delay');
    hM=double(squeeze(Obj.Data.IR(:,ch,:)));
    pos=Obj.SourcePosition;
    idx=find(abs(pos(:,1))>90);
    pos(idx,2)=180-pos(idx,2);
    pos(idx,1)=180-pos(idx,1);
    idx=find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
    M=(20*log10(abs(hM(idx,:))));
    pos=pos(idx,:);
    del=round(Obj.Data.Delay(idx,ch));
    M2=zeros(size(M)+[0 max(del)]);
    for ii=1:size(M,1)
      M2(ii,del(ii)+(1:Obj.API.N))=M(ii,:);
    end
    if flags.do_normalize
      M=M2-max(max(M2));
    else
      M = M2;
    end
    M(M<noisefloor)=noisefloor;
    [ele,i]=sort(pos(:,2));
    M=M(i,:);
    meta.time = 0:1/fs*1000:(size(M,2)-1)/fs*1000;
    meta.ele = ele;
    h=surface(meta.time,ele,M(:,:));
    set(gca,'FontName','Arial','FontSize',10);
    set(gca, 'TickLength', [0.02 0.05]);
    set(gca,'LineWidth',1);
    cmap=colormap(hot);
    cmap=flipud(cmap);
    shading flat
    colormap(cmap);
    box on;
    colorbar;
    xlabel('Time (ms)');
    ylabel('Elevation (deg)');
    title([Obj.GLOBAL_Title '; channel: ' num2str(ch)],'Interpreter','none');

  case 'magspectrum'
    noisefloor=-50;
    pos=round(Obj.SourcePosition,1);
    switch length(dir)
        case 1
            idx=find(pos(:,1)==dir(:,1));
        case 2
            idx=find(pos(:,1)==dir(:,1) & pos(:,2)==dir(:,2));
        otherwise
            idx=find(pos(:,1)==dir(:,1) & pos(:,2)==dir(:,2) & pos(:,3)==dir(:,3));
    end
    if isempty(idx), error('Position not found'); end
    IR=squeeze(Obj.Data.IR(idx,ch,:));

    if length(idx) > 1,
        M=20*log10(abs(fft(IR')))';
        M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
        h=plot(0:fs/2/size(M,2):(size(M,2)-1)*fs/2/size(M,2),M);
        for ii=1:length(idx)
            labels{ii}=['#' num2str(idx(ii)) ': (' num2str(pos(idx(ii),1)) ', ' num2str(pos(idx(ii),2)) ')'];
        end
        legend(labels);
    else
        hM=20*log10(abs(fft(IR)));
        M=hM(1:floor(length(hM)/2));
        hold on;
        h=plot(0:fs/2/length(hM):(length(M)-1)*fs/2/length(hM),M,color,...
            'DisplayName',['#' num2str(idx) ': (' num2str(pos(idx,1)) ', ' num2str(pos(idx,2)) ')']);
        leg=legend;
        if isempty(leg),
            legend(['#' num2str(idx) ': (' num2str(pos(idx,1)) ', ' num2str(pos(idx,2)) ')']);
        else
            leg=leg.String;
            leg{end+1}=['#' num2str(idx) ': (' num2str(pos(idx,1)) ', ' num2str(pos(idx,2)) ')'];
            legend('off');
            legend(leg);
        end
    end
    ylabel('Magnitude (dB)');
    xlabel('Frequency (Hz)');
    ylim([max(max(M))+noisefloor-10 max(max(M))+10]);
    xlim([0 fs/2]);
end


function f=myifftreal(c,N) % thanks goto the LTFAT <http://ltfat.sf.net>
if rem(N,2)==0
  f=[c; flipud(conj(c(2:end-1,:)))];
else
  f=[c; flipud(conj(c(2:end,:)))];
end;
f=real(ifft(f,N,1));
