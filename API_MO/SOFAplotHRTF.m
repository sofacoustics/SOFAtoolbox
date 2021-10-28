function [M,meta,h]=SOFAplotHRTF(Obj,type,varargin)
% SOFAplotHRTF(OBJ, TYPE, R, DIR, COLOR) plots the R receiver of HRTFs given in OBJ.
%  The following TYPEs are supported:
%  'EtcHorizontal'  energy-time curve in the horizontal plane (+/- THR)
%  'EtcMedian'      energy-time curve in the median plane (+/- THR)
%  'MagHorizontal'  magnitude spectra in the horizontal plane (+/- THR)
%  'MagMedian'      magnitude spectra in the median plane (+/- THR)
%  'MagSpectrum'    single magnitude spectrum for direction(s) DIR in COLOR
%  'MagSagittal'    magnitude spectra in a sagittal plane specified by OFFSET +/- THR
%  'ITDhorizontal'  interaural time delay in the horizontal plane (not
%  supported in Octave)
%
%  More options are available by SOFAplotHRTF(Obj,type,parameter,value)
%
%   Parameter
%     'R'      receiver  to be plotted. Default: 1
%     'dir'    fixes the positions to be plotted:
%              [azi]: shows all direction for that azimuth
%              [azi, ele]: shows all distances for that direction 
%              [azi, ele, distance]: shows only that position 
%              default: [0,0]
%     'offset' chooses a plane to be plotted. Default: 0 deg.
%     'thr'    threshold for selecting positions around a plane. Default: 2 deg.
%     'floor'  lowest amplitude shown (dB). Default: -50 dB.
%     'convert' convert to FIR and then to TF domain. Can be set to 0 if data is available in TF domain. Default: 1.
% 
%   Additionally, 'b', 'r', 'g', etc. can be used for plotting in color 
%   as used by PLOT.
%
%
%  Supported conventions: 
%    SimpleFreeFieldHRIR
%    SimpleFreeFieldHRSOS
%    SimpleFreeFieldHRTF
%    SHFreeFieldHRTF
%    some special cases of GeneralTF, GeneralTF-E.
%
% [M,meta,h]=SOFAplotHRTF... returns the matrix M and axes (meta) displayed in the figure.
%    h is the handle of the plot.
%

% #Author: Piotr Majdak
% #Author: Michael Mihocic: type ITDhorizontal added and updated (10.2021)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
%
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

% for backward compatibility (type as position-dependent input parameter)
if nargin == 3 && ischar(type) && isscalar(varargin{1})
%   varargin = flipud(varargin(:));
    R = varargin{1};
    flags.do_normalize=1;
    dir=[0,0];
    color='b';
    thr=2;
    offset=0;
    noisefloor=-50;             
    convert=1;
else
    definput.keyvals.R=1;
    definput.keyvals.dir=[0,0];
    definput.keyvals.thr=2;
    definput.keyvals.offset=0;
    definput.keyvals.floor=-50;
    definput.flags.color={'b','r','k','y','g','c','m'};
    definput.flags.level={'normalize','absolute'};          
    definput.keyvals.convert=1;
    argin=varargin;
    for ii=1:length(argin)
        if ischar(argin{ii}), argin{ii}=lower(argin{ii}); end
    end
    [flags,kv] = SOFAarghelper({'R','dir','thr','offset','floor'},definput,argin);
    R = kv.R;
    dir = kv.dir;
    thr=kv.thr;
    color = flags.color;
    offset = kv.offset;
    noisefloor=kv.floor;         
    
    if exist('OCTAVE_VERSION','builtin')
      % We're in Octave
       if ismember(type,{'MagHorizontal','MagMedian','MagSpectrum','MagSagittal'}) && ismember(lower(Obj.GLOBAL_SOFAConventions),{'freefielddirectivitytf','generaltf'})
          % frequency domain input data only; for Octave the list has to be extended manually because 'contains' is not available
          convert=kv.convert;
      else
          convert = 1;
      end       
    else
      % We're in Matlab
      if contains(lower(type),'mag') && ismember(lower(Obj.GLOBAL_SOFAConventions),{'freefielddirectivitytf','generaltf'})
          % frequency domain input data only 
          convert=kv.convert;
      else
          convert = 1;
      end    
    end
    
end

meta=[];

if convert == 1 
    %% Convert data to FIR
    Obj=SOFAconvertConventions(Obj); 
    fs=Obj.Data.SamplingRate;
    
    %% check if receiver selection is possible
    if R > size(Obj.Data.IR,2)
        error(['Choosen receiver out of range. Only ', num2str(size(Obj.Data.IR,2)), ' receivers recorded.'])
    end
    titlepostfix='';
else
    %% check if receiver selection is possible
    if R > size(Obj.Data.Real,2)
        error(['Choosen receiver out of range. Only ', num2str(size(Obj.Data.Real,2)), ' receivers recorded.'])
    end
    titlepostfix=' (unconverted)';
end
    
%% Convert to spherical if cartesian
if strcmp(Obj.SourcePosition_Type,'cartesian')
    for ii=1:Obj.API.M
        [Obj.SourcePosition(ii,1),Obj.SourcePosition(ii,2),Obj.SourcePosition(ii,3)]=cart2sph(Obj.SourcePosition(ii,1),Obj.SourcePosition(ii,2),Obj.SourcePosition(ii,3));
        Obj.SourcePosition(ii,2)=rad2deg(Obj.SourcePosition(ii,2));
        Obj.SourcePosition(ii,1)=rad2deg(Obj.SourcePosition(ii,1));
        Obj.SourcePosition(ii,1)=npi2pi(Obj.SourcePosition(ii,1),'degrees'); % requires Mapping toolbox in Matlab
    end
    Obj.SourcePosition_Type='spherical';
    Obj.SourcePosition_Units='degrees';
end

%% Plot according to the type
switch lower(type)
    % Energy-time curve (ETC) in the horizontal plane
  case 'etchorizontal'
    Obj=SOFAexpand(Obj,'Data.Delay');
    hM=double(squeeze(Obj.Data.IR(:,R,:)));
    pos=Obj.SourcePosition;
    pos(pos(:,1)>180,1)=pos(pos(:,1)>180,1)-360;
    idx=find(pos(:,2)<(offset+thr) & pos(:,2)>(offset-thr));
    M=(20*log10(abs(hM(idx,:))));
    pos=pos(idx,:);
    del=round(Obj.Data.Delay(idx,R));
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
    title([Obj.GLOBAL_Title '; receiver: ' num2str(R)],'Interpreter','none');

    % Magnitude spectrum in the horizontal plane
  case 'maghorizontal'
        pos=Obj.SourcePosition;   % copy pos to temp. variable
        pos(pos(:,1)>180,1)=pos(pos(:,1)>180,1)-360; % find horizontal plane
        idx=find(pos(:,2)<(offset+thr) & pos(:,2)>(offset-thr)); % find indices
        pos=pos(idx,:); % truncate pos
    if convert == 1  % converted
        hM=double(squeeze(Obj.Data.IR(:,R,:)));
        M=(20*log10(abs(fft(hM(idx,:)')')));
        M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
        if flags.do_normalize
          M=M-max(max(M)); % normalize
        end

        M(M<noisefloor)=noisefloor;
        [azi,i]=sort(pos(:,1));
        M=M(i,:);
        meta.freq = 0:fs/size(hM,2):(size(M,2)-1)*fs/size(hM,2);
        meta.azi = azi;
%         figure; 


        h=surface(meta.freq,azi,M(:,:));
        
    else
      M=20*log10(abs(sqrt(squeeze(Obj.Data.Real(idx,R,:)).^2 + squeeze(Obj.Data.Imag(idx,R,:)).^2)));
        if flags.do_normalize
          M=M-max(max(M)); % normalize
        end
        M(M<noisefloor)=noisefloor;
        [azi,i]=sort(pos(:,1));
        M=M(i,:);
%         figure; 

        h=surface(Obj.N',azi,M);

    end
    shading flat
    xlabel('Frequency (Hz)');
    ylabel('Azimuth (deg)');
    title([Obj.GLOBAL_Title '; receiver: ' num2str(R) titlepostfix],'Interpreter','none');
        
    % Magnitude spectrum in the median plane
  case 'magmedian'
      azi=0;
      pos=Obj.SourcePosition;
      idx=find(abs(pos(:,1))>90);
      pos(idx,2)=180-pos(idx,2);
      pos(idx,1)=180-pos(idx,1);
      idx=find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
      pos=pos(idx,:);
      
      if convert == 1  % converted
        
        hM=double(squeeze(Obj.Data.IR(:,R,:)));
        M=(20*log10(abs(fft(hM(idx,:)')')));
        M=M(:,1:floor(size(M,2)/2));  % only positive frequencies

        if flags.do_normalize
          M=M-max(max(M));
        end
        M(M<noisefloor)=noisefloor;
        [ele,i]=sort(pos(:,2));
        M=M(i,:);
        meta.freq = 0:fs/size(hM,2):(size(M,2)-1)*fs/size(hM,2);
        meta.ele = ele;

        h=surface(meta.freq,ele,M(:,:));
      else
        M=20*log10(abs(sqrt(squeeze(Obj.Data.Real(idx,R,:)).^2 + squeeze(Obj.Data.Imag(idx,R,:)).^2)));
        if flags.do_normalize
          M=M-max(max(M)); % normalize
        end
        M(M<noisefloor)=noisefloor;
        [ele,i]=sort(pos(:,2));
        M=M(i,:);
%         figure; 
        h=surface(Obj.N',ele,M);
        
      end
    shading flat
    xlabel('Frequency (Hz)');
    ylabel('Elevation (deg)');
    title([Obj.GLOBAL_Title '; receiver: ' num2str(R) titlepostfix],'Interpreter','none');
   
    % Magnitude spectrum in the median plane
  case 'magsagittal'
      
    [lat,pol]=sph2hor(Obj.SourcePosition(:,1),Obj.SourcePosition(:,2));
    pos=[lat pol];
    idx=find(pos(:,1)<(offset+thr) & pos(:,1)>(offset-thr));
    pos=pos(idx,:);
    
    if convert == 1  % converted
    
        hM=double(squeeze(Obj.Data.IR(:,R,:)));
        M=(20*log10(abs(fft(hM(idx,:)')')));
        M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
        if flags.do_normalize
          M=M-max(max(M));
        end
        M(M<noisefloor)=noisefloor;
        [ele,i]=sort(pos(:,2));
        M=M(i,:);
        meta.freq = 0:fs/size(hM,2):(size(M,2)-1)*fs/size(hM,2);
        meta.ele = ele;
        h=surface(meta.freq,ele,M(:,:));
    else
        M=20*log10(abs(sqrt(squeeze(Obj.Data.Real(idx,R,:)).^2 + squeeze(Obj.Data.Imag(idx,R,:)).^2)));
        if flags.do_normalize
          M=M-max(max(M)); % normalize
        end
        M(M<noisefloor)=noisefloor;
        [ele,i]=sort(pos(:,2));
        M=M(i,:);
        h=surface(Obj.N',ele,M(:,:));
        
    end
    shading flat
    xlabel('Frequency (Hz)');
    ylabel('Polar angle (deg)');
    title([Obj.GLOBAL_Title '; receiver: ' num2str(R) '; Lateral angle: ' num2str(offset) 'deg' titlepostfix],'Interpreter','none');
 

    % ETC in the median plane
  case 'etcmedian'
%     noisefloor=-50;
    azi=0;
    Obj=SOFAexpand(Obj,'Data.Delay');
    hM=double(squeeze(Obj.Data.IR(:,R,:)));
    pos=Obj.SourcePosition;
    idx=find(abs(pos(:,1))>90);
    pos(idx,2)=180-pos(idx,2);
    pos(idx,1)=180-pos(idx,1);
    idx=find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
    M=(20*log10(abs(hM(idx,:))));
    pos=pos(idx,:);
    del=round(Obj.Data.Delay(idx,R));
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
    title([Obj.GLOBAL_Title '; receiver: ' num2str(R)],'Interpreter','none');

  case 'magspectrum'
    pos=round(Obj.SourcePosition*10)/10;
    switch size(dir,2)
        case 1
            aziPos = pos(:,1);
            aziDir=dir(:,1);
            aziComp = intersect(aziPos,aziDir,'rows');
            idx= find(ismember(aziPos,aziComp,'rows'));
        case 2
            aziPos = pos(:,1);
            aziDir=dir(:,1);
            elePos = pos(:,2);
            eleDir=dir(:,2);
            aziComp = intersect(aziPos,aziDir,'rows');
            eleComp = intersect(elePos,eleDir,'rows');
            idx=find(ismember(aziPos,aziComp,'rows') & ...
                ismember(elePos,eleComp,'rows'));
        otherwise
            aziPos = pos(:,1);
            aziDir=dir(:,1);
            elePos = pos(:,2);
            eleDir=dir(:,2);
            rPos = pos(:,3);
            rDir=dir(:,3);
            aziComp = intersect(aziPos,aziDir,'rows');
            eleComp = intersect(elePos,eleDir,'rows');
            rComp = intersect(rPos,rDir,'rows');
            idx=find(ismember(aziPos,aziComp,'rows') & ...
                ismember(elePos,eleComp,'rows') & ismember(rPos,rComp,'rows'));
    end
    if isempty(idx), error('Position not found'); end
    
    if convert == 1  % converted
        IR=squeeze(Obj.Data.IR(idx,R,:));
        if length(idx) > 1
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
            h=plot(0:fs/2/length(M):(length(M)-1)*fs/2/length(M),M,color,...
                'DisplayName',['#' num2str(idx) ': (' num2str(pos(idx,1)) ', ' num2str(pos(idx,2)) ')']);
            legend;
        end
        xlim([0 fs/2]);
    else
        
        M=20*log10(abs(sqrt(squeeze(Obj.Data.Real(idx,R,:)).^2 + squeeze(Obj.Data.Imag(idx,R,:)).^2)));
        
        if length(idx) > 1
            h=plot(Obj.N',M);
            for ii=1:length(idx)
                labels{ii}=['#' num2str(idx(ii)) ': (' num2str(pos(idx(ii),1)) ', ' num2str(pos(idx(ii),2)) ')'];
            end
            legend(labels);
        else
            hold on;
            h=plot(Obj.N',M,color,...
                'DisplayName',['#' num2str(idx) ': (' num2str(pos(idx,1)) ', ' num2str(pos(idx,2)) ')']);
            legend;
        end        
        
    end
    ylabel('Magnitude (dB)');
    xlabel('Frequency (Hz)');
    ylim([max(max(M))+noisefloor-10 max(max(M))+10]);
    
    % Interaural time delay in the horizontal plane
    case 'itdhorizontal'
 
      if exist('OCTAVE_VERSION','builtin')
        warning("Command 'polarplot' not supported by Octave (yet)!")
      else
        [itd, ~] = SOFAcalculateITD(Obj, 'time');
        pos = Obj.SourcePosition;
        idx=find(pos(:,2)<(offset+thr) & pos(:,2)>(offset-thr));
        itd = itd(idx);
        [pos, idx_sort] = sort(pos(idx,1));
        itd = itd(idx_sort);
        angles = deg2rad(pos);   
        %figure('Renderer', 'painters', 'Position', [10 10 700 450]); 
        polarplot(angles, itd, 'linewidth', 1.2);
        ax = gca;
        ax.ThetaDir = 'counterclockwise'; 
        ax.ThetaZeroLocation = 'top';
        rticks([max(itd)*2/3, max(itd)]); 
        rticklabels({[num2str(round(max(itd)*2/3*1e6,1)) ' ' char(181) 's'],...
                       [num2str(round(max(itd)*1e6,1)) ' ' char(181) 's']});
        thetaticks(0:30:330)
        thetaticklabels({'0°', '30°', '60°', '90°', '120°', '150°', '180°', ...
                        '210°', '240°','270°', '300°', '330°'});
        grid on;        
      end
      
  otherwise
    error([type , ' no supported plotting type.'])
end


function f=myifftreal(c,N) % thanks goto the LTFAT <http://ltfat.sf.net>
if rem(N,2)==0
  f=[c; flipud(conj(c(2:end-1,:)))];
else
  f=[c; flipud(conj(c(2:end,:)))];
end
f=real(ifft(f,N,1));
