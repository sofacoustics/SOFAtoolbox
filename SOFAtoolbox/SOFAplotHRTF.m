function [M,meta,h]=SOFAplotHRTF(Obj,type,varargin)
%SOFAplotHRTF - Plot HRTFs 
%   Usage: SOFAplotHRTF(Obj, type)
%          SOFAplotHRTF(Obj, type, R) 
%          SOFAplotHRTF(Obj, type, key,value)
%
%   SOFAplotHRTF(Obj, type) plots HRTFs of the first receiver (left ear) 
%   provided in Obj as the figure described by type: 
%     'EtcHorizontal'  Energy-Time Curve (ETC) in the horizontal plane (+/- THR)
%     'EtcMedian'      ETC in the median plane (+/- THR)
%     'MagHorizontal'  Magnitude spectra in the horizontal plane (+/- THR)
%     'MagMedian'      Magnitude spectra in the median plane (+/- THR)
%     'MagSagittal'    Magnitude spectra in a sagittal plane specified by OFFSET +/- THR
%     'MagSpectrum'    A single magnitude spectrum for given direction(s) DIR in COLOR
%     'ITDhorizontal'  Interaural time delays (ITDs) in the horizontal plane (Matlab only)
%
%   SOFAplotHRTF(Obj,type,key,value) defines the plotting in more detail:
%     'receiver' : Plot only specific receiver. Default: 1
%     'dir'      : Plot only specific positions (specify up to 0.1 accuracy).                   
%                  [azi]:         show all directions described by the azimuth angles azi.
%                  [azi, ele]:    show all directions described by the pairs of 
%                                 azimuth and elevation angles azi an ele, respectively.
%                  [azi, ele, r]: show only that specific positions described by the 
%                                 spherical coordinates azi, ele, and r.
%                  Default: [0,0].
%     'offset'   : Shifts the plane for which the HRTFs will be plotted. The amount of the 
%                  shift is described by value in degrees. 
%                  For plots in the horizontal plane, offset shifts the plane to be above the eye level.
%                  For plots in the the sagittal plane, offset shifts the plane to the left. 
%                  Default: 0 deg.
%     'thr'      : Threshold for considering the data in the proximity of the selected plane. 
%                  Default: 2 deg.
%     'floor'    : Lowest amplitude (dB) shown in the plot. Default: -50 dB.
%
%   SOFAplotHRTF(Obj,type,..,flags) defines the following flags:
%     'b','r','k','y','g','c','m'       : Color of the plots (MagSpectrum only)
%     'normalization','nonormalization' : Normalize the data before applying the floor. Default: normalization.
%     'conversion2ir','noconversion2ir' : Convert to the IR domain before plotting. SOFAplotHRTF 
%                                         automatically selects if the convertion is required but this 
%                                         mechanism can be overwritten with this flag. 
%      'Threshold' (default),'Cen_e2','MaxIACCr', 'MaxIACCe','CenIACCr','CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD': 
%                                         To select an ITD estimation method 
%                                         (see SOFAcalculateITD for more details)
%                                         
%   SOFAplotHRTF supports the following conventions: 
%     SimpleFreeFieldHRIR
%     SimpleFreeFieldHRSOS
%     SimpleFreeFieldHRTF
%     SHFreeFieldHRTF
%     and some special cases of GeneralTF, GeneralTF-E.
%
%   [M,meta,h]=SOFAplotHRTF(..) returns the matrix M, the meta information about 
%   the displayed the figure, and the handle h of the plot. meta contains the field idx
%   which is the index to the vectors actually plotted.
%

% #Author: Piotr Majdak
% #Author: Michael Mihocic: type ITDhorizontal added and updated (10.2021)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: dependency on function 'npi2pi' removed (required toolbox in Matlab; in Octave not supported; outdated anyway) (08.02.2022)
% #Author: Michael Mihocic: keyvalue 'R'/'r' renamed to 'receiver' (10.05.2022)
% #Author: Michael Mihocic: 'do not convert' option enabled for SimpleFreeFieldHRTF convention; 
%                           global title only displayed if 'Obj.GLOBAL_Title' not empty (30.05.2022)
% #Author: Michael Mihocic: plotting simplefreefieldhrtf fixed (in Octave); titles fixed when plotting magspectrum (02.06.2022) 
% #Author: Michael Mihocic: plotting improved when data is available in TF format (more stable, no conversions by default);
%                           figure titles improved (04.07.2022) 
% #Author: Piotr Majdak: conversion to TF is a flag now. It's called convert2TF.
% #Author: Michael Mihocic: flag convert2TF renamed to conversion2ir and noconversion2ir. (01.09.2022)
% #Author: Robert Baumgartner: changed ITDhorizontal to plot absolute ITD values. (24.11.2022)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

% for backward compatibility (type as position-dependent input parameter)
if nargin == 3 && ischar(type) && isscalar(varargin{1})
%   varargin = flipud(varargin(:));
    R = varargin{1};
    flags.do_normalization=1;
    dir=[0,0];
    color='b';
    thr=2;
    offset=0;
    noisefloor=-50;             
%     convert=1; more comples differing below:

    if exist('OCTAVE_VERSION','builtin')
      % We're in Octave
       if ismember(type,{'MagHorizontal','MagMedian','MagSpectrum','MagSagittal'}) && ismember(lower(Obj.GLOBAL_SOFAConventions),{'freefielddirectivitytf','generaltf','simplefreefieldhrtf'})
          % In Octave 'contains' is not available, thus, the list has to be extended manually 
          do_conversion2ir = 0;
      else
          do_conversion2ir = 1;
      end       
    else
      % We're in Matlab
      if contains(lower(type),'mag') && ismember(lower(Obj.GLOBAL_SOFAConventions),{'freefielddirectivitytf','generaltf','simplefreefieldhrtf'})
          % frequency domain input data only 
          do_conversion2ir = 0;
      else
          do_conversion2ir = 1;
      end    
    end

else
    definput.keyvals.receiver=1;
    definput.keyvals.dir=[0,0];
    definput.keyvals.thr=2;
    definput.keyvals.offset=0;
    definput.keyvals.floor=-50;
    definput.flags.color={'b','r','k','y','g','c','m'};
    definput.flags.normalization={'normalization','nonormalization'};          
    definput.flags.conversion2ir={'conversion2ir','noconversion2ir'};
    definput.flags.itdestimator = {'Threshold','Cen_e2','MaxIACCr', 'MaxIACCe', 'CenIACCr', 'CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD'};
    argin=varargin;
%     for ii=1:length(argin)
%         if ischar(argin{ii}), argin{ii}=lower(argin{ii}); end
%     end
    [flags,kv] = SOFAarghelper({'receiver','dir','thr','offset','floor'},definput,argin);
    R = kv.receiver;
    dir = kv.dir;
    thr=kv.thr;
    color = flags.color;
    offset = kv.offset;
    noisefloor=kv.floor;      
    do_conversion2ir=flags.do_conversion2ir; % force convertion to TF (or not)
    
end

meta=[];

if do_conversion2ir == 1 
    %% Convert data to FIR
    Obj=SOFAconvertConventions(Obj);   
    fs=Obj.Data.SamplingRate;
    
    %% check if receiver selection is possible
    if R > size(Obj.Data.IR,2)
        error(['Choosen receiver out of range. Only ', num2str(size(Obj.Data.IR,2)), ' receivers recorded.'])
    end
    titlepostfix=' (converted to IR)';
else
    %% check if receiver selection is possible
    if R > size(Obj.Data.Real,2)
        error(['Choosen receiver out of range. Only ', num2str(size(Obj.Data.Real,2)), ' receivers recorded.'])
    end
    titlepostfix='';
end
if isfield(Obj, 'GLOBAL_Title') && isempty(Obj.GLOBAL_Title) == 0
    titleprefix = [Obj.GLOBAL_Title ': '];
else
    titleprefix = '';
end

    
%% Convert to spherical if cartesian
if strcmp(Obj.SourcePosition_Type,'cartesian')
% %     Obj2=Obj; % compare to old method (Obj2)
    for ii=1:Obj.API.M
        [Obj.SourcePosition(ii,1),Obj.SourcePosition(ii,2),Obj.SourcePosition(ii,3)]=cart2sph(Obj.SourcePosition(ii,1),Obj.SourcePosition(ii,2),Obj.SourcePosition(ii,3));
        Obj.SourcePosition(ii,2)=rad2deg(Obj.SourcePosition(ii,2));
        Obj.SourcePosition(ii,1)=rad2deg(Obj.SourcePosition(ii,1));
        Obj.SourcePosition(ii,1)=mywrapTo180(Obj.SourcePosition(ii,1));
    end
    Obj.SourcePosition_Type='spherical';
    Obj.SourcePosition_Units='degrees,degrees,metre';
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
    meta.idx=idx;
    M2=noisefloor*ones(size(M)+[0 max(del)]);
    for ii=1:size(M,1)
      M2(ii,del(ii)+(1:Obj.API.N))=M(ii,:);
    end
    [azi,i]=sort(pos(:,1));
    M=M2(i,:);
    if flags.do_normalization
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
    title([titleprefix 'receiver: ' num2str(R)],'Interpreter','none');

    % Magnitude spectrum in the horizontal plane
  case 'maghorizontal'
        pos=Obj.SourcePosition;   % copy pos to temp. variable
        pos(pos(:,1)>180,1)=pos(pos(:,1)>180,1)-360; % find horizontal plane
        idx=find(pos(:,2)<(offset+thr) & pos(:,2)>(offset-thr)); % find indices
        pos=pos(idx,:); % truncate pos
        meta.idx=idx;
    if do_conversion2ir == 1  % converted
        hM=double(squeeze(Obj.Data.IR(:,R,:)));
        M=(20*log10(abs(fft(hM(idx,:)')')));
        M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
        if flags.do_normalization
          M=M-max(max(M));
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
        if flags.do_normalization
          M=M-max(max(M)); 
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
    title([titleprefix 'receiver: ' num2str(R) titlepostfix],'Interpreter','none');
        
    % Magnitude spectrum in the median plane
  case 'magmedian'
      azi=0;
      pos=Obj.SourcePosition;
      idx0=find(abs(pos(:,1))>90);
      pos(idx0,2)=180-pos(idx0,2);
      pos(idx0,1)=180-pos(idx0,1);
      idx=find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
      pos=pos(idx,:);
      meta.idx=idx; % PM: TODO: check if the correct index
      
      if do_conversion2ir == 1  % converted
        
        hM=double(squeeze(Obj.Data.IR(:,R,:)));
        M=(20*log10(abs(fft(hM(idx,:)')')));
        M=M(:,1:floor(size(M,2)/2));  % only positive frequencies

        if flags.do_normalization
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
        if flags.do_normalization
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
    title([titleprefix 'receiver: ' num2str(R) titlepostfix],'Interpreter','none');
   
    % Magnitude spectrum in the median plane
  case 'magsagittal'
      
    [lat,pol]=sph2hor(Obj.SourcePosition(:,1),Obj.SourcePosition(:,2));
    pos=[lat pol];
    idx=find(pos(:,1)<(offset+thr) & pos(:,1)>(offset-thr));
    pos=pos(idx,:);
    meta.idx=idx;
    
    if do_conversion2ir == 1  % converted
    
        hM=double(squeeze(Obj.Data.IR(:,R,:)));
        M=(20*log10(abs(fft(hM(idx,:)')')));
        M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
        if flags.do_normalization
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
        if flags.do_normalization
          M=M-max(max(M));
        end
        M(M<noisefloor)=noisefloor;
        [ele,i]=sort(pos(:,2));
        M=M(i,:);
        h=surface(Obj.N',ele,M(:,:));
        
    end
    shading flat
    xlabel('Frequency (Hz)');
    ylabel('Polar angle (deg)');
    title([titleprefix 'receiver: ' num2str(R) '; Lateral angle: ' num2str(offset) 'deg' titlepostfix],'Interpreter','none');
 

    % ETC in the median plane
  case 'etcmedian'
%     noisefloor=-50;
    azi=0;
    Obj=SOFAexpand(Obj,'Data.Delay');
    hM=double(squeeze(Obj.Data.IR(:,R,:)));
    pos=Obj.SourcePosition;
    idx0=find(abs(pos(:,1))>90);
    pos(idx0,2)=180-pos(idx0,2);
    pos(idx0,1)=180-pos(idx0,1);
    idx=find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
    meta.idx=idx; % PM: TODO: Check if the correct index
    M=(20*log10(abs(hM(idx,:))));
    pos=pos(idx,:);
    del=round(Obj.Data.Delay(idx,R));
    M2=zeros(size(M)+[0 max(del)]);
    for ii=1:size(M,1)
      M2(ii,del(ii)+(1:Obj.API.N))=M(ii,:);
    end
    if flags.do_normalization
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
    title([titleprefix 'receiver: ' num2str(R)],'Interpreter','none');

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
    meta.idx=idx;
    
    if do_conversion2ir == 1  % convert
        IR=squeeze(Obj.Data.IR(idx,R,:));
        if length(idx) > 1
            M=20*log10(abs(fft(IR')))';
            M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
            h=plot(0:fs/2/size(M,2):(size(M,2)-1)*fs/2/size(M,2),M);
            for ii=1:length(idx)
                labels{ii}=['#' num2str(idx(ii)) ': (' num2str(pos(idx(ii),1)) ', ' num2str(pos(idx(ii),2)) ')'];
            end
            legend(labels);
        else % only one curve
            hM=20*log10(abs(fft(IR)));
            M=hM(1:floor(length(hM)/2));
            hold on;
            h=plot(0:fs/2/length(M):(length(M)-1)*fs/2/length(M),M,color,...
                'DisplayName',['#' num2str(idx) ': (' num2str(pos(idx,1)) ', ' num2str(pos(idx,2)) ')']);
            legend;
        end
        xlim([0 fs/2]);
        titlepostfix=' (converted to IR)';
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
        titlepostfix='';
    end
    ylabel('Magnitude (dB)');
    xlabel('Frequency (Hz)');
    ylim([max(max(M))+noisefloor-10 max(max(M))+10]);
    title([titleprefix 'receiver: ' num2str(R) titlepostfix],'Interpreter','none');

    % Interaural time delay in the horizontal plane
    case 'itdhorizontal'
 
      if exist('OCTAVE_VERSION','builtin')
        warning('Command ''polarplot'' not supported by Octave (yet)!')
      else
        [itd, ~] = SOFAcalculateITD(Obj, 'time',flags.itdestimator);
        pos = Obj.SourcePosition;
        idx=find(pos(:,2)<(offset+thr) & pos(:,2)>(offset-thr));
        itd = itd(idx);
        meta.idx=idx;
        [pos, idx_sort] = sort(pos(idx,1));
        itd = itd(idx_sort);
        angles = deg2rad(pos);   
        %figure('Renderer', 'painters', 'Position', [10 10 700 450]); 
        polarplot(angles, abs(itd), 'linewidth', 1.2);
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


% function f=myifftreal(c,N) % thanks goto the LTFAT <http://ltfat.sf.net>
%     if rem(N,2)==0
%       f=[c; flipud(conj(c(2:end-1,:)))];
%     else
%       f=[c; flipud(conj(c(2:end,:)))];
%     end
%     f=real(ifft(f,N,1));
% end

function newangle = mywrapTo180(angle)
    % transfer to range -180:180
    newangle = mod(angle+360, 360);
    if newangle > 180
        newangle = newangle-360;
    end
