function M=SOFAplotHRTF(Obj,type,ch)
% SOFAplotHRTF(OBJ, TYPE, CH) plots the CH channel of HRTFs given in OBJ as
% SimpleFreeFieldHRIR. The following TYPEs are supported:
%  'EtcHorizontal'  energy-time curve in the horizontal plane (+/- 5 deg)
%  'EtcMedian'      energy-time curve in the median plane (+/- 2 deg)
%  'MagMedian'      magnitude spectrum in the median plane (+/- 2 deg)
%
% M=SOFAplotHRTF... returns the matrix M displayed in the figure.
%

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


if ~exist('ch','var')
    ch=1;
end
fs=Obj.Data.SamplingRate;

switch lower(type)
  case 'etchorizontal'
    noisefloor=-50;
    ele=0;
    thr=5;

    hM=double(squeeze(Obj.Data.IR(:,ch,:)));
    pos=Obj.SourcePosition;
    pos(pos(:,1)>180,1)=pos(pos(:,1)>180,1)-360;
    idx=find(pos(:,2)<(ele+thr) & pos(:,2)>(ele-thr));
    M=(20*log10(abs(hM(idx,:))));
    pos=pos(idx,:);
    M=M-max(max(M));
    M(M<noisefloor)=noisefloor;
    [azi,i]=sort(pos(:,1));
    M=M(i,:);
    surface(0:1/fs*1000:(size(M,2)-1)/fs*1000,azi,M(:,:));
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
    
  case 'magmedian'
    noisefloor=-50;
    azi=0;
    thr=2;
    hM=double(squeeze(Obj.Data.IR(:,ch,:)));
    pos=Obj.SourcePosition;
    idx=find(abs(pos(:,1))>90);
    pos(idx,2)=180-pos(idx,2);
    pos(idx,1)=180-pos(idx,1);    
    idx=find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
    M=(20*log10(abs(fft(hM(idx,:)')')));
    M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
    pos=pos(idx,:);
    M=M-max(max(M));
    M(M<noisefloor)=noisefloor;
    [ele,i]=sort(pos(:,2));
    M=M(i,:);
    surface(0:fs/size(hM,2):(size(M,2)-1)*fs/size(hM,2),ele,M(:,:));
    shading flat
    xlabel('Frequency (Hz)');
    ylabel('Elevation (deg)');
    title([Obj.GLOBAL_Title '; channel: ' num2str(ch)],'Interpreter','none');
  case 'etcmedian'
    noisefloor=-50;
    azi=0;
    thr=2;
    hM=double(squeeze(Obj.Data.IR(:,ch,:)));
    pos=Obj.SourcePosition;
    idx=find(abs(pos(:,1))>90);
    pos(idx,2)=180-pos(idx,2);
    pos(idx,1)=180-pos(idx,1);    
    idx=find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
    M=(20*log10(abs(hM(idx,:))));
    pos=pos(idx,:);
    M=M-max(max(M));
    M(M<noisefloor)=noisefloor;
    [ele,i]=sort(pos(:,2));
    M=M(i,:);
    surface(0:1/fs*1000:(size(M,2)-1)/fs*1000,ele,M(:,:));
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
end