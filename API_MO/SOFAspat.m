function [out, azi, ele, idx] = SOFAspat(in,Obj,azi,ele,flag)
% SOFAspat
% [out, azi, ele, idx] = SOFAspat(in,Obj,azi,ele) spatializes the sound IN using
% the HRTFs from OBJ according to the trajectory given in AZI and ELE.
% Input: 
%		in: vector with the sound
%		Obj: SOFA object containing the HRTFs
%		azi, ele: vectors with the trajectory (in degrees) independent for
%							azimuth and elevation
%   flag: 'interaural-polar': azi and ele are given in the 
%         interaural-polar coordinate system as lateral and polar angles, 
%         respectively.
% 
% Output: 
%		out: binaural signal
%		azi, ele: azimuth and elevation of the actual trajectory (degrees)
%		idx: index of the filters (corresponds to AZI and ELE)
%
% This is an example of how to use SOFA.
%
% Piotr Majdak, 2013
%

% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Define required parameters
hop=0.5;		% the hop size for the time-variant filtering (in fraction of the filter length)

%% Initial checks 
if ~strcmp(Obj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR')
	error('HRTFs must be saved in the SOFA conventions SimpleFreeFieldHRIR');
end
if min(azi)<0,	% Check for the required coordinate system
	Obj.SourcePosition(:,1)=sph2nav(Obj.SourcePosition(:,1)); % if negative azimuths are required, swith to -90/+90 system
end
N=Obj.API.N;

%% resize the input signal to be integer multiple of HRIR
L=length(in);
in=[in; zeros(N-mod(L,N),1)];
L=length(in);		% correct length of the input signal
S=L/N/hop;	% number of segments to filter

%% Resample the trajectory
if length(azi)>1, 
	azi= interp1(0:1/(length(azi)-1):1,azi,0:1/(S-1):1); 
else
	azi=repmat(azi,1,S);
end;
if length(ele)>1, 
	ele= interp1(0:1/(length(ele)-1):1,ele,0:1/(S-1):1); 
else
	ele=repmat(ele,1,S);
end;

%% Convert to spherical system if required
if exist('flag','var')
  if strcmp(flag,'horizontal-polar')
    x=SOFAconvertCoordinates([azi; ele; ones(size(azi))]','horizontal-polar','spherical');
    azi=mod(x(:,1)',360);
    ele=x(:,2)';
  end
end
%% create a 2D-grid with nearest positions of the moving source
idx=zeros(S,1);
for ii=1:S % find nearest point on grid (LSP)
    dist=(Obj.SourcePosition(:,1)-azi(ii)).^2+(Obj.SourcePosition(:,2)-ele(ii)).^2;
    [~,idx(ii)]=min(dist);
end

%% normalize HRTFs to the frontal, eye-level position
ii=find(Obj.SourcePosition(:,1)==0 & Obj.SourcePosition(:,2)==0);   % search for position 0°/0°
if isempty(ii)
	peak=max([sqrt(sum(Obj.Data.IR(:,1,:).*Obj.Data.IR(:,1,:))) sqrt(sum(Obj.Data.IR(:,2,:).*Obj.Data.IR(:,2,:)))]);   % not found - normalize to IR with most energy
else
	peak=([sqrt(sum(Obj.Data.IR(ii,1,:).*Obj.Data.IR(ii,1,:))) sqrt(sum(Obj.Data.IR(ii,2,:).*Obj.Data.IR(ii,2,:)))]);  % found - normalize to this position
end

%% Spatialize   
out=zeros(L+N/hop,2);
window=hanning(N);
ii=0;
jj=1;
iiend=L-N;
while ii<iiend    
		segT=in(ii+1:ii+N).*window;	% segment in time domain
		segF=fft(segT,2*N);	% segment in frequency domain with zero padding
		%-----------
		segFO(:,1)=squeeze(fft(Obj.Data.IR(idx(jj),1,:),2*N)).*segF;
		segFO(:,2)=squeeze(fft(Obj.Data.IR(idx(jj),2,:),2*N)).*segF;
		%-----------
		segTO=real(ifft(segFO));   % back to the time domain
		out(ii+1:ii+2*N,:)=out(ii+1:ii+2*N,:)+segTO;  % overlap and add
		ii=ii+ceil(N*hop);
		jj=jj+1;
end

%% Normalize
out=out./max(peak);

if exist('converted','var'), azi=nav2sph(azi,ele); end;