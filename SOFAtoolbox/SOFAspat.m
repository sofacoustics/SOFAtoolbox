function [out, azi, ele, idx] = SOFAspat(in,Obj,azi,ele)
%SOFAspat - Spatialize a sound source along a trajectory
%   Usage: out = SOFAspat(IN,Obj,AZI,ELE)
%          [out, azi, ele, idx] = SOFAspat(IN,Obj,AZI,ELE)
%
%   out = SOFAspat(IN, Obj, AZI, ELE) spatializes the sound IN using
%   the HRTFs from Obj along the trajectory given in AZI and ELE.
%
%   [out, azi, ele, idx] = SOFAspat(..) returns the actual trajectory 
%   and the index vector of the actually used filters from Obj. 
%   
%   Input parameters: 
%		IN:  vector with the sound
%		Obj: SOFA object containing the HRTFs
%		azi: vector with the azimuth angles (in degrees) of the trajectory
%       ele: vector with the elevation angles (in degrees) of the trajectory
%   The first and last element in azi and ele defines the beginning and end direction
%   of the spatialized source. Directions inbetween will be interpolated. 
% 
%   Output parameters: 
%		out: the spatialized binaural signal
%		azi: the azimuth angles of the actual trajectory (degrees)
%		ele: the elevation angles of the actual trajectory (degrees)
%		idx: index of the used filters corresponding to the actual trajectory
%

% #Author: Piotr Majdak (2013)
% #Author: Robert Baumgartner: adaptations (2016)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
%
% SOFA Toolbox
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
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

%% create a 2D-grid with nearest positions of the moving source
idx=zeros(S,1);
[target.x,target.y,target.z] = sph2cart(deg2rad(azi),deg2rad(ele),ones(1,S));
[pos.x,pos.y,pos.z] = sph2cart(deg2rad(Obj.SourcePosition(:,1)),...
  deg2rad(Obj.SourcePosition(:,2)),Obj.SourcePosition(:,3));
for ii=1:S % find nearest point on grid (LSP)
  dist = (pos.x-target.x(ii)).^2 + (pos.y-target.y(ii)).^2 + (pos.z-target.z(ii)).^2; 
  [~,idx(ii)]=min(dist);
end

%% normalize HRTFs to the frontal, eye-level position
% ii=find(Obj.SourcePosition(:,1)==0 & Obj.SourcePosition(:,2)==0);   % search for position 0°/0°
% if isempty(ii)
% 	peak=max([sqrt(sum(Obj.Data.IR(:,1,:).*Obj.Data.IR(:,1,:))) sqrt(sum(Obj.Data.IR(:,2,:).*Obj.Data.IR(:,2,:)))]);   % not found - normalize to IR with most energy
% else
% 	peak=([sqrt(sum(Obj.Data.IR(ii,1,:).*Obj.Data.IR(ii,1,:))) sqrt(sum(Obj.Data.IR(ii,2,:).*Obj.Data.IR(ii,2,:)))]);  % found - normalize to this position
% end

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
		ii=ii+N*hop;
		jj=jj+1;
end

%% Normalize
% out(:,1)=out(:,1)/peak(1);
% out(:,2)=out(:,2)/peak(2);

%% Output
% actually used angles
azi = Obj.SourcePosition(idx,1);
ele = Obj.SourcePosition(idx,2);
% upsampled for each sample
idup = floor(1:1/(N*hop):S+1-1/(N*hop));
azi = azi(idup);
ele = ele(idup);

end
