function [out, yi, yj] = SOFAspat(in,Obj,azi,ele)
% SOFAspat
% [OUT, A, E] = SOFAspat(IN, OBJ, AZI, ELE) spatializes the sound IN using
% the HRTFs from OBJ according to the trajectory given in AZI and ELE.
% Input: 
%		IN: vector with the sound
%		OBJ: SOFA object containing the HRTFs
%		AZI, ELE: vectors with the trajectory (in degrees) independent for
%							azimuth and elevation
% 
% Output: 
%		OUT: binaural signal
%		A, E: azimuth and elevation of the actual trajectory (degrees)
%
% This is an example of how to use SOFA.
%
% Piotr Majdak, 2013

if ~strcmp(Obj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR')
	error('HRTFs must be saved in the SOFA conventions SimpleFreeFieldHRIR');
end

winS=1024;
N=length(in);

%% Resample the trajectory
  % create the azimuth vector
f=length(azi)-length(find(azi(find(azi==360)+1)==0)); % number of steps between azi values
iiadd=0;
if f==1 % one azi value?
    yi(1:floor(N/winS))=azi(length(azi));
else % more than one azi values
    yi=[];     f=f-1;
    for ii=1:f 
        ii=ii+iiadd;
        if azi(ii)==360
            if azi(ii+1)==0, ii=ii+1; iiadd=iiadd+1; end
            yiadd=azi(ii):(azi(ii+1)-azi(ii))/N*winS*f:azi(ii+1)-(azi(ii+1)-azi(ii))/N*winS*f ;
            yi=[yi yiadd]; % "optimal" azi values
        else
            yiadd=azi(ii):(azi(ii+1)-azi(ii))/N*winS*f:azi(ii+1)-(azi(ii+1)-azi(ii))/N*winS*f;
            yi=[yi yiadd]; % "optimal" azi values
        end
    end
end
yi=[yi azi(length(azi))]; % all "optimal" azi values

  % create the elevation vector
if length(ele)==1 % one ele value?
     yj(1:floor(N/winS))=ele(length(ele));
else % more than one ele values
    yj=[]; g=length(ele)-1;
    for ii=1:g
        yjadd=ele(ii):(ele(ii+1)-ele(ii))/N*winS*g:ele(ii+1)-(ele(ii+1)-ele(ii))/N*winS*g ;
        yj=[yj yjadd]; % "optimal" ele values
    end
end
yj=[yj ele(length(ele))]; % all "optimal" ele values

%% create a 2D-grid with nearest positions of the moving source
if length(ele)==1 && length(azi)==1
    % stationary source
  dist=(Obj.ListenerRotation(:,1)-yi(1)).^2+(Obj.ListenerRotation(:,2)-yj(1)).^2;
  [y,idx]=min(dist);
else
    % moving source
  win=min(length(yi),length(yj));    % get the number of windows
  idx=zeros(win,1);
  for ii=1:win % find nearest point on grid (LSP)
    dist=(Obj.ListenerRotation(:,1)-yi(ii)).^2+(Obj.ListenerRotation(:,2)-yj(ii)).^2;
    [y,idx(ii)]=min(dist);
  end
end

%% normalize HRTFs
hM=Obj.Data.IR;
ii=find(Obj.ListenerRotation(:,1)==0 & Obj.ListenerRotation(:,2)==0);   % search for position 0°/0°
if isempty(ii)
	peak=max([sqrt(sum(hM(:,1,:).*hM(:,1,:))) sqrt(sum(hM(:,2,:).*hM(:,2,:)))]);   % not found - normalize to IR with most energy
else
	peak=([sqrt(sum(hM(ii,1,:).*hM(ii,2,:))) sqrt(sum(hM(ii,2,:).*hM(ii,2,:)))]);  % found - normalize to this position
end
	% get the necessary HRIRs and normalize them
hM1=double(squeeze(hM(idx,1,:))/peak(1));
hM2=double(squeeze(hM(idx,2,:))/peak(2));

%% Spatialize
M=size(hM,3)-1;  % length of impulse response
    
if length(azi)==1 && length(ele)==1
		% stationary sound
	out=zeros(N,2);
	out(:,1)=fftfilt(hM1(1,:),in);
	out(:,2)=fftfilt(hM2(1,:),in);
else
		% moving source
	win=size(idx,1);
	ovlap=round(winS/2);  % add an overlap between windows
	in=[zeros(M,1); in(:,1); zeros(winS*win-N+ovlap,1)]; % zero padding of input
	out=zeros(win*winS+M+ovlap,2);     
	window = hann(winS+ovlap+M);
	for ii=1:win
		y = conv ( hM1(ii,:)' , in((ii-1)*winS+1 : ii*winS+ovlap) );
		out( (ii-1)*winS+1 : ii*winS+M+ovlap,1) = out( (ii-1)*winS+1 : ii*winS+M+ovlap,1) + y.*window; % update
		y = conv ( hM2(ii,:)' , in((ii-1)*winS+1 : ii*winS+ovlap) );        
		out( (ii-1)*winS+1 : ii*winS+M+ovlap,2) = out( (ii-1)*winS+1 : ii*winS+M+ovlap,2) + y.*window; % update
	end
end

out=out./max(max(abs(out)));