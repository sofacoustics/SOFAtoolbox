function y=rms(x,dim)

if ~exist('dim','var'), dim=1; end
y=sqrt(mean(x.*conj(x),dim));
