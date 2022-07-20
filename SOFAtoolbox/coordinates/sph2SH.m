function [S, N_SH] = sph2SH(dirs, L)
%S = sph2SH(dirs, N);
%  sph2SH calculates real-valued spherical harmonics S_l,m
%  for directions dirs = [azi ele] (in degrees) up to order L. 
%  [S, N_SH] = sph2SH(dirs, N); returns the number of SH coefficients.
% 

% #Author: Piotr Majdak: adapted from getSH.m from
%   https://github.com/polarch/Spherical-Harmonic-Transform (24.07.2020)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)


  N_dirs = size(dirs, 1); % number of discrete positions
  N_SH = (L+1)^2;         % number of SH coefficients
  dirs=deg2rad(dirs);			% convert to radiant
  dirs(:,2) = pi/2 - dirs(:,2); % convert to zenith angle

  S = zeros(N_SH, N_dirs);% matrix with the coefficients

	  % n = 0, direction-independent component
	Llm = legendre(0, cos(dirs(:,2)'));
	Nlm = sqrt(1./(4*pi)) * ones(1,N_dirs);
	CosSin = zeros(1,N_dirs);
	CosSin(1,:) = ones(1,size(dirs,1));
	S(1, :) = Nlm .* Llm .* CosSin;
	
	  % n > 0, direction-dependent components
	idx = 1;
  for l=1:L

    m = (0:l)';

    Llm = legendre(l, cos(dirs(:,2)'));
    condon = (-1).^[m(end:-1:2);m] * ones(1,N_dirs);
    Llm = condon .* [Llm(end:-1:2, :); Llm];

    mag = sqrt( (2*l+1)*factorial(l-m) ./ (4*pi*factorial(l+m)) );
      % create the ACN ordering: [m<0; m=0; m>0]
    Nlm = [mag(end:-1:2) * ones(1,N_dirs); mag * ones(1,N_dirs)]; 

    CosSin = zeros(2*l+1,N_dirs);
    CosSin(l+1,:) = ones(1,size(dirs,1)); % m=0
    CosSin(m(2:end)+l+1,:) = sqrt(2)*cos(m(2:end)*dirs(:,1)'); % m>0
    CosSin(-m(end:-1:2)+l+1,:) = sqrt(2)*sin(m(end:-1:2)*dirs(:,1)'); % m<0

    S(idx+1:idx+(2*l+1), :) = Nlm .* Llm .* CosSin;
    idx = idx + 2*l+1;
  end
    
  S = S.';
