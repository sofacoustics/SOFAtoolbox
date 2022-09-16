function [S, N_SH] = sph2SH(dirs, L)
%sph2SH - Calculate real-valued spherical harmonics.
%  Usage: [S, N_SH] = sph2SH(dirs, L)
% 
%  sph2SH calculates real-valued spherical harmonics for directions dirs = [azi ele] (in degrees) up to order L.
% 
%   Input:
%       dirs : dirs = [azi ele]
%       L    : order L
%
%   Output:
%       S    : real-valued spherical harmonics
%       N_SH : number of SH coefficients
% 
%   See also HOR2SPH, SPH2HOR, SPH2NAV, SPH2HOR, NAV2SPH

% #Author: Piotr Majdak: adapted from getSH.m from https://github.com/polarch/Spherical-Harmonic-Transform (24.07.2020)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)

% SOFA Toolbox - function sph2SH
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 

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
