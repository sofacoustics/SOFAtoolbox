% Demonstrates the usage of SingleRoomMIMOSRIR. 
% Idea: demo_SingleRoomMIMOSRIR loads SRIRs from ???, 
% converts to SH as SingleRoomMIMOSRIR, does planewave decompositon, and stores
% the data as a ???.sofa file.

% SOFA API - demo script
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Data
% https://phaidra.kug.ac.at/view/o:104376, LS: 3 oder 5

%% Paper:
% https://www.mdpi.com/2076-3417/10/11/3747


%% Story line...
% First, the coefficients pnm up to order N are calculated using Eq. 5, 

% pnm = sum j=1 to M, a_j p(omega_j) Yn m_j.


% and then the waves directional amplitude
% density at each frequency, i.e., plane-wave decomposition, is
% computed at the desired directions using Eq. 9.
%
%  b_n (kr,ka) = 4pi i^n(j_n(kr) ? ...) Eq. 7
% 
% w (omega_l) = sum n=0 to N sum m=?n to n p_nm / b_n Y (omega_l)

% to be saved: p_nm and b_n(kr,ka) or simple a. 
% 

