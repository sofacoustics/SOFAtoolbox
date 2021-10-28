% Convert data from SingleFreeFieldHRIR to SingleFreeFieldHRTF, 
%   extract a few frequency bins and save.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% SOFA API - demo script
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Define parameters
% Subject index of the file to convert
subjectID='NH4';
% HRTF or DTF?
ARIfile='hrtf';
% Which frequency bins to store?
bins=[10, 20, 50, 70];
% Data compression (0..uncompressed, 9..most compressed)
compression=1; % results in a nice compression within a reasonable processing time

%% Load file in SimpleFreeFieldHRIR Conventions
f=filesep;
SOFAfn=fullfile(SOFAdbPath,'database','ari', [ARIfile '_' lower(subjectID) '.sofa']);

if isfile(SOFAfn)
    disp(['Loading:  ' SOFAfn]);
    IR=SOFAload(SOFAfn);
else
    warning(['File not existing: ' SOFAfn '  -->  Please download it from http://www.oeaw.ac.at/isf/hrtf and save it to: ' fullfile(SOFAdbPath,'database','ari', [ARIfile '_' lower(subjectID) '.sofa'])]);
    error(['Sorry.... ' mfilename ' cannot complete!']);
end

%% Get a new SimpleFreeFieldTF conventions
TF=SOFAgetConventions('SimpleFreeFieldHRTF');
disp('Converting SimpleFreeFieldHRIR to SimpleFreeFieldHRTF');

%% Copy variables and metadata
TFempty=rmfield(TF,fieldnames(SOFAgetConventions('SimpleFreeFieldHRTF','r')));  % skip all read-only metadata
Xf=fieldnames(rmfield(TFempty,{'API','Data'}));  % skip other internal
for ii=1:length(Xf)
  if isfield(IR, (Xf{ii})), TF.(Xf{ii})=IR.(Xf{ii}); end % copy if available
end

%% Transform data
TF.Data.Real=zeros(IR.API.M,IR.API.R,length(bins));
TF.Data.Imag=zeros(IR.API.M,IR.API.R,length(bins));
TF.N=(bins*IR.Data.SamplingRate/IR.API.N)';

for ii=1:IR.API.M
  cplx=fft((IR.Data.IR(ii,:,:)));
  TF.Data.Real(ii,:,:)=real(cplx(:,:,bins));
  TF.Data.Imag(ii,:,:)=imag(cplx(:,:,bins));
end

%% Update dimensions
TF=SOFAupdateDimensions(TF);

%% Save
SOFAfn=fullfile(SOFAdbPath,'sofa_api_mo_test',['ARI_' ARIfile '_' subjectID '_' num2str(length(bins)) '_freqs.sofa']);
disp(['Saving:   ' SOFAfn]);
SOFAsave(SOFAfn,IR,compression);
