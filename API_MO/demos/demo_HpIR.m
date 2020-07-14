% SOFA API - demo script
% load headphone IRs from a SOFA file from the ARI headphones database

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Define parameters
% Subject index of the file to convert
subjectID='NH5';

%% Load SOFA file
SOFAfn=fullfile(SOFAdbPath, 'headphones', 'ari', ['hpir_' lower(subjectID) '.sofa']);
disp(['Loading:  ' SOFAfn]);
X=SOFAload(SOFAfn);  

%% Plot amplitude spectra
figure;
hold on; box on;
cols=['bgrmky'];
meastime=[0; diff(X.MeasurementDate)];
for ii=1:X.API.M
  plot(20*log10(abs(fft(squeeze(X.Data.IR(ii,1,:)),X.Data.SamplingRate))),cols(ii));
  if ii>1, leg{ii}=['#' num2str(ii) ':' num2str(meastime(ii)) ' seconds later']; end
end
axis([0 X.Data.SamplingRate/2 -80 20]);
leg{1}='#1, first measurement';
legend(leg);

  