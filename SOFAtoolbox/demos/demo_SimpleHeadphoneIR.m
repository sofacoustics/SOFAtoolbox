%demo_SimpleHeadphoneIR - Load headphone IRs (SimpleHeadphoneIR convention) from a SOFA file from the ARI headphones database.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)

% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
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
figure('Name',mfilename);
hold on; box on;
cols='bgrmky';

%if ~isoctave
if exist('OCTAVE_VERSION','builtin') == 0
  if isfield(X, 'MeasurementDate')
      meastime=[0; diff(X.MeasurementDate)]; % diff not working in Octave
  else
      meastime=diff(X.GLOBAL_DateCreated); % diff not working in Octave
  end
end

for ii=1:X.API.M
  plot(20*log10(abs(fft(squeeze(X.Data.IR(ii,1,:)),X.Data.SamplingRate))),cols(ii));
  if ii>1
    % if ~isoctave; 
	if exist('OCTAVE_VERSION','builtin') == 0
      leg{ii}=['#' num2str(ii) ':' num2str(meastime(ii)) ' seconds later']; 
    else
      leg{ii}=['#' num2str(ii)]; 
    end
  end
end

for ii=1:X.API.M
  plot(20*log10(abs(fft(squeeze(X.Data.IR(ii,2,:)),X.Data.SamplingRate)))-20,cols(ii));  
end

xlim([-200 18200]);

axis([-200 18200 -65 15]);
leg{1}='#1, first measurement';
legend(leg);
title('Amplitude Spectra of Repeated Headphones Measurements (Left, Right)')
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')

  