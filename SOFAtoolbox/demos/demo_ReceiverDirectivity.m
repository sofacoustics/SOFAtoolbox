%demo_ReceiverDirectivity - Demonstrates the usage of the ReceiverDirectivity conventions.

% #Author: Piotr Majdak: adapted from FreeFieldDirectivityTF (21.11.2026)
%
% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

close all;

clear fn;
fn{1} = 'receiver/NT5_directivity_IR_01.sofa';
fn{end+1} = 'receiver/NT5_directivity_TF_01.sofa';
fn{end+1} = 'receiver/synthetic_directivity_IR_01.sofa';
fn{end+1} = 'receiver/synthetic_directivity_TF_01.sofa';

freqs = [125, 250, 500, 1000, 2000, 4000, 8000];  % Frequencies to plot (Hz)
M = 1; 
R = 1;

for ii=1:length(fn)
    disp([' Plotting (' num2str(ii) '/' num2str(length(fn)) '): ' fn{ii}]);
    try
      Obj=SOFAload(['db://database/Directivities/' fn{ii}]);
    catch
      warning('Skipping because loading failed.'); 
      continue;
    end

      % Plot the geometry because why not
    SOFAplotGeometry(Obj);
    set(gcf, 'Name', mfilename);
    title(strrep(fn{ii},'_','\_'));
    
      % If IR, convert to TF
    if strcmp(Obj.GLOBAL_SOFAConventions(end-1:end),'IR')
      hM=double(fft(Obj.Data.IR,Obj.API.N,3));
      Obj.Data.Real=real(hM);
      Obj.Data.Imag=imag(hM);
      Obj.N=0:Obj.Data.SamplingRate/Obj.API.N:Obj.Data.SamplingRate-Obj.Data.SamplingRate/Obj.API.N;
    end
      % Prepare data for a single M and R
    TF = squeeze(Obj.Data.Real(M, R, :, :) + 1i * Obj.Data.Imag(M, R, :, :))';
    freq = double(Obj.N);  % Frequency axis from the data
      % Get the positions from the EmitterPosition
    if ~strcmpi(Obj.EmitterPosition_Type,'cartesian'), warning('Skipping because EmitterPosition_Type is not cartesian'); continue; end
    pos = Obj.EmitterPosition(:,:);
    pos = SOFAconvertCoordinates(pos, 'cartesian', 'spherical');
    mask = abs(pos(:,2)) <= 2; % Filter positions for horizontal ring with an elevation tolerance
    if ~any(mask)
        min_elev = min(abs(pos(:,2)));
        mask = abs(abs(pos(:,2)) - min_elev) < 1e-3;
    end

    pos_filtered = pos(mask, :);
    TF_filtered = TF(mask, :);

    azi = mod(pos_filtered(:,1), 360);  % Azimuth in degrees (wrapped to 0–360)
    theta = deg2rad(azi);      % Convert to radians for polarplot
    [theta_sorted, idx] = sort(theta);
    figure;
    jj=0; leg={};
    for f = freqs
      [~, idxF] = min(abs(freq - f));
      if freq(idxF)>1.5*f || freq(idxF)<=0.5*f, continue; end % skip if the found frequency in `freq` would fall into other `freqs`

      mag = 20 * log10(abs(TF_filtered(:, idxF))+eps);
      if prod(isnan(mag(idx)-min(mag)))>0
          % all entries are NaN, no valid data found, skip this frequency
          continue;
      else
        if exist('OCTAVE_VERSION','builtin') == 0 
          polarplot(theta_sorted, (mag(idx)-min(mag)));
        else          
          polar(theta_sorted, (mag(idx)-min(mag)));
        end
        hold on;
        jj=jj+1;
      end

      fTitle = freq(idxF);
      if mod(fTitle,1) == 0, digitTitle = '0'; else, digitTitle = '1'; end
      if fTitle == f
        leg{jj} = sprintf(['%.' digitTitle 'f Hz'], fTitle);
      else
        leg{jj} = sprintf(['nominal %.0f Hz, actual %.' digitTitle 'f Hz'], f, fTitle);
      end
    end
    legend(leg);
    title(strrep(fn{ii},'_','\_'));
end

%%
disp('    ');
disp(['###   ' mfilename ': done   ###']);
