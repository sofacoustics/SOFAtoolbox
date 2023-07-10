%demo_FreeFieldDirectivityTF - Demonstrates the usage of the FreeFieldDirectivityTF conventions.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Michael Mihocic: '(demo_)plot_trumpet_directivity' added to this script (19.02.2022)
% #Author: Michael Mihocic: check for AK tools improved, link updated (10.07.2023)
% 
% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% close all; 

clear fn;
fn{1}='ITA_Dodecahedron.sofa';
fn{2}='Trumpet_modern_et_ff_all_tensorData.sofa'; 
% fn{3}='Trumpet_modern_et_ff_a4_rawData.sofa';  % (21 MB!!! uncomment if you don't mind, or have downloaded it already)


for ii=1:length(fn)
    Obj=SOFAload(['db://database/tu-berlin (directivity)/' fn{ii}]);

    %   % update dimensions to see if it works
    % D=SOFAupdateDimensions(D, 'verbose',1);

    % plot the geometry because why not
    SOFAplotGeometry(Obj);
    set(gcf, 'Name', mfilename);
    % move every figure a little to the right from previous one
    H=gcf;
    % if ~isoctave; if ii>1; movegui(H,[(H.Position(1)+(ii-1)*300) H.Position(2)]); end; end
	if exist('OCTAVE_VERSION','builtin') == 0; if ii>1; movegui(H,[(H.Position(1)+(ii-1)*300) H.Position(2)]); end; end
   
    % set title
    title(strrep(char(fn(ii)),'_',' '));
    disp([' Figure ' num2str(ii) ' of ' num2str(length(fn))+1 ' plotted: ' strrep(char(fn(ii)),'_',' ')]);
end


%%  plot_trumpet_directivity
% original code from Fabian Brinkmann (12.2021), adapted by Michael Mihocic (2021-2022)
% requires AKp.m, check dependencies, otherwise skip this part
if exist('AKp.m', 'file') && exist('AKf.m', 'file')
    
    % load Data
    H = SOFAload('db://database/tu-berlin%20(directivity)/Trumpet_modern_a4_fortissimo.sofa');
    
    p = H.Data.Real(:,:,1) + 1j * H.Data.Imag(:,:,1);
    p_log = 20*log10(p/max(abs(p)));
    
    % plot
    AKf(20)
    AKp(p_log, 'x2', 'g', H.ReceiverPosition(:, 1:2), 'dr', [-10 0], ...
        'hp_view', [30 45], 'cm', 'RdBu_flip', 'cb', 0, ...
        'sph_proc', 'interpSpline1');
    set(gcf, 'Name', mfilename);
    
    AKtightenFigure;
    disp(' Figure 3 of 3 plotted: Trumpet_modern_a4_fortissimo.sofa');
else
    warning(' Figure 3 of 3 skipped: Trumpet_modern_a4_fortissimo.sofa');
    disp('   To plot this figure you need to:');
    disp('   - download AKtools.zip (release) from: https://www.tu.berlin/ak/forschung/publikationen/open-research-tools/aktools');
    disp('   - extract the zip file');
    disp('   - run AKtoolsStart.m');
    disp('   - add the path of AKtoolsStart.m to Matlab (including subdirectories!)');
end

%%
disp('    ');
disp('###   demo_FreeFieldDirectivityTF: done   ###');
