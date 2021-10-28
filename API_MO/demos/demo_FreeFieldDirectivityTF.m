% Demonstrates the usage of the FreeFieldDirectivityTF conventions.

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
    
    % move every figure a little to the right from previous one
    H=gcf;
    if ~isoctave; if ii>1; movegui(H,[(H.Position(1)+(ii-1)*300) H.Position(2)]); end; end
   
    % set title
    title(strrep(char(fn(ii)),'_',' '));
    disp([' Figure ' num2str(ii) ' of ' num2str(length(fn)) ' plotted: ' strrep(char(fn(ii)),'_',' ')]);
end

disp('    ');
disp('###   demo_FreeFieldDirectivityTF: done   ###');
