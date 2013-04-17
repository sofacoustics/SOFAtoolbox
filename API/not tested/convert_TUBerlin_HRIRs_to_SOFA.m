% convert_TUBerlin_HRIRs_to_SOFA
%
% This script downloads the HRIRs from the web and converts them to SOFA.
% The HRIRs measurements are available for free, can be downloaded and are
% descripted here:
% https://dev.qu.tu-berlin.de/projects/measurements/wiki/2010-11-kemar-anechoic
%
% Reference:
% Hagen Wierstorf, Matthias Geier, Alexander Raake and Sascha Spors. A Free
% Database of Head-Related Impulse Response Measurements in the Horizontal Plane
% with Multiple Distances. In 130th Convention of the Audio Engineering Society,
% May 2011.

% AUTHOR: Hagen Wierstorf


%% ===== Configuration ===================================================
% Files to convert
HRIR_files = { ...
    %'QU_KEMAR_anechoic_0.5m.mat'; ...
    %'QU_KEMAR_anechoic_1m.mat'; ...
    %'QU_KEMAR_anechoic_2m.mat'; ...
    'QU_KEMAR_anechoic_3m.mat'; ...
};
% Compression of data: 0 (no compression) .. 9 (maximum compression)
compression = 1;

for ii=1:length(HRIR_files)
    % Download HRIR files if they are not present in the current directory
    if ~exist(HRIR_files{ii},'file')
        download_cmd = sprintf('wget https://dev.qu.tu-berlin.de/projects/measurements/repository/raw/2010-11-kemar-anechoic/mat/%s', ...
            HRIR_files{ii});
        system(download_cmd);
    end
    % Load HRIR file
    irs = read_irs(HRIR_files{ii});
    % Create output name, by replacing the file ending
    outfile = strrep(HRIR_files{ii},'mat','sofa');
    % Convert HRIR
    TUBerlin2SOFA(outfile,irs,compression);
end
