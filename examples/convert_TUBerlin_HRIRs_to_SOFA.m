%% ===== Configuration ===================================================
HRIR_files = { ...
    'QU_KEMAR_anechoic_3m.mat'; ...
};
% Compression of data: 0 (no compression) .. 9 (maximum compression)
compression = 1;

for ii=1:length(HRIR_files)
    % Load HRIR file
    irs = read_irs(HRIR_files{ii});
    % Create output name, by replacing the file ending
    outfile = strrep(HRIR_files{ii},'mat','sofa');
    % Convert HRIR
    TUBerlin2SOFA(outfile,irs,compression);
end
