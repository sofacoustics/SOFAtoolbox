close all; clear; clc

% #Author: Fabian Brinkmann (12.2021)
% #Author: Michael Mihocic: adapted to SOFA demos, database (04.01.2022)

% check dependencies
if ~exist('AKp.m', 'file')
    error('This script needs the AKTools which are available from: https://www.ak.tu-berlin.de/menue/publications/open_research_tools/aktools/');
end

% SOFAstart;

% load Data
% H = SOFAload('Trumpet_modern_a4_fortissimo.sofa');
H = SOFAload('db://database/tu-berlin%20(directivity)/Trumpet_modern_a4_fortissimo.sofa');

p = H.Data.Real(:,:,1) + 1j * H.Data.Imag(:,:,1);
p_log = 20*log10(p/max(abs(p)));

% plot
AKf(20)
AKp(p_log, 'x2', 'g', H.ReceiverPosition(:, 1:2), 'dr', [-10 0], ...
    'hp_view', [30 45], 'cm', 'RdBu_flip', 'cb', 0, ...
    'sph_proc', 'interpSpline1');
title ''

AKtightenFigure