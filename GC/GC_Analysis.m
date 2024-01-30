%% GC_Analysis
% This script visualizes the GC-MS results of the blackberry leaf
% headspace.
%
% Version:
% 05-Jan-2024 (R2023a)

% Prepare
clc; clear all; close all
warning('off')

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig)
addpath(genpath('...\GitHub\export_fig'))
mkdir('Output')

%% Settings

% Set paths
SET.main_path = '...';
SET.files = {'CHROMATOGRAM_berryleaves.txt'; 'CHROMATOGRAM_zhae.txt'};
SET.x_range = [2, 19]; %{s]

% Give annotation file
SET.annotation_file = 'berryleaves_peaks.txt';

%% Plot everything
% Iterate over chromatograms and plot the spectra
hFig = figure("Units","normalized", "Position",[0 0 1 0.5], "Color", 'w');
for iFile = 1:length(SET.files)
    % Load data
    curr_data  = readtable([SET.main_path, SET.files{iFile}]);
    % Plot data
    subplot(1,2,iFile); hold on
    plot(curr_data.Time, curr_data.Intensity/sum(curr_data.Intensity), 'color','k')
    title(SET.files{iFile}, 'Interpreter','none')
    xlim(SET.x_range)
    xlabel('time (minutes)')
    ylabel('relative abundance')
    plot([0 0],SET.x_range,'k')
end

% Add the peak annotation
curr_data  = readtable([SET.main_path, SET.files{find(contains(SET.files, 'berry'))}]);
subplot(1,2,find(contains(SET.files, 'berry'))); hold on
annotation = readtable([SET.main_path, SET.annotation_file]);
for iPeak = 1:size(annotation,1)
    % Get the corresponding intensity
    [~,idx] = min(abs(curr_data.Time-annotation.Time(iPeak)));
    val = curr_data.Intensity(idx)/sum(curr_data.Intensity);
    % Plot
    plot(annotation.Time(iPeak), val, 'ro', 'MarkerFaceColor','r', 'MarkerEdgeColor','none')
    text(annotation.Time(iPeak), val*1.01, annotation.Chemical(iPeak), "HorizontalAlignment","center")
end
export_fig('Output\GC_Chroma', '-pdf', '-painters')
close(hFig)

