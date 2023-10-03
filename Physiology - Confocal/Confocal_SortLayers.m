%% ConfocalAnalysis_OdorSpace
% This script pools all imaging data and sorts each ROIs activity in two
% ways. First, we use k-means clustering to determine different response
% profiles, i.e., stereotyped time courses ROIs show in response to the
% different stimuli. For this, we vertically concatenate all ROIs from all
% animals and all their responses to all stimuli.
% Second, we use k-means clustering again to determine different response
% vectros, i.e., stereotyped responses to a set of stimuli (e.g. active
% with 2 out of 4; or specificallz active to one). For this, we concatenate
% all responses (time courses and average values during the active window)
% of a ROI horizontally (i.e., one row per ROI) and all ROIs of all animals
% vertically.
% At the end of the script, we plot the resulting PCA 'olfaction' spaces,
% response profiles, response vectors and other diagrams describing the
% combinatorial coding of locust antennal lobes.
%
% Version:
% 29-June-2023 (R2023a) Yannick Günzel

% Prepare
clc; clear all; close all
warning('off')

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig)
addpath(genpath('...\GitHub\export_fig'))


%% Settings

% Set paths
SET.main_path = '...\Data\Physiology\Cal520_Confocal\';

% Give list of phases
SET.phases = {'gregarious'; 'solitarious'};

% Give layer names
SET.layers = {'layer01_top', 'layer02_middle', 'layer03_bottom'};

% Give list of stimuli (We have this in the meta data, but good to double
% check)
SET.StimList = {...
    'ZHAE';...
    'COL';...
    'ZHAECOL'};

% Clusterng & Dimensionality reduction
SET.kmeans.MaxIter = 1e4;
SET.kmeans.OnlinePhase = 'on';
SET.kmeans.Replicates = 5e1;
SET.maxk_eval = 15;
SET.maxK_threshold = 1/exp(1);
SET.kmeans.dist_RespProfile = 'correlation';
SET.kmeans.dist_AvgValue = 'sqeuclidean';
SET.dist_pca = 'correlation';
SET.pca_explained = 99;


%% Pool data

% Iterate over both phases
for iPhase = 1:length(SET.phases)

    % Prepare
    all_data.(SET.phases{iPhase}) = [];
    all_info.(SET.phases{iPhase}) = {};
    cnt = 1;

    % And iterate over all imaging layers
    for iLayer = 1:length(SET.layers)

        % Make a new folder for each layer
        mkdir(['ConfocalAnalysis_OdorSpace\', SET.layers{iLayer}])

        % Contruct the path
        basepath = [SET.main_path, SET.phases{iPhase},'\', SET.layers{iLayer},'\'];
        all_animals = dir(basepath);

        figure('units', 'normalized', 'Position',[0 0 1 1], 'name', SET.layers{iLayer})
        % Iterate over all animals
        for iAni = 1:size(all_animals)

            % Check whether the current folder is belonging to an animal
            if contains(all_animals(iAni).name, '_Animal')

                % Get the correct path
                path2data = [basepath, all_animals(iAni).name, '\03_Data_Processed\'];

                % Load meta and segmentation
                load([path2data, 'meta_info.mat'])

                % Get data
                stack = [];
                for iStim = 1:length(SET.StimList)
                    % Load data
                    curr_data = load([path2data, SET.StimList{iStim}, '.mat']);
                    curr_data = curr_data.ImageStream.(SET.StimList{iStim});
                    stack = cat(3,stack, curr_data);
                end

                % Get std-projection
                stack = nanstd(stack,[],3);

                % Reshape to row vector
                stack = reshape(stack, [1, numel(stack)]);

                % Normalize
                stack = stack - min(stack);
                stack = sqrt(stack);
                stack = stack - nanmean(stack);
                stack = stack / nanstd(stack);

                nexttile
                imagesc(reshape(stack, [256 256]))
                title(all_animals(iAni).name, 'interpreter', 'none')
                
                axis equal tight off

            end%if Animal
        end%iAni
        colormap gray
        export_fig(['BFs\',SET.phases{iPhase},'_',SET.layers{iLayer}], '-pdf')
    end%iLayer
end%iPhase
















