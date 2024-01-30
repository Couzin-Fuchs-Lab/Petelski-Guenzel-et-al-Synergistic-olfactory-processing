%% Calcium Imaging Analysis Operator (manually select ROIs)
% CIAO_02b_selectROIs allows you add a manual selevtion regions of
% interest.
% We assume that there is a 03_Data_Processed folder within the main animal
% folder that contains a mat-file for each stimulus, the meta_info.mat
% file, and the img_segmentation.mat file.
%
% Version: 03-Mar-22 (R2023a)

%% Clean and prepare environment
clc; clear all; close all
disp('----- Ciao, welcome back! -----')
% Add paths:
% --- Current folder and subfolders
addpath(genpath(pwd))

%% Iterate over all animals and their corresponding trials
% Get path to raw data
path2animal = uigetdir(pwd,'Select folder listing all animals');

% If there is a refined version that you want to overwrite, set the
% following variable to true, if it is set to false the previous selection
% can be further refined.
overwrite_previous = true;

% Get overview of animal folders
curr.dir.all = dir(path2animal);

% Load classifier
classifier = load('locust_AL_glomeruli.mat');

%% Get distribution of zScores
relativeScore = [];
hFig = figure('units', 'normalized', 'Position', [0 0 1 1]);
% Iterate over all animals
for iAni = 1:size(curr.dir.all,1)
    clearvars -except iAni hFig relativeScore curr overwrite* path2animal classifier

    % Check whether this is what we are looking for, i.e. whether the
    % string "Animal" is present
    if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir

        % Load meta and segmentation
        curr.path.segmentation = [path2animal,'\',curr.dir.all(iAni).name,'\03_Data_Processed\'];
        load([curr.path.segmentation,'img_segmentation.mat']);

        % Get data to predict selection
        D = [...
            Segmentation.summary_stats.pocket_Avg_img(:),...
            Segmentation.summary_stats.pocket_Std_img(:),...
            Segmentation.summary_stats.pocket_Max_img(:),...
            Segmentation.summary_stats.pocket_Corr_img(:)];
        D = log(sqrt(abs(D)));
        img_label = classifier.netClassifier(D')>0;
        img_label = reshape(img_label, size(Segmentation.pockets_labeled));

        Segmentation.autoROI = img_label>0;
        save([curr.path.segmentation,'img_segmentation.mat'], 'Segmentation');

        figure(hFig)
        clf
        imagesc(sqrt(Segmentation.summary_stats.pocket_Std_img))
        hold on
        [B, L] = bwboundaries(Segmentation.autoROI>0);
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:, 2), boundary(:, 1), 'w:', 'LineWidth', 1)
        end%k
        axis equal off
        title(curr.dir.all(iAni).name)
        colormap gray
        drawnow

    end%if animal
end%iAni