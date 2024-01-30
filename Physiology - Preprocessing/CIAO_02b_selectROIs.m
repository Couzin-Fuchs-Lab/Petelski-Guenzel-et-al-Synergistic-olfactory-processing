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
overwrite_previous = false;

% Get overview of animal folders
curr.dir.all = dir(path2animal);

%% Get distribution of zScores
relativeScore = [];
hFig = figure('units', 'normalized', 'Position', [0 0 1 1]);
% Iterate over all animals
for iAni = 1:size(curr.dir.all,1)
    clearvars -except iAni hFig relativeScore curr overwrite* path2animal

    % Check whether this is what we are looking for, i.e. whether the
    % string "Animal" is present
    if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir

        % Load meta and segmentation
        curr.path.segmentation = [path2animal,'\',curr.dir.all(iAni).name,'\03_Data_Processed\'];
        load([curr.path.segmentation,'img_segmentation.mat']);
        load([curr.path.segmentation,'meta_info.mat']);

        % Get image to click on.
        img = mean(Segmentation.summary_stats.pocket_Corr_img,3);
        img = cat(3, cat(3, img, zeros(size(Segmentation.pockets_labeled))), zeros(size(Segmentation.pockets_labeled)));

        % Get brightfield
        BF = zeros(size(img,1),size(img,2), length(meta.unique_stim));
        for iStim = 1:length(meta.unique_stim)
            if isfile([curr.path.segmentation, meta.unique_stim{iStim},'.mat'])
                dat = load([curr.path.segmentation, meta.unique_stim{iStim},'.mat']);
                BF(:,:,iStim) = std(dat.ImageStream.(meta.unique_stim{iStim}), [], 3);
                clear dat
            end%if isfile
        end%iStim
        % Get mx projection
        BF = nanmean(BF,3); clear dat
        % --- Normalize
        BF = BF - nanmin(BF(:));
        BF = log(sqrt(BF));
        % BF = BF - nanmean(BF(:));
        % BF = BF / nanstd(BF(:));

        % Check whether to overwrite or not
        if ~overwrite_previous && isfield(Segmentation, 'manualROI')
            img_label = Segmentation.manualROI;
        else
            % --- Empty matrix to keep track of annotation
            img_label = zeros(size(Segmentation.pockets_labeled,1), size(Segmentation.pockets_labeled,2));
            % Get coarse selection
            hFig2 = figure('units', 'normalized', 'Position', [0.15 0.15 0.7 0.7]);
            set(gcf, 'name', curr.dir.all(iAni).name)
            imagesc(BF)
            colormap gray
            axis equal off
            roi = drawpolygon;
            currRegion = roi.createMask;
            currRegion = find(currRegion);
            for iPoint = 1:length(currRegion)
                % Get current label
                curr_label = Segmentation.pockets_labeled(currRegion(iPoint));
                curr_label_idx = find(Segmentation.pockets_labeled == curr_label);
                % Update image
                if img_label(curr_label_idx(1)) == 1 % --------- Do nothing
                else % ------------------------------------------ Add a certain region
                    img_label(curr_label_idx) = 1;
                    BF(curr_label_idx) = BF(curr_label_idx)*2;
                end
            end
            img(:,:,3) = img_label;
            close(hFig2)
        end%if overwrite
        img(:,:,3) = img_label;

        figure(hFig); clf; set(gcf, 'name', curr.dir.all(iAni).name)
        subplot(1,2,1); hold on
        [B, L] = bwboundaries(img_label>0);
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:, 2), boundary(:, 1), 'w:', 'LineWidth', 1)
        end%k
        imagesc(BF)
        axis equal off
        subplot(1,2,2); hold on
        [B, L] = bwboundaries(img_label>0);
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:, 2), boundary(:, 1), 'w:', 'LineWidth', 1)
        end%k
        imagesc(img)
        axis equal off

        while true
            % Click on region
            curr_point = round(ginput(1), 0);

            % Check whether to stop
            if any(curr_point>max(size(img))) || any(curr_point<min(size(img)))
                break
            end%if stop

            % Get current label
            curr_label = Segmentation.pockets_labeled(curr_point(2), curr_point(1));
            curr_label_idx = find(Segmentation.pockets_labeled == curr_label);

            % Update image
            if img_label(curr_label_idx(1)) == 1 % --------- Take out a certain region
                img_label(curr_label_idx) = 0;
                BF(curr_label_idx) = BF(curr_label_idx)/2;
            else % ------------------------------------------ Add a certain region
                img_label(curr_label_idx) = 1;
                BF(curr_label_idx) = BF(curr_label_idx)*2;
            end

            % FILL GAPS
            img(:,:,3) = img_label;

            % Update figure
            figure(hFig)
            clf

            subplot(1,2,1)
            imagesc(BF); hold on
            [B, L] = bwboundaries(img_label>0);
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:, 2), boundary(:, 1), 'w:', 'LineWidth', 1)
            end%k
            axis equal off

            subplot(1,2,2)
            imagesc(img); hold on
            [B, L] = bwboundaries(img_label>0);
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:, 2), boundary(:, 1), 'w:', 'LineWidth', 1)
            end%k
            axis equal off

        end%while

        % Save manual annotation
        Segmentation.manualROI = img_label>0;
        save([curr.path.segmentation,'img_segmentation.mat'], 'Segmentation');

    end%if animal
end%iAni
clc; clear all; close all
disp('----- Ciao, see you next time! -----')