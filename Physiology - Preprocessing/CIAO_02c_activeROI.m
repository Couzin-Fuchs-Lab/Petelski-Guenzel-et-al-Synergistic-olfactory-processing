%% Calcium Imaging Analysis Operator (get active pocktes)
% CIAO_02e_activeROI allows you to identify active pockets based on a
% global image threshold using Otsu's method. We determine a threshold and
% create a corresponding mask for each stimulus.
% We assume that there is a 03_Data_Processed folder within the main animal
% folder that contains a mat-file for each stimulus, the meta_info.mat
% file, and the img_segmentation.mat file.
%
% Version: 18-April-23 (R2023a)

%% Clean and prepare environment
clc; clear all; close all
disp('----- Ciao, welcome back! -----')
% Add paths:
% --- Current folder and subfolders
addpath(genpath(pwd))

%% Iterate over all animals and their corresponding trials
% Get path to raw data
path2animal = uigetdir(pwd,'Select folder listing all animals');

% Get overview of animal folders
curr.dir.all = dir(path2animal);

%% Iterate over all animals and pool data
for iAni = 1:size(curr.dir.all,1)

    % Check whether this is what we are looking for, i.e. whether the
    % string "Animal" is present
    if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir

        % Load meta and segmentation
        curr.path.segmentation = [path2animal,'\',curr.dir.all(iAni).name,'\03_Data_Processed\'];
        if isfile([curr.path.segmentation,'img_segmentation.mat'])
            load([curr.path.segmentation,'img_segmentation.mat']);
            load([curr.path.segmentation,'meta_info.mat']);

            % Get a list of all pockets in overall region
            if isfield(Segmentation, 'manualROI')
                pocket_list = unique(Segmentation.pockets_labeled.*double(Segmentation.manualROI));
            else
                if isfield(Segmentation, 'autoROI')
                    pocket_list = unique(Segmentation.pockets_labeled.*Segmentation.autoROI);
                else
                    pocket_list = unique(Segmentation.pockets_labeled.*Segmentation.pockets_bw);
                end%if auto
            end%if manual
            pocket_list(pocket_list==0) = [];


            % Determine threshold
            Segmentation.active_threshold = zeros(length(meta.unique_stim), 1);
            % Prepare results
            Segmentation.pocket_active_img = zeros([size(Segmentation.pockets_labeled), length(meta.unique_stim)]);

            % Preallocation
            pool_tc = zeros(length(pocket_list), meta.duration);
            pool_avg = zeros(length(pocket_list), length(meta.unique_stim));

            % Iterate over stimuli and get each pocket's average value.
            for iStim = 1:length(meta.unique_stim)
                % Check whether everything is there
                if isfile([curr.path.segmentation, meta.unique_stim{iStim},'.mat'])

                    % Load data
                    dat = load([curr.path.segmentation, meta.unique_stim{iStim},'.mat']);
                    dat = dat.ImageStream.(meta.unique_stim{iStim});
                    dat_reshape = reshape(dat, [size(dat,1) * size(dat,2), size(dat,3)]);

                    % Iterate over all pockets and get their time courses and
                    % the average value during the active window.
                    for iPocket = 1:length(pocket_list)
                        % Get time course
                        idx = find(Segmentation.pockets_labeled == pocket_list(iPocket));
                        pool_tc(iPocket,:) = mean(dat_reshape(idx,:));
                        pool_avg(iPocket,iStim) = mean(pool_tc(iPocket, :));
                    end%iPocket

                end%if is file
            end%iStim

            % Get a threshold that sorts pockets into active and not
            % active. Make sure to not include the reference (if it was
            % subtracted)
            if ~isempty(meta.ref_stim)
                Segmentation.active_threshold = graythresh(pool_avg(:, [find(~strcmp(meta.unique_stim, meta.ref_stim))]));
            else
                Segmentation.active_threshold = graythresh(pool_avg(:));
            end

            % Iterate over stimuli again and apply the threshold
            for iStim = 1:length(meta.unique_stim)
                % Check whether everything is there
                if isfile([curr.path.segmentation, meta.unique_stim{iStim},'.mat'])

                    % % Get all pockets that are above the threshold
                    % idx_active = find(pool_avg(:,iStim)>Segmentation.active_threshold);

                    % Get all pockets that are above the threshold
                    idx_active = find(pool_avg(:,iStim)>graythresh(pool_avg(:,iStim)));

                    % Iterate over all pockets and create a mask
                    pocket_active_img = zeros(size(Segmentation.pockets_labeled));
                    for iP = 1:length(idx_active)
                        % Get time course
                        idx = find(Segmentation.pockets_labeled == pocket_list(idx_active(iP)));
                        pocket_active_img(idx) = 1;
                    end%iP
                    Segmentation.pocket_active_img(:,:,iStim) = pocket_active_img;
                end%if is file
            end%iStim

            % Save identification of active pockets
            save([curr.path.segmentation,'img_segmentation.mat'], 'Segmentation');
        end%if segmented
    end%if animal
end%iAni
clc; clear all; close all
disp('----- Ciao, see you next time! -----')


