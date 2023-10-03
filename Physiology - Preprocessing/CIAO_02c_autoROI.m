%% Calcium Imaging Analysis Operator (automatically determin overall ROI)
% CIAO_02c_autoROI applies k-medoids clustering to all granules across all 
% stimuli to determine which granule was active (high z-score) and which 
% was not.
% We assume that there is a 03_Data_Processed folder within the main animal
% folder that contains a mat-file for each stimulus, the meta_info.mat
% file, and the img_segmentation.mat file.
%
% Version: 03-Mar-23 (R2022a)

%% Clean and prepare environment
clc; clear all; close all
disp('----- Ciao, welcome back! -----')
% Add paths:
% --- Current folder and subfolders
addpath(genpath(pwd))

%% Iterate over all animals and their corresponding trials
% Get path to raw data
path2animal = uigetdir(pwd,'Select folder listing all animals');

% Set into how many clusters the data should be separated. Note that only
% the one with the highest z-score will be used
nClusters = 2;

% Get overview of animal folders
curr.dir.all = dir(path2animal);

%% Get distribution of zScores
% Iterate over all animals
hFig = figure('Units','normalized','Position',[0.25,0.25,0.5,0.5]);
for iAni = 1:size(curr.dir.all,1)

    % Check whether this is what we are looking for, i.e. whether the
    % string "Animal" is present
    if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir

        % Load meta and segmentation
        curr.path.processed = [path2animal,'\',curr.dir.all(iAni).name,'\03_Data_Processed\'];
        load([curr.path.processed,'meta_info.mat']);
        load([curr.path.processed,'img_segmentation.mat']);

        % Get list of all granules
        pocket_list = unique(Segmentation.pockets_labeled);
        pocket_list(pocket_list==0) = [];

        % Preallocate a lot od space to get traces for all granules from
        % all stimuli
        granulePool_info = nan(size(meta.unique_stim,1)*length(pocket_list), 2);
        granulePool_data = nan(size(meta.unique_stim,1)*length(pocket_list), meta.duration);

        % Iterate over all stimuli
        cnt = 1;
        for iStim = 1:size(meta.unique_stim,1)
            % Load data
            curr_data = load([curr.path.processed, meta.unique_stim{iStim}, '.mat']);
            % Reshape data for easy indexing
            curr_data = reshape(curr_data.ImageStream.(meta.unique_stim{iStim}), [numel(meta.Mask), meta.duration]);
            %Iterate over all granules
            for iGran = 1:length(pocket_list)
                % Get avg time course for the current granule
                idx = find(Segmentation.pockets_labeled == pocket_list(iGran));
                avg = nanmean(curr_data(idx,:));
                % Normalize
                avg = normalize(avg,'range');
                % Store in the big table
                granulePool_data(cnt,:) = avg;
                granulePool_info(cnt,:) = [iStim, pocket_list(iGran)];
                cnt = cnt+1;
            end%iGran
            clear avg idx curr_data
        end%iStim

        % Perform kmedoids clustering to sort granules into active and
        % not-active
        idx_cluster = kmedoids(granulePool_data, 2,...
            'Distance','correlation',...
            'Replicates',10,...
            'Options',statset('MaxIter',1e6,'UseParallel',true));

        % Get which cluster has the highest z-score
        zScores = nan(nClusters,1);
        for iC = 1:nClusters
            if sum(idx_cluster==iC) > 1
                avg = mean(granulePool_data(idx_cluster==iC,:));
                zScores(iC) = (mean(avg(meta.activeRegion)) - mean(avg(meta.baseline))) / std(avg(meta.baseline));
            end
        end
        [~,best_cluster] = max(zScores);

        % Kick out all clusters that where not active
        invalid_cluster = 1:nClusters;
        invalid_cluster(best_cluster) = [];
        idx = [];
        for iC = 1:length(invalid_cluster)
            idx = [idx; find(idx_cluster == invalid_cluster(iC))];
        end
        granulePool_info(idx,:) = [];

        % Now, create a mask that is only valid for granules that were
        % active at least once
        Segmentation.ROI = nan([size(meta.Mask),size(meta.unique_stim,1)]);
        for iStim = 1:size(meta.unique_stim,1)
            ROI = zeros(numel(meta.Mask),1);
            granulePool_valid = unique(granulePool_info(granulePool_info(:,1)==iStim,2));
            for iGran = 1:length(granulePool_valid)
                idx = find(Segmentation.pockets_labeled == granulePool_valid(iGran));
                ROI(idx) = 1;
            end%iGran
            Segmentation.ROI(:,:,iStim) = reshape(ROI, size(meta.Mask));
        end%iStim

        % Show what is going on
        figure(hFig); clf;
        subplot(1,2,1)
        imagesc(nanmean(Segmentation.ROI,3))
        axis equal off
        subplot(1,2,2)
        imagesc(sum(Segmentation.ROI,3)>=size(meta.unique_stim,1)-1)
        axis equal off
        colormap gray
        drawnow

        % Save ROI
        Segmentation.autoROI = sum(Segmentation.ROI,3)>=size(meta.unique_stim,1)-1;
        save([curr.path.processed,'img_segmentation.mat'], 'Segmentation');

        clearvars -except iAni hFig curr nClusters path2animal
    end%if animal
end%iAni


%% Clean and prepare environment
clc; clear all; close all
disp('----- Ciao, see you next time! -----')