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
% response motifs, response vectors and other diagrams describing the
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
SET.layers = {'layer01_top', 'layer02_middle'};

% Hardcode active region (1:18)
SET.active_region = [7 13];

% Set framerate
SET.fps = 1;

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
SET.maxk_eval = 25;
SET.kmeans.dist_RespProfile = 'sqeuclidean';
SET.kmeans.dist_AvgValue = 'sqeuclidean';
SET.dist_pca = 'correlation';
SET.pca_explained = 99;

% Set how to sort response motifs
% 'manual' or 'rise_time'
SET.sort_respprofiles = 'manual';

% Cosmetics
SET.Colors = [...
    102,166,030;...BERRYLEAF
    217,095,002;...Col
    231,041,138;...ZBERRYLEAF_COL
    ]/255;

% Colormap for heatmaps
SET.HeatmapColors.gregarious = Confocal_SubFcn.ColMapPlasma(1000);
SET.HeatmapColors.solitarious = [...
    SET.HeatmapColors.gregarious(:,2),...
    SET.HeatmapColors.gregarious(:,1),...
    SET.HeatmapColors.gregarious(:,3)];
SET.HeatmapColors.other = gray(1000);


%% Pool data

% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % And iterate over all imaging layers
    for iLayer = 1:length(SET.layers)

        % Make a new folder for each layer
        mkdir(['ConfocalAnalysis_OdorSpace\', SET.layers{iLayer}])

        % Contruct the path
        basepath = [SET.main_path, SET.phases{iPhase},'\', SET.layers{iLayer},'\'];
        all_animals = dir(basepath);

        % Prepare data structures
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).stim_dist_val = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).stim_dist_act = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_resp = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_RespVec = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_RespProfile = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_value = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).active_mask = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).dist2center = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).VennData = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_RespProfile = [];
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_roi_mask = {};
        cnt = 1;

        % Iterate over all animals
        for iAni = 1:size(all_animals)

            % Check whether the current folder is belonging to an animal
            if contains(all_animals(iAni).name, '_Animal')

                % Get the correct path
                path2data = [basepath, all_animals(iAni).name, '\03_Data_Processed\'];

                % Load meta and segmentation
                load([path2data, 'meta_info.mat'])
                load([path2data, 'img_segmentation.mat'])

                % Get list of pockets. Take all pockets that were active
                % at least once
                mask_list = zeros(length(meta.unique_stim),1);
                for iStim = 1:length(SET.StimList)
                    mask_list = mask_list+ (meta.unique_stim == SET.StimList{iStim});
                end
                mask = sum(Segmentation.pocket_active_img(:,:,find(mask_list)),3);
                mask = mask>0;
                pocket_list = unique(Segmentation.pockets_labeled.*double(mask));
                pocket_list(pocket_list==0) = [];

                % Reshape mask for easy access
                mask_reshape = reshape(Segmentation.pocket_active_img(:,:,find(mask_list)), [size(Segmentation.pocket_active_img,1)*size(Segmentation.pocket_active_img,2),sum(mask_list)]);

                % Get each pocket's centroid
                curr_centroid = nan(length(pocket_list),3);
                for iP = 1:length(pocket_list)
                    % Get all pixels that belong to the pocket
                    idx_k = find(Segmentation.pockets_labeled==pocket_list(iP));
                    C = regionprops(Segmentation.pockets_labeled==pocket_list(iP), 'Centroid');
                    curr_centroid(iP,1:2) = C.Centroid;
                end
                % Get the distance to the overall centroid
                C = regionprops(Segmentation.manualROI, 'Centroid'); C = C.Centroid;
                helper = curr_centroid(:,1:2)-C;
                curr_centroid(:,3) = sqrt(sum(helper'.*helper'))';
                % Normalize
                curr_centroid(:,3) = curr_centroid(:,3)/max(curr_centroid(:,3));

                % Prepare data structures
                curr_tc = [];
                curr_active = nan(length(pocket_list), length(SET.StimList));
                curr_avg = nan(length(pocket_list), length(SET.StimList));

                % Iterate over stimuli and get data
                for iStim = 1:length(SET.StimList)
                    % Load data
                    curr_data = load([path2data, SET.StimList{iStim}, '.mat']);
                    curr_data = curr_data.ImageStream.(SET.StimList{iStim});
                    % Reshape
                    curr_data_reshape = reshape(curr_data, [size(curr_data,1)*size(curr_data,2), size(curr_data,3)]);
                    % Cut
                    curr_data_reshape = curr_data_reshape(:,meta.baseline(1):meta.activeRegion(end)+length(meta.baseline)-1);
                    % Preallocation
                    avg_RespVec = nan(length(pocket_list), size(curr_data_reshape,2));
                    % Iterate over all pockets
                    for iP = 1:length(pocket_list)
                        % Get all pixels that belong to the pocket
                        idx_k = find(Segmentation.pockets_labeled==pocket_list(iP));
                        % Get mask whether it was active
                        curr_active(iP,iStim) = nanmean(mask_reshape(idx_k,iStim));
                        % Get the average time course
                        avg_RespVec(iP,:) = nanmean(curr_data_reshape(idx_k,:));
                        % Get avg activity during active window
                        curr_avg(iP,iStim) = nanmean(avg_RespVec(iP,SET.active_region(1):SET.active_region(2)));
                        % Keep track of some infos
                        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_RespProfile = [PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_RespProfile; iPhase, iAni, iStim, iP];
                    end%iP
                    curr_tc = [curr_tc, avg_RespVec];
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_RespProfile = [PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_RespProfile; avg_RespVec];
                end%iStim

                % Get distances between stimuli
                % --- Based on  avg value
                dist = tril(squareform(pdist(curr_avg')));
                dist(dist==0)=NaN;
                dist = dist/nanmean(dist(:));
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).stim_dist_val = cat(3,...
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).stim_dist_val,...
                    dist);
                % --- Based on active mask
                dist = triu(squareform(pdist(curr_active')));
                dist(dist==0)=NaN;
                dist = dist/nanmean(dist(:));
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).stim_dist_act = cat(3,...
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).stim_dist_act,...
                    dist);

                % Extract info on how uniquely ROIs respond: Number of stimuli
                % a ROI is responding to
                curr_n_resp = hist(sum(curr_active,2),1:length(SET.StimList));
                curr_n_resp = curr_n_resp/sum(curr_n_resp);
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_resp = [...
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_resp;...
                    curr_n_resp];

                % Extract info on how uniquely ROIs respond: Number of ROIs
                % uniquely responding to one stimulus class
                stim_yes = {'ZHAE',         'COL',            'ZHAE',              'COL',                'ZHAECOL',        'ZHAE',                  'COL'};
                stim_no = {{'COL';'COL'},  {'ZHAE';'ZHAE'},  {'COL'; 'ZHAECOL'},  {'ZHAE'; 'ZHAECOL'},  {'ZHAE'; 'COL'},  {'ZHAECOL'; 'ZHAECOL'},  {'ZHAECOL'; 'ZHAECOL'}};
                curr_unique = nan(1,length(stim_yes));
                curr_unique_info = cell(2,length(stim_yes));
                % Iterate over all stimuli
                for iStim = 1:length(stim_yes)
                    % Get index position of the + stimulus
                    idx_yes = find(strcmp(SET.StimList, stim_yes{iStim}));
                    idx_no1 = find(strcmp(SET.StimList, stim_no{iStim}{1}));
                    idx_no2 = find(strcmp(SET.StimList, stim_no{iStim}{2}));
                    % Get the proportion of regions that responded to the current
                    % stimulus but not to the other two
                    curr_unique(1,iStim) = length(find(curr_active(:,idx_yes)==1 & curr_active(:,idx_no1)==0 & curr_active(:,idx_no2)==0)) / size(curr_active,1);
                    % Keep track
                    curr_unique_info{1,iStim} = stim_yes{iStim};
                    curr_unique_info{2,iStim} = stim_no{iStim};
                end

                % Extract data for Venn diagram. Similar to data exraction
                % above
                curr_Venn = zeros(1,7);
                % --- only zhae
                idx_yes = find(strcmp(SET.StimList, 'ZHAE'));
                idx_no = find(~strcmp(SET.StimList, 'ZHAE'));
                curr_Venn(1,1) = length(find(curr_active(:,idx_yes)==1 & curr_active(:,idx_no(1))==0 & curr_active(:,idx_no(2))==0));
                % --- only colony
                idx_yes = find(strcmp(SET.StimList, 'COL'));
                idx_no = find(~strcmp(SET.StimList, 'COL'));
                curr_Venn(1,2) = length(find(curr_active(:,idx_yes)==1 & curr_active(:,idx_no(1))==0 & curr_active(:,idx_no(2))==0));
                % --- only mixture
                idx_yes = find(strcmp(SET.StimList, 'ZHAECOL'));
                idx_no = find(~strcmp(SET.StimList, 'ZHAECOL'));
                curr_Venn(1,3) = length(find(curr_active(:,idx_yes)==1 & curr_active(:,idx_no(1))==0 & curr_active(:,idx_no(2))==0));
                % --- only zhae and colony
                idx_yes = find(~strcmp(SET.StimList, 'ZHAECOL'));
                idx_no = find(strcmp(SET.StimList, 'ZHAECOL'));
                curr_Venn(1,4) = length(find(curr_active(:,idx_yes(1))==1 & curr_active(:,idx_yes(2))==1 & curr_active(:,idx_no(1))==0));
                % --- only zhae and mixture
                idx_yes = find(~strcmp(SET.StimList, 'COL'));
                idx_no = find(strcmp(SET.StimList, 'COL'));
                curr_Venn(1,5) = length(find(curr_active(:,idx_yes(1))==1 & curr_active(:,idx_yes(2))==1 & curr_active(:,idx_no(1))==0));
                % --- only colony and mixture
                idx_yes = find(~strcmp(SET.StimList, 'ZHAE'));
                idx_no = find(strcmp(SET.StimList, 'ZHAE'));
                curr_Venn(1,6) = length(find(curr_active(:,idx_yes(1))==1 & curr_active(:,idx_yes(2))==1 & curr_active(:,idx_no(1))==0));
                % --- zhae and colony and mixture
                idx_yes = [find(strcmp(SET.StimList, 'ZHAE')), find(strcmp(SET.StimList, 'COL')), find(strcmp(SET.StimList, 'ZHAECOL'))];
                curr_Venn(1,7) = length(find(curr_active(:,idx_yes(1))==1 & curr_active(:,idx_yes(2))==1 & curr_active(:,idx_yes(3))==0));
                % Normalize
                curr_Venn = curr_Venn/sum(curr_Venn);

                % Pool data
                % --- number of unique stimuli a ROI is responding to
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique = [...
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique;...
                    curr_unique];
                % --- info on the number of unique stimuli a ROI is responding to
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique_info = curr_unique_info;
                % --- avg tc
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_RespVec = [...
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_RespVec;...
                    curr_tc];
                % --- avg activity
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_value = [...
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_value;...
                    curr_avg];
                % --- activity mask
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).active_mask = [...
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).active_mask;...
                    curr_active];
                % --- information
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info = [...
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info; ...
                    ones(length(curr_avg),1)*iAni, pocket_list(:)];
                % --- pocket distance to center
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).dist2center = [...
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).dist2center; curr_centroid];
                % --- Venn diagram proportions
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).VennData = [...
                    PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).VennData; curr_Venn];
                % --- info on the mask of ROIs
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_roi_mask{cnt,1} = iAni;
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_roi_mask{cnt,2} = all_animals(iAni).name;
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_roi_mask{cnt,3} = pocket_list;
                PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_roi_mask{cnt,4} = Segmentation.pockets_labeled;
                cnt = cnt+1;

                % Clean
                clearvars -except PoolData iStim iAni iPhase iLayer SET basepath all_animals cnt

            end%if animal
        end%iAni
    end%iLayer
end%iPhase

%% Pool phases and get response motif IDs
% -------------------------------------------------------------------------
for iLayer = 1:length(SET.layers)
    active_mask.(SET.layers{iLayer}) = [];
    dat_RespProfile.(SET.layers{iLayer}) = [];
    dat_AvgVal.(SET.layers{iLayer}) = [];
    dat_TC.(SET.layers{iLayer}) = [];
    anis_RespProfile.(SET.layers{iLayer}) = [];
    anis_AvgVal.(SET.layers{iLayer}) = [];
    phase_info.(SET.layers{iLayer}) = [];
    respprofile_info.(SET.layers{iLayer}) = [];
    dist2center.(SET.layers{iLayer}) = [];
    for iPhase = 1:length(SET.phases)
        active_mask.(SET.layers{iLayer}) =      [active_mask.(SET.layers{iLayer});      PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).active_mask];
        dat_RespProfile.(SET.layers{iLayer}) =  [dat_RespProfile.(SET.layers{iLayer});  PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_RespProfile];
        dat_AvgVal.(SET.layers{iLayer}) =       [dat_AvgVal.(SET.layers{iLayer});       PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_value];
        dat_TC.(SET.layers{iLayer}) =           [dat_TC.(SET.layers{iLayer});           PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_RespVec];
        anis_RespProfile.(SET.layers{iLayer}) = [anis_RespProfile.(SET.layers{iLayer}); PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_RespProfile(:,2)*(2*((iPhase==1)-0.5))];
        anis_AvgVal.(SET.layers{iLayer}) =      [anis_AvgVal.(SET.layers{iLayer});      PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info(:,1)*(2*((iPhase==1)-0.5))];
        phase_info.(SET.layers{iLayer}) =       [phase_info.(SET.layers{iLayer});       ones(size(PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).avg_value,1),1)*iPhase];
        respprofile_info.(SET.layers{iLayer}) = [respprofile_info.(SET.layers{iLayer}); PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_RespProfile];
        dist2center.(SET.layers{iLayer}) =      [dist2center.(SET.layers{iLayer});      PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).dist2center(:,3)];
    end%iPhase
    % Keep original version (not normalized) for later
    dat_TC_original.(SET.layers{iLayer}) = dat_TC.(SET.layers{iLayer});
    dat_AvgVal_original.(SET.layers{iLayer}) = dat_AvgVal.(SET.layers{iLayer});
    dat_RespProfile_original.(SET.layers{iLayer}) = dat_RespProfile.(SET.layers{iLayer});
end%iLayer



%% Normalize data from each animal
for iLayer = 1:length(SET.layers)
    % --- response motif
    ani_list = unique(anis_RespProfile.(SET.layers{iLayer}));
    for iAni = 1:length(ani_list)
        idx_k = find(anis_RespProfile.(SET.layers{iLayer}) == ani_list(iAni));
        curr_dat = dat_RespProfile.(SET.layers{iLayer})(idx_k,:);
        curr_dat = curr_dat - min(curr_dat(:));
        curr_dat = sqrt(curr_dat);
        curr_dat = curr_dat - nanmean(curr_dat(:));
        curr_dat = curr_dat / nanstd(curr_dat(:));
        dat_RespProfile.(SET.layers{iLayer})(idx_k,:) = curr_dat;
    end%iAni

    % --- Average value
    ani_list = unique(anis_AvgVal.(SET.layers{iLayer}));
    for iAni = 1:length(ani_list)
        idx_k = find(anis_AvgVal.(SET.layers{iLayer}) == ani_list(iAni));
        curr_dat = dat_AvgVal.(SET.layers{iLayer})(idx_k,:);
        curr_dat = curr_dat - min(curr_dat(:));
        curr_dat = sqrt(curr_dat);
        curr_dat = curr_dat - nanmean(curr_dat(:));
        curr_dat = curr_dat / nanstd(curr_dat(:));
        dat_AvgVal.(SET.layers{iLayer})(idx_k,:) = curr_dat;
        curr_dat = dat_TC.(SET.layers{iLayer})(idx_k,:);
        curr_dat = curr_dat - min(curr_dat(:));
        curr_dat = sqrt(curr_dat);
        curr_dat = curr_dat - nanmean(curr_dat(:));
        curr_dat = curr_dat / nanstd(curr_dat(:));
        dat_TC.(SET.layers{iLayer})(idx_k,:) = curr_dat;
    end%iAni
end%iLayer


%% Cluster response motifs
% Here we cluster all time courses from all ROIs responsing to all stimuli
% to check whether there are unique response motifs (e.g., slow rise,
% long peak, inhibition, etc)

% Iterate over all layers
for iLayer = 1:length(SET.layers)

    % Normalize each row
    dat = normalize(dat_RespProfile.(SET.layers{iLayer}),2);

    % PCA preprocessing
    opts = statset('pca'); opts.MaxIter = 1e6; opts.TolFun = 1e-8; opts.TolFunX = 1e-8;
    [~,dat,~,~,explained,~] = pca(dat, 'Economy', false, 'Options', opts);
    ind = find(cumsum(explained)>SET.pca_explained, 1);
    dat = dat(:,1:ind);

    % Run multiple times to identify the optimal number of clusters
    VRC = nan(SET.maxk_eval,1);
    hWait = waitbar(1, 'Evaluating clusters ...');
    for k = 2:SET.maxk_eval
        waitbar(k/SET.maxk_eval, hWait)
        rng(1234)
        [idx,C,sumd] = kmeans(dat, k,...
            "Distance", SET.kmeans.dist_RespProfile,...
            'MaxIter', SET.kmeans.MaxIter,...
            'OnlinePhase', SET.kmeans.OnlinePhase,...
            'Replicates', SET.kmeans.Replicates,...
            'Options', statset('UseParallel',1));
        % Get the variance ratio criterion (VRC)
        VRC(k) = Confocal_SubFcn.VarianceRatioCriterion(dat, k, idx, C);
    end%k
    % Get the optimal number of clusters based on the ellbow point of the
    % exponentially decaying VRC
    % First, fit the exponential function to the data
    [xData, yData] = prepareCurveData(2:SET.maxk_eval, VRC(2:end));
    ft = fittype( 'exp2' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    fitresult = fit(xData, yData, ft, opts);
    % Next, Interpolate the x values to get a smooth curve and use this to get
    % the fitted line
    xData_interp = linspace(xData(1), xData(end), 1000);
    VRC_fit = fitresult.a*exp(fitresult.b*xData_interp) + fitresult.c*exp(fitresult.d*xData_interp);
    % Now, get the elbow point which is defined as the point on the curve
    % that is farthest from the line segment connecting the first and last
    %  points on the curve.
    elbow_x_value = Confocal_SubFcn.find_elbow_point(xData_interp, VRC_fit);
    % For the optimal number of clusters, add 1 (to besure)
    kmeans_result.best_k = round(elbow_x_value,0);
    % Repeat clustering with optimal number and more iteration
    rng(1234)
    [ID, C] = kmeans(dat, kmeans_result.best_k,...
        "Distance", SET.kmeans.dist_RespProfile,...
        'MaxIter', SET.kmeans.MaxIter,...
        'OnlinePhase', SET.kmeans.OnlinePhase,...
        'Replicates', SET.kmeans.Replicates,...
        'Options', statset('UseParallel',1));
    close(hWait)

    % Sort centroids either manually or by rise time
    switch SET.sort_respprofiles
        case 'manual'
            % Depict all clusters
            hFig1 = figure('units', 'normalized', 'Position', [0 0 0.5 1], 'Color', 'w');
            hFig2 = figure('units', 'normalized', 'Position', [0.5 0 0.5 1], 'Color', 'w');
            for iK = 1:kmeans_result.best_k
                idx = find(ID==iK);
                figure(hFig1)
                nexttile
                plot(mean(dat_RespProfile_original.(SET.layers{iLayer})(idx,:)), 'LineWidth',3)
                title(iK)
                figure(hFig2); hold on
                plot(mean(dat_RespProfile_original.(SET.layers{iLayer})(idx,:)), 'LineWidth',3)
            end
            figure(hFig2); legend;
            % Ask user
            clc
            disp('Give response motif cluster order as, e.g., [1 3 2]')
            order_idx = input('Order: ');
            % --- Sort centroids
            temp_ID = nan(size(ID));
            for iK = 1:kmeans_result.best_k
                temp_ID(ID==order_idx(iK)) = iK;
            end%iK
            ID = temp_ID;
            close(hFig1); close(hFig2)
        case 'rise_time'
            % --- Get location of maximum
            avg_val = zeros(kmeans_result.best_k,1);
            for iM = 1:kmeans_result.best_k
                [~,avg_val(iM)] = max(C(iM,:));
            end%iDec
            % --- Sort centroids
            temp_ID = nan(size(ID));
            id_list = unique(ID);
            [~,order_idx] = sort(avg_val);
            for iK = 1:kmeans_result.best_k
                temp_ID(ID==order_idx(iK)) = iK;
            end%iK
            ID = temp_ID;
    end

    % Get the corresponding "real" centroids
    for iPhase = 1:length(SET.phases)
        % Get list of animal and IDs
        ani_list = unique(respprofile_info.(SET.layers{iLayer})(respprofile_info.(SET.layers{iLayer})(:,1) == iPhase,2));
        id_list = unique(ID(respprofile_info.(SET.layers{iLayer})(:,1) == iPhase));
        avg_val = nan(kmeans_result.best_k,1);
        for iK = 1:kmeans_result.best_k
            kmeans_result.centroid{iPhase}{iK} = nan(length(ani_list), size(dat_RespProfile_original.(SET.layers{iLayer}),2));
            for iAni = 1:length(ani_list)
                idx_ani = find(respprofile_info.(SET.layers{iLayer})(:,2) == ani_list(iAni) & ID == id_list(iK));
                kmeans_result.centroid{iPhase}{iK}(iAni,:) = nanmean(dat_RespProfile_original.(SET.layers{iLayer})(idx_ani,:),1);
            end%iAni
            avg_val(iK) = nanmean(nanmean(kmeans_result.centroid{iPhase}{iK}(:,SET.active_region(1):SET.active_region(2))));
        end%iK
    end%iPhase

    % Save
    PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).best_k = kmeans_result.best_k;
    PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).ID = ID;
    PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).dat_RespProfile = dat_RespProfile_original.(SET.layers{iLayer});
    PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).centroid = kmeans_result.centroid;
    PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info = respprofile_info.(SET.layers{iLayer});

end%iLayer
% Clean up
clearvars -except PoolData SET dat_* anis_* phase_info respprofile_info active_mask dist2center


%% Cluster Response Vectors
% Here, we cluster ROIs based on their response motifs to each stimulus.
% We do not take average values only since the profile offers considerable
% information here.

% Iterate over all layers
for iLayer = 1:length(SET.layers)

    % Normalize each row
    dat = normalize(dat_AvgVal.(SET.layers{iLayer}),2);
    dat_avg = normalize(dat_AvgVal.(SET.layers{iLayer}),2);

    % PCA preprocessing
    opts = statset('pca'); opts.MaxIter = 1e6; opts.TolFun = 1e-8; opts.TolFunX = 1e-8;
    [~,dat,~,~,explained,~] = pca(dat, 'Economy', false, 'Options', opts);
    ind = find(cumsum(explained)>SET.pca_explained, 1);
    dat = dat(:,1:ind);

    % Run multiple times to identify the optimal number of clusters
    VRC = nan(SET.maxk_eval,1);
    hWait = waitbar(1, 'Evaluating clusters ...');
    % Test several cluster sizes
    for k = 2:SET.maxk_eval
        waitbar(k/SET.maxk_eval, hWait)
        rng(1234)
        [idx,C,sumd] = kmeans(dat, k,...
            "Distance", SET.kmeans.dist_AvgValue,...
            'MaxIter', SET.kmeans.MaxIter,...
            'OnlinePhase', SET.kmeans.OnlinePhase,...
            'Replicates', SET.kmeans.Replicates,...
            'Options', statset('UseParallel',1));
        % Get the variance ratio criterion (VRC)
        VRC(k) = Confocal_SubFcn.VarianceRatioCriterion(dat, k, idx, C);
    end%k
    % Get the optimal number of clusters based on the ellbow point of the
    % exponentially decaying VRC
    [xData, yData] = prepareCurveData(2:SET.maxk_eval, VRC(2:end));
    ft = fittype( 'exp2' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    fitresult = fit(xData, yData, ft, opts);
    % Next, Interpolate the x values to get a smooth curve and use this to get
    % the fitted line
    xData_interp = linspace(xData(1), xData(end), 1000);
    VRC_fit = fitresult.a*exp(fitresult.b*xData_interp) + fitresult.c*exp(fitresult.d*xData_interp);
    % Now, get the elbow point which is defined as the point on the curve
    % that is farthest from the line segment connecting the first and last
    %  points on the curve.
    elbow_x_value = Confocal_SubFcn.find_elbow_point(xData_interp, VRC_fit);
    % For the optimal number of clusters, add 1 (to besure)
    kmeans_result.best_k = round(elbow_x_value,0);
    % Repeat clustering with optimal number and more iteration
    rng(1234)
    [ID, C] = kmeans(dat, kmeans_result.best_k,...
        "Distance", SET.kmeans.dist_AvgValue,...
        'MaxIter', SET.kmeans.MaxIter,...
        'OnlinePhase', SET.kmeans.OnlinePhase,...
        'Replicates', SET.kmeans.Replicates,...
        'Options', statset('UseParallel',1));
    close(hWait)

    % Sort centroids by their correlation
    % --- Get a dendrogram based on the correlation
    D = squareform(pdist(C, SET.dist_pca));
    tree = linkage(D,'single', SET.dist_pca);
    hD = figure;
    [~,~,order_idx]=dendrogram(gca,tree);
    close (hD)
    % --- Sort centroids
    temp_ID = nan(size(ID));
    for iK = 1:kmeans_result.best_k
        temp_ID(ID==order_idx(iK)) = iK;
    end%iK
    ID = temp_ID;
    clear D tree hD temp_ID id_list iK

    % Get each animal's centroid
    % Iterate over all phases
    for iPhase = 1:length(SET.phases)
        % Get list of animal
        ani_list = unique(anis_AvgVal.(SET.layers{iLayer})(phase_info.(SET.layers{iLayer})==iPhase));
        id_list = unique(ID(phase_info.(SET.layers{iLayer})==iPhase));
        % Iterate over all clusters
        for iK = 1:kmeans_result.best_k
            % Preallocation
            kmeans_result.centroid{iPhase}{iK} = nan(length(ani_list), size(dat_AvgVal_original.(SET.layers{iLayer}),2));
            kmeans_result.centroid_pca{iPhase}{iK} = nan(length(ani_list), 2);
            % Iterate over all animal
            for iAni = 1:length(ani_list)
                idx_ani = find(anis_AvgVal.(SET.layers{iLayer}) == ani_list(iAni) & ID == id_list(iK));
                kmeans_result.centroid{iPhase}{iK}(iAni,:) =       nanmean(dat_AvgVal_original.(SET.layers{iLayer})(idx_ani,:),1);
                kmeans_result.centroid_norm{iPhase}{iK}(iAni,:) =  nanmean(dat_avg(idx_ani,:),1);
                kmeans_result.centroid_pca{iPhase}{iK}(iAni,:) =   nanmean(dat(idx_ani,1:2),1);
                dist2center_sorted{iPhase}{iK}(iAni,1) =           nanmean(dist2center.(SET.layers{iLayer})(idx_ani));
            end%iAni
        end%iK
    end%iPhase

    % Save
    % --- The extracted data
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dat_active_mask = active_mask.(SET.layers{iLayer});      % whether the ROI was actvive or not
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dat_TC = dat_TC_original.(SET.layers{iLayer});           % time traces (not normalized)
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dat_AvgVal = dat_AvgVal_original.(SET.layers{iLayer});   % average values (not normalized)
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dat_pca = dat;                      % the PCA pre-proceed data used for clustering
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dat_pca_explained = explained;      % the PCA pre-proceed data used for clustering
    % --- Each ROIs distance to the respective AL's center
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dat_dist2center = dist2center.(SET.layers{iLayer});
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dist2center_sorted = dist2center_sorted;
    % --- The clustring results
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k = kmeans_result.best_k;                % optimal number of clusters
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_ID = ID;                                      % each ROI's assigned ID
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid = kmeans_result.centroid;            % each animal's avg value centroid (not normalized)
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_norm = kmeans_result.centroid_norm;  % each animal's avg value centroid (not normalized)
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_pca = kmeans_result.centroid_pca;    % each animal's centroid in PCA space (pca-1 & pca-2)
    % --- Information on each ROI's phase (col-1) and animal (col-2)
    PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info = [phase_info.(SET.layers{iLayer}), anis_AvgVal.(SET.layers{iLayer})];

    clearvars -except PoolData SET dat_* anis_* phase_info respprofile_info active_mask dist2center iLayer
end%iLayer

% Save and clean up
save('ConfocalAnalysis_OdorSpace\OdorSpace_kmeans_PoolData.mat', 'PoolData', 'SET')
clearvars -except PoolData SET


%% Visualize the results.
% response motif vectors

% Define the range of the y-axis
ylim_range = [-4 8; -1.5 7.0];

% Iterate over both layers
for iLayer = 1:length(SET.layers)
    % One figure per layer
    hFig = figure('units', 'normalized', 'Position', [0 0 0.5 0.5]);
    % Iterate over both phases and plot the average motifs for each phase
    % in a subplot
    for iPhase = 1:length(SET.phases)
        % One subplot per phase
        subplot(1,3,iPhase) ; hold on
        % Get colors
        cols = SET.HeatmapColors.(SET.phases{iPhase});
        cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).best_k+2)),:);
        cols = cols(2:end-1,:);
        % Get time vector
        xvec = linspace(0, size(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).dat_RespProfile,2)/SET.fps, size(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).dat_RespProfile,2));
        % Get position of active window
        pos = [xvec(SET.active_region(1)), -1, xvec(SET.active_region(2))-xvec(SET.active_region(1)),3];
        % Plot active window as gray box
        rectangle('position', pos, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
        % Iterate over all motifs
        for iK = 1:PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).best_k
            % Get the average
            avg = nanmean(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).centroid{iPhase}{iK},1);
            % Get the corresponding confidence interval
            ci = bootci(5000, {@nanmean, PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).centroid{iPhase}{iK}});
            % Plot both. Turn confidence intervals into shaded areas during
            % figure post-processing
            plot(xvec, avg, 'LineWidth', 2, 'Color', cols(iK,:))
            plot(xvec, ci', 'LineWidth', 1, 'Color', cols(iK,:))
            drawnow
        end%iK
        % Cosmetics
        ylabel('\DeltaF/F (%)')
        xlabel('time (s)')
        title(SET.phases{iPhase})
    end%iPhase
    %Also, plot the motifs for both phases combined in the third subplot
    subplot(1,3,3) ; hold on
    % Use a different colormap
    cols = Confocal_SubFcn.ColMapInferno(1000);
    cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).best_k+2)),:);
    cols = cols(2:end-1,:);
    % Get the time vector
    xvec = linspace(0, size(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).dat_RespProfile,2)/SET.fps, size(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).dat_RespProfile,2));
    % Get the position of the active region
    pos = [xvec(SET.active_region(1)), -1, xvec(SET.active_region(2))-xvec(SET.active_region(1)),3];
    % Plot the active region as a gray box
    rectangle('position', pos, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
    % Iterate over all motifs
    for iK = 1:PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).best_k
        % Pool gregarious and solitarious data
        currDat = [];
        for iPhase = 1:length(SET.phases)
            currDat = [currDat; PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).centroid{iPhase}{iK}];
        end%iPhase
        % Get the average
        avg = nanmean(currDat,1);
        % And the corresponding confidence interval
        ci = bootci(5000, {@nanmean, currDat});
        % Plot both. Turn confidence intervals into shaded areas during
        % figure post-processing
        plot(xvec, avg, 'LineWidth', 2, 'Color', cols(iK,:))
        plot(xvec, ci', 'LineWidth', 1, 'Color', cols(iK,:))
        drawnow
    end%iK
    % Cosmetics
    ylim(ylim_range(iLayer,:))
    ylabel('\DeltaF/F (%)')
    xlabel('time (s)')
    title('both')
    % Save figure and close it
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\kmeans_ResponseProfile'], '-pdf', '-painters')
    close(hFig)
end%iLayer
clearvars -except PoolData SET


%% Visualize the results.
% Number stimuli a ROI responds to

% Iterate over both layers
for iLayer = 1:length(SET.layers)
    % One figure per layer
    hFig = figure('Name', 'n_resp'); hold on
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % Get colors
        col = SET.HeatmapColors.(SET.phases{iPhase}); col = col(floor(size(col,1)/2),:);
        % Get average and CI across animals
        avg = nanmean(bootstrp(5000, @nanmean, PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_resp));
        ci = bootci(5000, @nanmean, PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_resp);
        % Iterate over the number of stimuli
        for iStim = 1:length(SET.StimList)
            % Get position of current bar
            pos = [iStim+iPhase/2-1+0.25*(iPhase==1), 0, 0.25, avg(iStim)];
            x = iStim+iPhase/2-1+0.25*(iPhase==1) + 0.25/2;
            % Show bar
            rectangle('position', pos, 'FaceColor', col, 'EdgeColor', 'none')
            plot([1 1]*x, ci(:,iStim), 'k')
            % Show CIs
            plot([-0.25/2 0.25/2]+x, [1 1]*ci(1,iStim), 'k')
            plot([-0.25/2 0.25/2]+x, [1 1]*ci(2,iStim), 'k')
        end%iStim
    end%iPhase
    % Cosmetics
    xlim([0.5 length(SET.StimList)+0.5])
    ylim([0 0.5])
    set(gca, 'xtick',1:7)
    ylabel('prop. ROIs')
    xlabel('number of stimuli')
    title('Proportion of ROIs responding to a given number of stimuli')
    text(1,0.4,SET.StimList, "Interpreter","none")
    % Save
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\kmeans_n_resp'], '-pdf', '-painters')
    close(hFig)
end%iLayer
clearvars -except PoolData SET


%% Visualize the results.
% Motif combinations and how they relate to synergy

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    
    % Iterate over both phases and get all motif combinations
    for iPhase = 1:length(SET.phases)
        % Get corresponding grouping factors (response motifs)
        motifs = PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).ID(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(:,1) == iPhase,:);
        stims = PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(:,1) == iPhase,3);
        data_grouping.(SET.phases{iPhase}) = [];
        % Iterate ove all stimuli and populate the grouping variable
        for iStim = 1:length(SET.StimList)
            idx = find(stims == iStim);
            data_grouping.(SET.phases{iPhase}) = [data_grouping.(SET.phases{iPhase}), motifs(idx)];
        end%iStim
        data_grouping.(SET.phases{iPhase})(std(data_grouping.(SET.phases{iPhase}), [], 2) == 0, :) = [];
    end%iPhase  
    

    % Get all motifs
    for iC = 1: PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).best_k
        motif.(['m', num2str(iC)]) = mean(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).dat_RespProfile(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).ID==iC,:))/100;
    end%iC

    % Get the overall Bliss score resulting from all motifs
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        syn.(SET.phases{iPhase}) = nan(size(data_grouping.(SET.phases{iPhase}),1), length(motif.m1));
        for iGran = 1:size(data_grouping.(SET.phases{iPhase}),1)
            % Get motifs
            a = motif.(['m', num2str(data_grouping.(SET.phases{iPhase})(iGran,1))]);
            b = motif.(['m', num2str(data_grouping.(SET.phases{iPhase})(iGran,2))]);
            c = motif.(['m', num2str(data_grouping.(SET.phases{iPhase})(iGran,3))]);
            % Get expected and observed
            f_expected = (a/2 + b/2) - ((a/2).*(b/2));
            f_observed = c;
            % Get BS
            syn.(SET.phases{iPhase})(iGran,:) = 100*(f_observed - f_expected);
        end%iGran
    end%iPhase

    % Get the largest difference between phases in terms of the number of
    % unique motif combinations
    grouping_all = [data_grouping.(SET.phases{1}); data_grouping.(SET.phases{2})];
    grouping_info = [ones(size(data_grouping.(SET.phases{1}),1),1); ones(size(data_grouping.(SET.phases{2}),1),1)*2];
    [grouping_unique, ~, grouping_id] = unique(grouping_all, 'rows');
    % Count
    grouping_count = zeros(size(grouping_unique, 1),3);
    for iGrp = 1:size(grouping_unique, 1)
        expected_greg = size(data_grouping.(SET.phases{1}),1)/length(grouping_unique);
        expected_soli = size(data_grouping.(SET.phases{2}),1)/length(grouping_unique);
        grouping_count(iGrp,1) = length(find(grouping_id == iGrp & grouping_info == 1)) / expected_greg;
        grouping_count(iGrp,2) = length(find(grouping_id == iGrp & grouping_info == 2)) / expected_soli;

        % grouping_count(iGrp,1) = 100 * length(find(grouping_id == iGrp & grouping_info == 1)) / size(data_grouping.(SET.phases{1}),1);
        % grouping_count(iGrp,2) = 100 * length(find(grouping_id == iGrp & grouping_info == 2)) / size(data_grouping.(SET.phases{2}),1);
        grouping_count(iGrp,3) = (grouping_count(iGrp,1) - grouping_count(iGrp,2));
    end%iGrp

    % Get the nine most substantial differences and their corresponding
    [sort_grouping_diff, sort_grouping_ind] = sort(grouping_count(:,3), 'descend', 'ComparisonMethod', 'abs');
    sort_grouping_diff = sort_grouping_diff(1:9);
    sort_grouping_ind = sort_grouping_ind(1:9);

    % Depict results - overall avg Bliss score
    hFig = figure('units', 'normalized', 'Position', [0 0 1 1], 'Color', 'w');
    % Overall Bliss score
    subplot(3,6, [1:3, 7:8, 13:15]); hold on
    cols = [191 0 0; 0 0 191]/255;    
    for iPhase = 1:length(SET.phases)
        avg = mean(syn.(SET.phases{iPhase}));
        ci = bootci(5000, @mean, syn.(SET.phases{iPhase}));
        xvec = linspace(0, length(avg)/SET.fps, length(avg));
        plot(xvec, avg, 'LineWidth',3, 'Color', cols(iPhase,:))
        plot(xvec, ci(1,:), 'LineWidth',1, 'Color', cols(iPhase,:))
        plot(xvec, ci(2,:), 'LineWidth',1, 'Color', cols(iPhase,:))
    end%iPhase
    xlim([xvec(1), xvec(end)])
    ylim([-0.5 1.5])
    xlabel('time (s)')
    ylabel('Bliss score')

    % Depict results - substantial differences
    plot_pos = [4:6, 10:12, 16:18];
    motif_groups = grouping_unique(sort_grouping_ind,:);
    for iDiff = 1:length(sort_grouping_diff)
        % Calculate synergy
        a = motif.(['m', num2str(motif_groups(iDiff,1))]);
        b = motif.(['m', num2str(motif_groups(iDiff,2))]);
        c = motif.(['m', num2str(motif_groups(iDiff,3))]);
        % Get expected and observed
        f_expected = (a/2 + b/2) - ((a/2).*(b/2));
        f_observed = c;
        % Get BS
        curr_syn = 100*(f_observed - f_expected);
        % Plot
        subplot(3,6, plot_pos(iDiff))
        plot(xvec, curr_syn, 'Color', [191*(sort_grouping_diff(iDiff)>0), 0, 191*(sort_grouping_diff(iDiff)<0)]/255, 'LineWidth', 2)
        xlim([xvec(1), xvec(end)])
        ylim([-4 8])
        yticks([-4 0 8])
        xticks([])
        box on
        title({num2str(motif_groups(iDiff,:)); num2str(round(sort_grouping_diff(iDiff),2))})
    end%iDiff
    % Save
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\motif_combinations_BS'], '-pdf', '-painters')
    close(hFig)
end%iLayer
clearvars -except PoolData SET


%% Visualize the results.
% Proportion of ROIS responding to unique stimuli

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % One figure per layer
    hFig = figure('Name', 'n_unique'); hold on
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % Get color
        col = SET.HeatmapColors.(SET.phases{iPhase}); col = col(floor(size(col,1)/2),:);
        % Get average and CI across animals
        avg = nanmean(bootstrp(5000, @nanmean, PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique));
        ci = bootci(5000, @nanmean, PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique);
        % Iterate over everything
        for iStim = 1:size(PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique,2)
            % Get position of current bar
            pos = [iStim+iPhase/2-1+0.25*(iPhase==1), 0, 0.25, avg(iStim)];
            x = iStim+iPhase/2-1+0.25*(iPhase==1) + 0.25/2;
            % Show bar
            rectangle('position', pos, 'FaceColor', col, 'EdgeColor', 'none')
            % Add CIs
            plot([1 1]*x, ci(:,iStim), 'k')
            plot([-0.25/2 0.25/2]+x, [1 1]*ci(1,iStim), 'k')
            plot([-0.25/2 0.25/2]+x, [1 1]*ci(2,iStim), 'k')
            % Indicate what has been checked
            text(iStim, 0.175, PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique_info{1,iStim}, 'horizontalAlignment', 'center', 'FontSize', 4, 'Color',[0 1 0],"Interpreter","none")
            text(iStim, 0.17,   PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique_info{2,iStim}, 'horizontalAlignment', 'center', 'FontSize', 4, 'Color',[1 0 0],"Interpreter","none")
        end%iStim
    end%iPhase
    % Cosmetics
    xlim([0.5 size(PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique,2)+0.5])
    ylim([0 0.5])
    set(gca, 'xtick',1:size(PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).n_unique,2) )
    ylabel('prop. ROIs')
    title('Proportion of ROIs uniquely responding to given stimulus')
    % Save
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\kmeans_n_unique'], '-pdf', '-painters')
    close(hFig)

    % Statistics. Note that this can take considerable time.
    PoolData.Stats.(SET.layers{iLayer}).seed=1234;
    PoolData.Stats.(SET.layers{iLayer}).n = 5e7;
    PoolData.Stats.(SET.layers{iLayer}).nComp = 1;
    % Preallocation
    PoolData.Stats.(SET.layers{iLayer}).kmeans_n_unique.p = nan(size(PoolData.(SET.phases{1}).(SET.layers{iLayer}).n_unique,2),1);
    PoolData.Stats.(SET.layers{iLayer}).kmeans_n_unique.s = nan(size(PoolData.(SET.phases{1}).(SET.layers{iLayer}).n_unique,2),1);
    PoolData.Stats.(SET.layers{iLayer}).kmeans_n_unique.c = nan(size(PoolData.(SET.phases{1}).(SET.layers{iLayer}).n_unique,2),1);
    % Iterate over all stimuli
    for iStim = 1:size(PoolData.(SET.phases{1}).(SET.layers{iLayer}).n_unique,2)
        [...
            PoolData.Stats.(SET.layers{iLayer}).kmeans_n_unique.p(iStim),...
            PoolData.Stats.(SET.layers{iLayer}).kmeans_n_unique.s(iStim),...
            ~,...
            PoolData.Stats.(SET.layers{iLayer}).kmeans_n_unique.c(iStim)] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample', PoolData.(SET.phases{1}).(SET.layers{iLayer}).n_unique(:,iStim), PoolData.(SET.phases{2}).(SET.layers{iLayer}).n_unique(:,iStim), PoolData.Stats.(SET.layers{iLayer}).n, PoolData.Stats.(SET.layers{iLayer}).seed, PoolData.Stats.(SET.layers{iLayer}).nComp);
    end%iStim

    %Save
    save('ConfocalAnalysis_OdorSpace\OdorSpace_kmeans_PoolData.mat', 'PoolData', 'SET')

end%iLayer
clearvars -except PoolData SET


%% Visualize the results.
% Venn diagram for proportion of ROIS responding to combinations of stimuli

% Give names for the different regions (hard coded)
venn_names = {'Laa only', 'Lct only', 'LaaLct only', 'Laa and Lct', 'Laa and LaaLct', 'Lct and LaaLct', 'all'};

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % One figure per layer
    hFig = figure('Name', 'venn'); hold on
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % Plot diagram in one subplot. Provide info further down (subplot
        % below)
        subplot(2,2,iPhase); hold on
        % Get the average
        avg = mean(PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).VennData);
        % And the corresponding confidence interval
        ci = bootci(5000, @mean, PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).VennData);
        % Make sure everything adds up to one
        avg_norm = avg/sum(avg);
        % Plot the Venn diagram. Here, provide the proportions for every
        % region (only_A, only_B, only_C, A_and_B, A_and_C, B_and_C, A_B_C)
        Confocal_SubFcn.draw_venn_diagram(avg_norm(1), avg_norm(2), avg_norm(3), avg_norm(4), avg_norm(5), avg_norm(6), avg_norm(7),...
            SET.Colors,...               Colors
            {'Laa', 'Lct', 'LaaLct'}...  Labels for the main regions
            );
        % Cosmetics
        axis equal tight off
        title((SET.phases{iPhase}))

        % Provide info on the resulting proportions and confidence
        % intervals in the subplot below the corresponding Venn diagram
        subplot(2,2,iPhase+2);
        % Iterate over all regions. Use the names defined above.
        for iVenn = 1:length(venn_names)
            text(1,iVenn,[venn_names{iVenn}, ': ', num2str(round(avg(iVenn),2)), ' [',num2str(round(ci(1,iVenn)',2)), ',', num2str(round(ci(2,iVenn)',2)),']'])
        end%iVenn
        % Cosmetics
        ylim([0, length(venn_names)+1])
        xlim([0 10])
        axis off
    end%iPhase
    % Save
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\venn_diagrams'], '-pdf', '-painters')
    close(hFig)

    % Do statistics. Note that this can take some time
    PoolData.Stats.(SET.layers{iLayer}).seed=1234;
    PoolData.Stats.(SET.layers{iLayer}).n = 5e7;
    PoolData.Stats.(SET.layers{iLayer}).nComp = 1;
    for iVenn = 1:length(venn_names)
        [...
            PoolData.Stats.(SET.layers{iLayer}).VennDiagram.p(iVenn),...
            PoolData.Stats.(SET.layers{iLayer}).VennDiagram.s(iVenn),...
            ~,...
            PoolData.Stats.(SET.layers{iLayer}).VennDiagram.c(iVenn)] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample', PoolData.(SET.phases{1}).(SET.layers{iLayer}).VennData(:,iVenn), PoolData.(SET.phases{2}).(SET.layers{iLayer}).VennData(:,iVenn), PoolData.Stats.(SET.layers{iLayer}).n, PoolData.Stats.(SET.layers{iLayer}).seed, PoolData.Stats.(SET.layers{iLayer}).nComp);
    end%iVenn

    %Save
    save('ConfocalAnalysis_OdorSpace\OdorSpace_kmeans_PoolData.mat', 'PoolData', 'SET')

end%iLayer
clearvars -except PoolData SET


%% Visualize the results.
% Distance between stimuli based on amplitude (top left) or on active
% regions (bottom right)

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        hFig = figure('Name', 'stim_dist', 'units', 'normalized', 'Position', [0 0 0.5 0.5], 'Color', 'w');

        % Distance between stimuli based on their amplitudes
        subplot(1,2,1); hold on
        % Get the average
        avg1 = nanmean(PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).stim_dist_val,3);
        % Plot heatmap
        imagesc(avg1, "AlphaData", ~isnan(avg1))
        % Cosmetics
        axis equal tight
        colorbar
        title('stim dist val')
        set(gca, 'XTick', 1:length(SET.StimList), 'XTickLabel', SET.StimList, 'TickLabelInterpreter', 'none')
        set(gca, 'YTick', 1:length(SET.StimList), 'YTickLabel', SET.StimList, 'TickLabelInterpreter', 'none')
        colormap(SET.HeatmapColors.(SET.phases{iPhase}))

        % Distance between stimuli based on their activation profile
        subplot(1,2,2); hold on
        % Get the average
        avg2 = nanmean(PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).stim_dist_act,3);
        % Plot heatmap
        imagesc(avg2, "AlphaData", ~isnan(avg2))
        % Cosmetics
        axis equal tight
        colorbar
        title('stim dist act')
        set(gca, 'XTick', 1:length(SET.StimList), 'XTickLabel', SET.StimList, 'TickLabelInterpreter', 'none')
        set(gca, 'YTick', 1:length(SET.StimList), 'YTickLabel', SET.StimList, 'TickLabelInterpreter', 'none')

        % Set the clim
        subplot(1,2,1)
        caxis([floor(min([avg2(:);avg2(:)])*10)/10, ceil(max([avg2(:);avg2(:)])*10)/10])
        subplot(1,2,2)
        caxis([floor(min([avg2(:);avg2(:)])*10)/10, ceil(max([avg2(:);avg2(:)])*10)/10])
        colormap(SET.HeatmapColors.(SET.phases{iPhase}))

        % Save
        export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\kmeans_stim_dist_',SET.phases{iPhase}], '-pdf', '-painters')
        close(hFig)

    end%iPhase
    clearvars -except PoolData SET
end%iLayer


%% Visualize the results.
% Animals per response motif

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % One figure per layer
    hFig = figure('Units','normalized','Position',[0 0 1 1]);
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % Get list of respective cluster IDs for this phase
        id_list = unique( PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).ID(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(:,1) == iPhase));
        % Get list of all possible animals for this phase
        ani_list = unique(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(:,1) == iPhase, 2));
        % Iterate over all clusters
        for iK = 1:length(id_list)
            % Get index positions and observations
            idx = find(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).ID == id_list(iK));
            anis = PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(idx,2);
            % Count how many times each animal was present
            cnts = histcounts(anis, [ani_list(:); max(ani_list)+1]);
            % Normalize by number of ROIs per animal
            cnts = cnts./histcounts(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(:,2), [ani_list(:); max(ani_list)+1]);
            % Get percentage
            cnts = cnts*100;
            % Depict as pie chart
            nexttile
            ax = gca;
            pie(ax, cnts);
            colormap(ax, SET.HeatmapColors.(SET.phases{iPhase}))
            title([SET.phases{iPhase}(1), num2str(iK)])
        end%iC
        % Depict centroids of the cluster
        % First get colors
        cols = SET.HeatmapColors.(SET.phases{iPhase});
        cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), length(id_list)+2)),:);
        cols = cols(2:end-1,:);
        % Then, iterate over all clusters
        for iK = 1:length(id_list)
            % Depict centroid
            nexttile
            plot(nanmean(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).centroid{iPhase}{iK}), 'color', cols(iK,:), 'LineWidth', 2)
            title(['cluster ',num2str(iK)])
            axis tight off
        end%iC
    end%iPhase
    % Save
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\pie_RespProfile'], '-pdf', '-painters')
    close(hFig)
    clearvars -except PoolData SET
end%iLayer


%% Visualize the results.
% Animals per response vector cluster

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % One figure per layer
    hFig = figure('Units','normalized','Position',[0 0 1 1]);
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % Get list of respective cluster IDs
        id_list = unique(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_ID(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(:,1) == iPhase));
        % Get list of all possible animals
        ani_list = unique(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(:,1) == iPhase,2));
        % Iterate over all clusters
        for iK = 1:length(id_list)
            % Get index positions and observations
            idx = find(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_ID == id_list(iK));
            anis = PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(idx,2);
            % Count how many times each animal was present
            cnts = histcounts(anis, [ani_list(:); max(ani_list)+1]);
            % Normalize by number of ROIs per animal
            cnts = cnts./histcounts(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(:,2), [ani_list(:); max(ani_list)+1]);
            % Get percentage
            cnts = cnts*100;
            % Depict as pie chart
            nexttile
            ax = gca;
            pie(ax, cnts);
            colormap(ax, SET.HeatmapColors.(SET.phases{iPhase}))
            title([SET.phases{iPhase}(1), num2str(iK)])
        end%iC
    end%iPhase
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\pie_RespVec'], '-pdf', '-painters')
    close(hFig)
    clearvars -except PoolData SET
end%iLayer


%% Visualize the results.
% Show PCA olfaction space and clustering

% Iterate over layers
for iLayer = 1:length(SET.layers)

    % Create one big figure per layer -------------------------------------
    hFig = figure('Units','normalized','Position',[0 0 1 1]);

    % Bootstrap PCA data to generate p-confidence ellipse around the 
    % resulting distribution
    % --- Keep track of the data's range
    xaxis_lim = []; yaxis_lim = [];
    dist2center_vec = [];
    % --- Iterate over all clusters
    for iK = 1:PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k
        % Pool data from both phases
        curr_dat = [];
        for iPhase = 1:length(SET.phases)
            curr_dat = [curr_dat; PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_pca{iPhase}{iK}];
        end%iPhase
        %Kick out NaNs
        curr_dat = curr_dat(~all(isnan(curr_dat), 2),:);
        % Bootstrap data for the ellipse
        dat_pca_boot = bootstrp(5000, @nanmean, curr_dat);
        % Get a 99%-confidence ellipse around the bootstrapped data
        r_ellipse{iK} = Confocal_SubFcn.dataEllipse(dat_pca_boot, 0.99);
        % Get the range of the data
        xaxis_lim = [xaxis_lim; min(curr_dat(:,1)), max(curr_dat(:,1))];
        yaxis_lim = [yaxis_lim; min(curr_dat(:,2)), max(curr_dat(:,2))];
        dist2center_vec = [dist2center_vec; min(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dist2center_sorted{iPhase}{iK}(:)), max(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dist2center_sorted{iPhase}{iK}(:))];
    end%iK
    % Cosmetics
    xaxis_lim = [min(xaxis_lim(:)), max(xaxis_lim(:))]*1.15;
    yaxis_lim = [min(yaxis_lim(:)), max(yaxis_lim(:))]*1.15;
    dist2center_vec = [min(dist2center_vec(:)), max(dist2center_vec(:))];

    % Raw PCA space (both phases combined) --------------------------------
    nexttile; hold on
    % Get nice colors
    cols = Confocal_SubFcn.ColMapInferno(size(SET.HeatmapColors.(SET.phases{iPhase}),1));
    cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1),sum(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k)+2)),:);
    cols = cols(2:end-1,:);
    % Iterate over all clusters
    for iK = 1:PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k
        % Pool data from both phases
        curr_dat = [];
        for iPhase = 1:length(SET.phases)
            curr_dat = [curr_dat; PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_pca{iPhase}{iK}];
        end%iPhase
        % Plot elipse
        plot(r_ellipse{iK}(:,1), r_ellipse{iK}(:,2), 'Color', cols(iK,:), 'LineWidth', 1)
        % Plot scatter
        plot(curr_dat(:,1), curr_dat(:,2), 'o', 'MarkerFaceColor', cols(iK,:), 'MarkerEdgeColor', 'none')
    end%iK
    % Cosmetics
    ax = gca;
    colormap(ax, cols)
    colorbar
    clim([1, size(cols,1)])
    title('cluster')
    axis equal
    xlim(xaxis_lim)
    ylim(yaxis_lim)
    xticks([])
    yticks([])
    box on

    % Raw PCA space (phases color-coded) ----------------------------------
    nexttile; hold on
    % Get nice colors
    cols = Confocal_SubFcn.ColMapInferno(size(SET.HeatmapColors.(SET.phases{iPhase}),1));
    cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1),sum(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k)+2)),:);
    cols = cols(2:end-1,:);
    % Iterate over all clusters
    for iK = 1:PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k
        % Iterate over both phases and get data
        curr_dat = []; curr_phase = []; col = [];
        for iPhase = 1:length(SET.phases)
            % Get data
            curr_dat = [curr_dat; PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_pca{iPhase}{iK}];
            curr_phase = [curr_phase; ones(size(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_pca{iPhase}{iK},1),1)*iPhase];
            % Get nice colors
            col = [col;...
                ones(size(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_pca{iPhase}{iK},1),1)*0.75*(iPhase==1),...
                zeros(size(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_pca{iPhase}{iK},1),1),...
                ones(size(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_pca{iPhase}{iK},1),1)*0.75*(iPhase==2)];
        end%iPhase
        % Plot all points and connect them to the controid
        C = nanmean(curr_dat,1);
        plot(C(1), C(2), 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none')
        for iP = 1:size(curr_dat,1)
            % Plot line to controid
            plot([C(1), curr_dat(iP,1)], [C(2), curr_dat(iP,2)], 'k')
            % Plot scatter
            plot(curr_dat(iP,1), curr_dat(iP,2), 'o', 'MarkerFaceColor', col(iP,:), 'MarkerEdgeColor', 'none')
        end%iPhase
        % Plot elipse
        plot(r_ellipse{iK}(:,1), r_ellipse{iK}(:,2), 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        % Plot convhull
        curr_dat = curr_dat(~all(isnan(curr_dat), 2),:);
        k = convhull(curr_dat(:,1),curr_dat(:,2));
        plot(curr_dat(k,1),curr_dat(k,2), 'color', cols(iK,:))
    end%iK
    % Cosmetics
    ax = gca;
    title('cluster')
    axis equal
    xlim(xaxis_lim)
    ylim(yaxis_lim)
    xticks([])
    yticks([])
    box on

    % Raw PCA space (dist2center color-coded) -----------------------------
    nexttile; hold on
    % Get colors and get a matching distance vector to assign colors correctly
    cols = Confocal_SubFcn.ColMapInferno(size(SET.HeatmapColors.(SET.phases{iPhase}),1));
    dist2center_vec = linspace(min(dist2center_vec), max(dist2center_vec), size(cols,1));
    % Iterate over all cluster
    for iK = 1:PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k
        % Pool phases
        curr_dat = [];
        curr_dist = [];
        for iPhase = 1:length(SET.phases)
            curr_dat = [curr_dat; PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_pca{iPhase}{iK}];
            curr_dist = [curr_dist; PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dist2center_sorted{iPhase}{iK}];
        end%iPhase
        % Plot elipse
        [~,d_index] = min(abs(dist2center_vec - nanmean(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dist2center_sorted{iPhase}{iK})));
        plot(r_ellipse{iK}(:,1), r_ellipse{iK}(:,2), 'Color', cols(d_index,:), 'LineWidth', 1)
        % Plot scatter
        for iAni = 1:size(curr_dat,1)
            if ~isnan(curr_dat(iAni,1))
                [~,d_index] = min(abs(dist2center_vec - curr_dist(iAni)));
                plot(curr_dat(iAni,1), curr_dat(iAni,2), 'o', 'MarkerFaceColor', cols(d_index,:), 'MarkerEdgeColor', 'none')
            end%if
        end%iAni
    end%iK
    % Cosmetics
    ax = gca;
    colormap(ax, cols)
    colorbar
    clim([dist2center_vec(1), dist2center_vec(end)])
    title('dist2center')
    axis equal
    xlim(xaxis_lim)
    ylim(yaxis_lim)
    xticks([])
    yticks([])
    box on

    % Violin plots (dist2center) ------------------------------------------
    nexttile; hold on
    % Pool data across both phases
    curr_dat = [];
    for iPhase = 1:length(SET.phases)
        curr_dat = [curr_dat; cell2mat(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dist2center_sorted{iPhase})];
    end%iPhase
    % Do stats
    % --- ANOVA
    [p,~,aov] = anova1(curr_dat, 1:PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k,'off');
    % --- Multiple comparison
    comp_results = multcompare(aov, 'CriticalValueType', 'hsd', 'Display', 'off');
    comp_results = sortrows(comp_results);
    comp_results = array2table(comp_results, 'VariableNames', {'Group1', 'Group2', 'MeanDifferenceLower', 'MeanDifference', 'MeanDifferenceUpper', 'pValue'});
    % --- Extract the comparison group indices and significant differences
    CLD = Confocal_SubFcn.multcompLetters(comp_results);
    PoolData.Stats.(SET.layers{iLayer}).(SET.phases{iPhase}).dist2center_vector.anova1_p = p;
    PoolData.Stats.(SET.layers{iLayer}).(SET.phases{iPhase}).dist2center_vector.comp_results = comp_results;
    PoolData.Stats.(SET.layers{iLayer}).(SET.phases{iPhase}).dist2center_vector.CLD = CLD;

    % Properties for the violin plot
    properties.MinVal =             [];                 % Smallest possible value (e.g. errors = 0)
    properties.MaxVal =             [];                 % Biggest possible value
    properties.AvgType =            'mean';             %(Set which measure should be plotted: 'median', 'mean' or 'both')
    properties.MeanWidth =          2;                  %(Set width of the median's line and/or symbol)
    properties.SeparateOutliers =   0;                  %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
    properties.MeanCol =            [0.251 0.251 0.251];
    % Iterate over all clusters
    for iK = 1:PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k
        Confocal_SubFcn.violinplot_advanced(curr_dat(:,iK), iK, 0.75, properties)
        text(iK, 0.1, CLD{iK}, 'color', [1 0 0], 'HorizontalAlignment', 'center')
    end%iK
    % Get the distance of an uniformly distributed sample
    dat = rand(1e6,2)*2;
    dat = dat-(mean(dat));
    % Get each point's distance
    dat = sqrt(sum(dat'.*dat'))';
    dat(dat>1) = [];
    plot([0, PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k+1], [1 1]*mean(dat))
    clear dat
    % Cosmetics
    title(['avg dist2center | ',SET.phases{iPhase}])
    set(gca, 'XTick', 1:PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k)
    ylabel('norm. distance to center')
    xlabel('response vector ID')
    ylim([0 1])
    xlim([0, PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k+1])

    % Export & Save -----------------------------------------------------------
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\kmeans_PCA_space'], '-pdf', '-painters')
    save('ConfocalAnalysis_OdorSpace\OdorSpace_kmeans_PoolData.mat', 'PoolData', 'SET')
    close(hFig)

end%iLayer


%% Visualize the results.
% Show how the clusters are different from each other (response vectors)

% Set the ylim (zscored)
ylim_range = [-1.5 1.5];

% Iterate over layers
for iLayer = 1:length(SET.layers)

    % Prepare a figure for the response vectors
    hFig_vec = figure('Units','normalized','Position',[0 0 0.5 1]);
    % Prepare a figure for the allocation of vectors per phase
    hFig_pi = figure('Units','normalized','Position',[0 0 0.5 1]);
    % Prepare a figure to show the location of all (sub)cluster
    hFig_cluster = figure('Units','normalized','Position',[0 0 0.5 1]);

    % Create some variables
    PCA_centroid = [];
    PCA_centroid_pca = [];
    PCA_label = {}; cnt = 1;

    % Iterate over all clusters
    for iK = 1:PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k

        % Pool data across phases
        resp_vec = [];
        temp_PCA_centroid = [];
        temp_PCA_centroid_pca = [];
        for iPhase = 1:length(SET.phases)
            % resp_vec = [resp_vec; PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid{iPhase}{iK}];
            resp_vec =              [resp_vec;              PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_norm{iPhase}{iK}];
            temp_PCA_centroid =     [temp_PCA_centroid;     nanmean(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid{iPhase}{iK})];
            temp_PCA_centroid_pca = [temp_PCA_centroid_pca; nanmean(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_centroid_pca{iPhase}{iK})];
        end%iPhase
        resp_vec = resp_vec(~all(isnan(resp_vec), 2),:);
        PCA_centroid = [PCA_centroid; nanmean(temp_PCA_centroid,1)];
        PCA_centroid_pca = [PCA_centroid_pca; nanmean(temp_PCA_centroid_pca,1)];
        PCA_label{cnt} = ['cluster ',num2str(iK)];

        % Do stats
        % --- get average response
        avg = nanmean(resp_vec,1);
        ci = bootci(5000, @nanmean, resp_vec);
        % --- ANOVA
        aov = anova(resp_vec);
        % --- Multiple comparison
        comp_results = multcompare(aov, CriticalValueType="hsd");
        comp_results = sortrows(comp_results);
        % --- Extract the comparison group indices and significant differences
        CLD = Confocal_SubFcn.multcompLetters(comp_results);

        % Get colors
        cols = SET.HeatmapColors.(SET.phases{iPhase});
        col_vec = linspace(min(avg), max(avg), size(cols,1));

        % Show bar plot
        figure(hFig_vec)
        nexttile; hold on
        % Iterate over all stimuli
        for iStim = 1:length(avg)
            % Check whether the bar is pointing in pos. or neg. direction
            if avg(iStim) >= 0
                pos = [iStim-0.5 0 1 avg(iStim)];
            else
                pos = [iStim-0.5 avg(iStim) 1 abs(avg(iStim))];
            end
            % Check whether zero lies withing the CI
            if sum(ci(:,iStim)>0)==2
                rectangle('position', pos, 'EdgeColor', 'none', 'FaceColor', [1 0 1])
            elseif sum(ci(:,iStim)<0)==2
                rectangle('position', pos, 'EdgeColor', 'none', 'FaceColor', [0 1 1])
            else
                rectangle('position', pos, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5])
            end
            % Plot CI
            plot([iStim iStim],ci(:,iStim),'k')
            % Add info on stats
            text(iStim, avg(iStim), CLD{iStim}, 'HorizontalAlignment', 'center', 'Color', [1 0 0])
        end%iStim
        % Cosmetics
        xlim([0.5, length(avg)+0.55])
        ylim(ylim_range)
        title(PCA_label{cnt})
        box on
        set(gca, 'XTick', 1:length(SET.StimList), 'XTickLabel', SET.StimList, 'TickLabelInterpreter', 'none')

        % Show the allocation of vectors per phase
        figure(hFig_pi)
        % Iterate over both phases
        for iPhase = 1:length(SET.phases)
            % Get data and colors
            col = SET.HeatmapColors.(SET.phases{iPhase}); col = col(floor(size(col,1)/2),:);
            col = [0.251, 0.251, 0.251; col];
            col = [0.75, 0.75, 0.75; 0.75*(iPhase==1) 0 0.75*(iPhase==2)];
            % Get the proportion of ROIs of the current phase that were put
            % into the current cluster
            prop = length(find(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_ID == iK & PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(:,1) == iPhase)) / sum(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(:,1) == iPhase);
            prop = [1- prop; prop];
            % Plot pie chart
            nexttile; hold on
            ax = gca;
            pie(ax, prop);
            colormap(ax, col)
            axis equal tight off
            title([SET.phases{iPhase}(1), num2str(iK)])
        end%iPhase

        % Show cluster's location on map
        figure(hFig_cluster); hold on
        plot(r_ellipse{iK}(:,1), r_ellipse{iK}(:,2), 'k', 'LineWidth', 1)
        text(PCA_centroid_pca(cnt,1)+0.05, PCA_centroid_pca(cnt,2), ['c',num2str(iK)], "HorizontalAlignment", "left")
        cnt=cnt+1;

    end%iC

    % Save figures
    % --- allocation
    figure(hFig_pi)
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\resp_vec_porp_pie'], '-pdf', '-painters')
    close(hFig_pi)
    % --- cluster location
    figure(hFig_cluster)
    axis equal
    xlim(xaxis_lim)
    ylim(yaxis_lim)
    xticks([])
    yticks([])
    box on
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\kmeans_cluster_locations'], '-pdf', '-painters')
    close(hFig_cluster)
    % --- response vectors
    figure(hFig_vec)
    nexttile
    D = squareform(pdist(PCA_centroid, SET.dist_pca));
    tree = linkage(D,'single', SET.dist_pca);
    dendrogram(tree, 'Labels', PCA_label)
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\kmeans_cluster_responses'], '-pdf', '-painters')
    close(hFig_vec)

end %iLayer

% Clean
clearvars -except PoolData SET


%% Visualize the results.
% Show how response motif assignments change between stimuli

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % Get avg values for each ROI (row) and each stimulus (column)
        data = PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dat_AvgVal(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(:,1) == iPhase, :);
        % Get corresponding grouping factors (response motifs)
        motifs = PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).ID(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(:,1) == iPhase,:);
        stims = PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(:,1) == iPhase,3);
        data_grouping = [];
        % Iterate ove all stimuli and populate the grouping variable
        for iStim = 1:length(SET.StimList)
            idx = find(stims == iStim);
            data_grouping = [data_grouping, motifs(idx)];
        end
        % Plot response hive (circular stacked bar plot)
        hFig = figure;
        Confocal_SubFcn.HivePlot(data(:,1:3), data_grouping(:,1:3), [0 0 0], [], 0.5, length(unique(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).ID(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(:,1) == iPhase,:)))*0.5)
        % Add a colorbar for the different levels
        cols = SET.HeatmapColors.(SET.phases{iPhase});
        cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).best_k+2)),:);
        cols = cols(2:end-1,:);
        colormap(cols); colorbar; clim([1 PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).best_k])
        % Add labels
        set(gca,'ThetaTickLabel', SET.StimList, 'TickLabelInterpreter', 'none')
        % Save
        export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\response_hive_',SET.phases{iPhase}], '-pdf', '-painters')
        close(hFig)
        
        % Iterate over all stimulus and motif-change combinations to get an
        % overview for the dynamics of motif transitions.
        % --- Get list of animals
        ani_list = unique(PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info(:,1));
        % --- Get list of motifs
        motif_list = unique(data_grouping); 
        % --- Get all possible combinations of motifs
        motif_change_pairs = repmat(motif_list, 1, 2);
        motif_change_list = perms(motif_list);  motif_change_list = motif_change_list(:,1:2);
        motif_change_list = [motif_change_list; motif_change_pairs];
        motif_change_list = fliplr(unique(motif_change_list, 'rows'));
        clear motif_change_pairs
        % --- Get all possible combinations of stimuli
        stim_change_list = perms(1:length(SET.StimList));  stim_change_list = stim_change_list(:,1:2);
        stim_change_list = fliplr(unique(stim_change_list, 'rows'));
        clear stim_change_pairs
         
        % % Kick motif 1 as it contains noise only
        % motif_change_list(find(sum(motif_change_list==1, 2)>0),:) = [];

        % Get data
        % Preallocation
        motif_dynamics = nan(size(stim_change_list,1), size(motif_change_list,1), length(ani_list));
        motif_dynamics_features = nan(length(ani_list), size(stim_change_list,1)*size(motif_change_list,1));
        stim_change_list_labels = cell(1,size(stim_change_list,1));
        motif_change_list_labels = cell(1,size(motif_change_list,1));
        motif_dynamics_features_labels = cell(1, size(stim_change_list,1)*size(motif_change_list,1));
        cnt = 1;
        % Iterate over all motif changes
        for iMotif = 1:size(motif_change_list,1)
            % Iterate over all stimulus changes
            for iStim = 1:size(stim_change_list,1)
                % Iterate over all animals
                for iAni = 1:length(ani_list)
                    % Get how many granules made the current motif
                    % transition at the given stimulus change
                    motif_dynamics(iStim, iMotif, iAni) = sum(...
                        data_grouping(:,stim_change_list(iStim,1))==motif_change_list(iMotif,1) & ...
                        data_grouping(:,stim_change_list(iStim,2))==motif_change_list(iMotif,2) & ...
                        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info(:,1) == ani_list(iAni));
                    % Normalize by the number of origin granules
                    origin_cnt = sum(...
                        data_grouping(:,stim_change_list(iStim,1))==motif_change_list(iMotif,1) & ...
                        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info(:,1) == ani_list(iAni));
                    if origin_cnt > 0
                        motif_dynamics(iStim, iMotif, iAni) = ...
                            motif_dynamics(iStim, iMotif, iAni) / origin_cnt;
                    else
                        motif_dynamics(iStim, iMotif, iAni) = NaN;
                    end% if motif transition was possible in the first place
                    % Also pool it as features for each animal
                    motif_dynamics_features(iAni, cnt) = motif_dynamics(iStim, iMotif, iAni);
                    motif_dynamics_features_labels{cnt} = [...
                        'S',num2str(stim_change_list(iStim,1)),'>S',num2str(stim_change_list(iStim,2)),...
                        '|',...
                        'M',num2str(motif_change_list(iMotif,1)),'>M',num2str(motif_change_list(iMotif,2))];                    
                end%iAni
                cnt = cnt+1;
                stim_change_list_labels{iStim} = [num2str(stim_change_list(iStim,1)),'>',num2str(stim_change_list(iStim,2))];
            end%iStim
            motif_change_list_labels{iMotif} = [num2str(motif_change_list(iMotif,1)),'>',num2str(motif_change_list(iMotif,2))];
        end%iMotif
        % Store
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).stim_change_list = stim_change_list;
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).motif_change_list = motif_change_list;
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).motif_dynamics = motif_dynamics;
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).motif_dynamics_features = motif_dynamics_features;
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).motif_change_list_labels = motif_change_list_labels;
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).stim_change_list_labels = stim_change_list_labels;
        PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).motif_dynamics_features_labels = motif_dynamics_features_labels;
    end%iPhase

    % Compare gregarious and solitarious motif dynamics
    PoolData.Stats.motif_dynamics.seed=1234;
    PoolData.Stats.motif_dynamics.n = 5e3;
    PoolData.Stats.motif_dynamics.nComp = 1;
    PoolData.Stats.motif_dynamics.pValues = nan(length(stim_change_list), length(motif_change_list));
    for iMotif = 1:size(motif_change_list,1)
        % Iterate over all stimulus changes
        for iStim = 1:size(stim_change_list,1)
            greg = squeeze(PoolData.gregarious.(SET.layers{iLayer}).motif_dynamics(iStim, iMotif, :));
            soli = squeeze(PoolData.solitarious.(SET.layers{iLayer}).motif_dynamics(iStim, iMotif, :));
            greg = greg(~isnan(greg)); soli = soli(~isnan(soli));
            PoolData.Stats.(SET.phases{iPhase}).motif_dynamics.pValues(iStim, iMotif) = ranksum(greg, soli);
        end%iStim
    end%iMotif

    % Show heat maps for both phases and the resulting difference
    hFig = figure('color', 'w', 'units', 'normalized', 'Position', [0 0 1 1 ]);
    % Get data
    img_greg = nanmean(PoolData.gregarious.(SET.layers{iLayer}).motif_dynamics, 3);
    img_soli = nanmean(PoolData.solitarious.(SET.layers{iLayer}).motif_dynamics, 3);
    img_diff = img_greg-img_soli;
    % Get colormap for differences
    cols = [...
        interp1([min(img_diff(:)) 0 max(img_diff(:))], [0, 1, 0.75], linspace(min(img_diff(:)),max(img_diff(:)),1000))',...
        interp1([min(img_diff(:)) 0 max(img_diff(:))], [0, 1, 0], linspace(min(img_diff(:)),max(img_diff(:)),1000))',...
        interp1([min(img_diff(:)) 0 max(img_diff(:))], [0.75, 1, 0], linspace(min(img_diff(:)),max(img_diff(:)),1000))'];
    % --- Gregarious
    ax = subplot(4,1,1);
    imagesc(img_greg)
    axis equal tight
    colormap(ax, SET.HeatmapColors.gregarious)
    clim([0, round(max(img_greg(:)),2)])
    cb = colorbar;
    cb.Ticks = [0, round(max(img_greg(:)),2)];
    yticks(1:length(stim_change_list));
    xticks(1:length(motif_change_list));
    set(ax, 'XTickLabels', motif_change_list_labels, 'YTickLabels', stim_change_list_labels)
    xtickangle(90)
    title('gregarious')
    % --- Solitarious
    ax = subplot(4,1,2);
    imagesc(img_soli)
    axis equal tight
    colormap(ax, SET.HeatmapColors.solitarious)
    clim([0, round(max(img_soli(:)),2)])
    cb = colorbar;
    cb.Ticks = [0, round(max(img_soli(:)),2)];
    yticks(1:length(stim_change_list));
    xticks(1:length(motif_change_list));
    set(ax, 'XTickLabels', motif_change_list_labels, 'YTickLabels', stim_change_list_labels)
    xtickangle(90)
    title('solitarious')
    % --- Difference
    ax = subplot(4,1,3);
    imagesc(img_diff)
    axis equal tight
    colormap(ax, cols)
    clim([round(min(img_diff(:)),2), round(max(img_diff(:)),2)])
    cb = colorbar;
    cb.Ticks = [round(min(img_diff(:)),2), 0, round(max(img_diff(:)),2)];
    yticks(1:length(stim_change_list));
    xticks(1:length(motif_change_list));
    set(ax, 'XTickLabels', motif_change_list_labels, 'YTickLabels', stim_change_list_labels)
    xtickangle(90)
    title('greg-soli')
    % ---- Significance
    ax = subplot(4,1,4);
    imagesc(PoolData.Stats.(SET.phases{iPhase}).motif_dynamics.pValues<0.05)
    axis equal tight
    colormap(ax, gray(2))
    cb = colorbar;
    yticks(1:length(stim_change_list));
    xticks(1:length(motif_change_list));
    set(ax, 'XTickLabels', motif_change_list_labels, 'YTickLabels', stim_change_list_labels)
    xtickangle(90)    
    title('p-value (ranksum)')
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\motif_transition_dynamics'], '-pdf', '-painters')
    close(hFig)


    % Try to predict the phase of the animal based on its motif
    % transitions
    all_motifs = [...
        PoolData.gregarious.(SET.layers{iLayer}).motif_dynamics_features;...
        PoolData.solitarious.(SET.layers{iLayer}).motif_dynamics_features];
    all_labels = [...
        zeros(size(PoolData.gregarious.(SET.layers{iLayer}).motif_dynamics_features,1),1);...
        ones(size(PoolData.solitarious.(SET.layers{iLayer}).motif_dynamics_features,1),1)];
    label_list = unique(all_labels);    
    % Iterate over all animals, leave one out, train an LDA and predict the
    % phase of the excluded animal
    model_classification = nan(length(all_labels),1);
    featureImportance = zeros(length(all_labels), size(all_motifs,2));
    for iRep = 1:length(all_labels)
        % Get the indices for the training data
        idx_train = 1:length(all_labels); idx_train(idx_train==iRep)=[];
        % Get the index for the test data
        idx_test = iRep;
        % Get weights
        w = all_labels(idx_train,:);
        idx_0 = find(w==0);
        idx_1 = find(w==1);
        w(idx_0) = 1/length(idx_0);
        w(idx_1) = 1/length(idx_1);
        w = w/sum(w);
        c = [0, length(idx_1); length(idx_0),0];
        c = c/min(c(c>0));
        % Train the model
        ldaModel = fitcdiscr(all_motifs(idx_train,:),all_labels(idx_train,:), 'DiscrimType', 'diaglinear', 'Cost', c, 'Weights', w);
        % Predict label for held-out data
        model_classification(iRep) =  predict(ldaModel, all_motifs(idx_test,:));
        % Extract the Coefficients for Each Class Pair
        coeffs = ldaModel.Coeffs;
        % Calculate the Aggregate Importance of Each Feature
        for iStim = 1:length(label_list) % Loop over classes
            for jStim = iStim+1:length(label_list)
                featureImportance(iRep,:) = featureImportance(iRep,:) + abs(coeffs(iStim,jStim).Linear');
            end%iStim
        end%iStim
    end%iRep

    
    % Get the performance based on the AIC
    AICc = @(n, SS, k) n*log(SS/n) + 2*k + ((2*k*(k+1))/(n-k-1));
    model_performance = AICc(length(all_labels), sum((all_labels-model_classification).^2), size(all_motifs,2));
    % Get a list of all features. In the next step we will gradually
    % decrease the feature space until the performance decreases
    good_features = 1:length(motif_dynamics_features_labels);
    % Sort the features by their importance
    [~,feature_idx] = sort(mean(featureImportance), 'descend');
    % Gradually decrease feature space
    for iF = 1:length(good_features)-1
        % Get current set of features
        curr_features = feature_idx(1:end-iF);       
        % Iterate over all animals, leave one out, train an LDA and predict the
        % phase of the excluded animal
        model_classification = nan(length(all_labels),1);
        for iRep = 1:length(all_labels)
            % Get the indices for the training data
            idx_train = 1:length(all_labels); idx_train(idx_train==iRep)=[];
            % Get the index for the test data
            idx_test = iRep;
            % Get weights
            w = all_labels(idx_train,:);
            idx_0 = find(w==0);
            idx_1 = find(w==1);
            w(idx_0) = 1/length(idx_0);
            w(idx_1) = 1/length(idx_1);
            w = w/sum(w);
            c = [0, length(idx_1); length(idx_0),0];
            c = c/min(c(c>0));
            % Train the model
            ldaModel = fitcdiscr(all_motifs(idx_train, curr_features),all_labels(idx_train,:), 'DiscrimType', 'diaglinear', 'Cost', c, 'Weights', w);
            % Predict label for held-out data
            model_classification(iRep) =  predict(ldaModel, all_motifs(idx_test, curr_features));            
        end%iRep
        % Check model
        curr_performance = AICc(length(all_labels), sum((all_labels-model_classification).^2), size(all_motifs,2));
        if curr_performance>model_performance
            best_features = feature_idx(1:end-iF+1);
            break
        elseif curr_performance<model_performance
            model_performance = curr_performance;
        end% if decrease in performance
    end%iF

    % Good, now we know which features are crucial for the discrimination
    % of solitarious and gregarious animals
    % Run the model again with the most important features
    % Iterate over all animals, leave one out, train an LDA and predict the
    % phase of the excluded animal
    model_classification = nan(length(all_labels),1);
    featureImportance = zeros(length(all_labels), length(best_features));
    for iRep = 1:length(all_labels)
        % Get the indices for the training data
        idx_train = 1:length(all_labels); idx_train(idx_train==iRep)=[];
        % Get the index for the test data
        idx_test = iRep;
        % Get weights
        w = all_labels(idx_train,:);
        idx_0 = find(w==0);
        idx_1 = find(w==1);
        w(idx_0) = 1/length(idx_0);
        w(idx_1) = 1/length(idx_1);
        w = w/sum(w);
        c = [0, length(idx_1); length(idx_0),0];
        c = c/min(c(c>0));
        % Train the model
        ldaModel = fitcdiscr(all_motifs(idx_train,best_features),all_labels(idx_train,:), 'DiscrimType', 'diaglinear', 'Cost', c, 'Weights', w);
        % Predict label for held-out data
        model_classification(iRep) =  predict(ldaModel, all_motifs(idx_test,best_features));
        % Extract the Coefficients for Each Class Pair
        coeffs = ldaModel.Coeffs;
        % Calculate the Aggregate Importance of Each Feature
        for iStim = 1:length(label_list) % Loop over classes
            for jStim = iStim+1:length(label_list)
                featureImportance(iRep,:) = featureImportance(iRep,:) + coeffs(iStim,jStim).Linear';
            end%iStim
        end%iStim
    end%iRep   
    
    % Show the  most important feature for each phase for the change to the mixture
    cropped_motif_dynamics_features_labels = motif_dynamics_features_labels(best_features);
    avg_featureImportance = nanmean(featureImportance);
    ci_featureImportance = bootci(5000, {@nanmean, featureImportance});
    % Plot everything
    hFig = figure('color', 'w', 'units', 'normalized', 'Position',[0 0 1 1]);  
    cnt = [1 1 1];
    labels1 = {};
    labels2 = {};
    labels3 = {};
    for iF = 1:length(best_features)
        % Differentiate between target stimuli
        if contains(cropped_motif_dynamics_features_labels{iF}, '>S1')
            col = SET.Colors(1,:);
            subplot(1,3,1); hold on
            % --- Avg as bar
            if avg_featureImportance(iF) > 0
                rectangle('position', [cnt(1)-0.75/2 0, 0.75, avg_featureImportance(iF)], 'FaceColor', col, 'EdgeColor', 'none')
            else
                rectangle('position', [cnt(1)-0.75/2 avg_featureImportance(iF), 0.75, abs(avg_featureImportance(iF))], 'FaceColor', col, 'EdgeColor', 'none')
            end
            % --- Ci as line
            plot([1 1]*cnt(1), ci_featureImportance(:,iF), 'color', [0.5 0.5 0.5])
            labels1{cnt(1)} = cropped_motif_dynamics_features_labels{iF};
            cnt(1) = cnt(1)+1;
        elseif contains(cropped_motif_dynamics_features_labels{iF}, '>S2')
            col = SET.Colors(2,:);
            subplot(1,3,2); hold on
            % --- Avg as bar
            if avg_featureImportance(iF) > 0
                rectangle('position', [cnt(2)-0.75/2 0, 0.75, avg_featureImportance(iF)], 'FaceColor', col, 'EdgeColor', 'none')
            else
                rectangle('position', [cnt(2)-0.75/2 avg_featureImportance(iF), 0.75, abs(avg_featureImportance(iF))], 'FaceColor', col, 'EdgeColor', 'none')
            end
            % --- Ci as line
            plot([1 1]*cnt(2), ci_featureImportance(:,iF), 'color', [0.5 0.5 0.5])
            labels2{cnt(2)} = cropped_motif_dynamics_features_labels{iF};
            cnt(2) = cnt(2)+1;
        elseif contains(cropped_motif_dynamics_features_labels{iF}, '>S3')
            col = SET.Colors(3,:);
            subplot(1,3,3); hold on
            % --- Avg as bar
            if avg_featureImportance(iF) > 0
                rectangle('position', [cnt(3)-0.75/2 0, 0.75, avg_featureImportance(iF)], 'FaceColor', col, 'EdgeColor', 'none')
            else
                rectangle('position', [cnt(3)-0.75/2 avg_featureImportance(iF), 0.75, abs(avg_featureImportance(iF))], 'FaceColor', col, 'EdgeColor', 'none')
            end
            % --- Ci as line
            plot([1 1]*cnt(3), ci_featureImportance(:,iF), 'color', [0.5 0.5 0.5])
            labels3{cnt(3)} = cropped_motif_dynamics_features_labels{iF};
            cnt(3) = cnt(3)+1;
        end%color      
    end%iF
    % Cosmetics    
    % Transition to Laa
    subplot(1,3,1)
    xlim([0, max([length(labels1), length(labels2), length(labels3)])+1])
    ylim([-25 25])
    plot([0, max([length(labels1), length(labels2), length(labels3)])+1], [0 0], 'k')
    ylabel('feature importance')
    set(gca, 'XTick', 1:length(labels1), 'XTickLabel', labels1)  
    view([90 90])   
    % Transition to Lct
    subplot(1,3,2)
    xlim([0, max([length(labels1), length(labels2), length(labels3)])+1])
    ylim([-25 25])
    plot([0, max([length(labels1), length(labels2), length(labels3)])+1], [0 0], 'k')
    ylabel('feature importance')
    set(gca, 'XTick', 1:length(labels2), 'XTickLabel', labels2)  
    view([90 90])
    % Transition to LaaLct
    subplot(1,3,3)
    xlim([0, max([length(labels1), length(labels2), length(labels3)])+1])
    ylim([-25 25])
    plot([0, max([length(labels1), length(labels2), length(labels3)])+1], [0 0], 'k')
    ylabel('feature importance')
    set(gca, 'XTick', 1:length(labels3), 'XTickLabel', labels3)  
    view([90 90])
    % Save
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\motif_transition_features'], '-pdf', '-painters')
    close(hFig)

    hFig = figure;
    plotconfusion(all_labels',model_classification')
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\motif_transition_prediction'], '-pdf', '-painters')
    close(hFig)

end%iLayer
save('ConfocalAnalysis_OdorSpace\OdorSpace_kmeans_PoolData.mat', 'PoolData', 'SET')
clearvars -except PoolData SET


%% Visualize the results.
% Show how how the proportion of response motifs changes between stimuli

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % Get avg values for each ROI (row) and each stimulus (column)
        data = PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).dat_AvgVal(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(:,1) == iPhase, :);
        % Get each ROI's response motif ID for each stimulus and the
        % corresponding grouping factors (response motifs)
        motifs = PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).ID(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(:,1) == iPhase,:);
        stims = PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).info(:,1) == iPhase,3);
        data_grouping = [];
        % Iterate ove all stimuli and populate the grouping variable
        for iStim = 1:length(SET.StimList)
            idx = find(stims == iStim);
            data_grouping = [data_grouping, motifs(idx)];
        end
        % Get info on animals
        info = PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info(:,1);
        
        % Get list of motifs
        id_list = unique(data_grouping(:));
        % Get the average proportion of motif IDs per stimulus across all 
        % animals
        prop_data{iPhase} = zeros(PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).best_k, size(data,2));
        % Get list of animals
        ani_list = unique(PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info(:,1));
        % Iterate over all stimui
        for iStim = 1:length(SET.StimList)
            % Iterate over all motifs
            for iK = 1:length(id_list)
                % Iterate over all animals and get the proportion of
                % motifs for the current stimulus
                for iAni = 1:length(ani_list)
                    p = sum(data_grouping(:,iStim) == id_list(iK) & info== ani_list(iAni));
                    p = p / sum(info == ani_list(iAni));
                    prop_data{iPhase}(iK, iStim) = prop_data{iPhase}(iK, iStim) + p;
                end%iAni
                prop_data{iPhase}(iK, iStim) = prop_data{iPhase}(iK, iStim) / length(ani_list);
            end%iK
        end%iStim
        % Get nice colormap
        cols = SET.HeatmapColors.(SET.phases{iPhase});
        cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), length(id_list)+3)),:);
        cols = cols(2:end-1,:);
        % Plot response hive
        hFig = figure;
        Confocal_SubFcn.HivePropPlot(prop_data{iPhase}(:,1:3), 'connected', cols(2:end,:), [], 0.5, 0.05)
        % Add a colorbar for the different levels
        colormap(cols(2:end,:)); colorbar; clim([1 PoolData.kmeans_result.ResponseProfile.(SET.layers{iLayer}).best_k])
        % Add labels
        set(gca, 'ThetaTickLabel', SET.StimList, 'TickLabelInterpreter', 'none')
        % Save
        export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\response_profile_prop_hive_',SET.phases{iPhase}], '-pdf', '-painters')
        close(hFig)
    end%iPhase

    % Show the difference between the two phases
    cols = [...
        interp1([1 0 -1], [0, 1, 0.75], linspace(1,-1,1000))',...
        interp1([1 0 -1], [0, 1, 0], linspace(1,-1,1000))',...
        interp1([1 0 -1], [0.75, 1, 0], linspace(1,-1,1000))'];
    currData = prop_data{1} - prop_data{2};
    % Plot response hive
    hFig = figure;
    hold on
    Confocal_SubFcn.HiveCircPropPlot(currData(:,1:3), cols, [], 0.5, 0.25)
    % Save
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\response_profile_prop_hive_diff'], '-pdf', '-painters')
    close(hFig)

end%iLayer
clearvars -except PoolData SET


%% Visualize the results.
% Show the response vector locations within the AL

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % One figure per phase
        hFig = figure('units','normalized','Position',[0 0 1 1], 'color', 'w');
        % Iterte over all animals
        ani_list = [PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_roi_mask{:,1}];
        for iAni = 1:length(ani_list)
            % Get all you need
            idx_ani = find(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(:,1)==iPhase & abs(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).info(:,2))==ani_list(iAni));
            pocket_list = PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_roi_mask{iAni,3};
            [~,~,id_list] = unique(PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_ID(idx_ani));
            pockets_labeled = PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_roi_mask{iAni,4};
            % Prepare image
            curr_pocket_img = zeros(size(pockets_labeled));
            mask = zeros(size(pockets_labeled));
            % Iterate over all pockets
            for iP = 1:length(pocket_list)
                idx_p = find(pockets_labeled==pocket_list(iP));
                curr_pocket_img(idx_p) = id_list(iP);
                mask(idx_p) = 1;
            end%iP
            % Plot
            nexttile
            imagesc(curr_pocket_img, 'AlphaData', mask)
            % Cosmetics
            clim([1, PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k])
            colorbar
            cols = SET.HeatmapColors.(SET.phases{iPhase});
            cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), PoolData.kmeans_result.AverageValue.(SET.layers{iLayer}).kmeans_best_k+2)),:);
            cols = cols(2:end-1,:);
            colormap(cols)
            axis equal tight
            xticks([])
            yticks([])
            box on
            title(PoolData.(SET.phases{iPhase}).(SET.layers{iLayer}).info_roi_mask{iAni,2}, 'interpreter', 'none')
        end%iAni
        print(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\resp_vec_loc_',SET.phases{iPhase}], '-dpdf')
        close(hFig)
    end%iPhase
end%iLayer
clearvars -except PoolData SET
