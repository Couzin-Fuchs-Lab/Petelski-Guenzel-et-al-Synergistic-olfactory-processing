%% PureOdorAnalysis_OdorSpace
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
mkdir('PureOdorAnalysis_OdorSpace')


%% Settings

% Set paths
SET.main_path = '...\Data\Physiology\Cal520_Widefield\COL_HEX3_MOL_N2_OC3L2_ZHAE2\';

% Seth the two phases
SET.phases = {'gregarious', 'solitarious'};

% Set frames (1:100) of active region. This info comes from the other
% scripts
SET.active_region = [36 65];
SET.fps = 10;

% Clusterng & Dimensionality reduction
SET.kmeans.MaxIter = 1e4;
SET.kmeans.OnlinePhase = 'on';
SET.kmeans.Replicates = 5e1;
SET.maxk_eval = 25;
SET.kmeans.dist_RespProfile = 'sqeuclidean';
SET.kmeans.dist_AvgValue = 'sqeuclidean';
SET.dist_pca = 'correlation';
SET.pca_explained = 99;

% Set how to sort response profiles
% 'manual' or 'rise_time'
SET.sort_respprofiles = 'manual';

% Give list of stimuli (We have this in the meta data, but good to double
% check)
SET.StimList = {...
    'ZHAE2';...
    'COL';...
    'ZHAE2_COL';...
    'HEX3';...
    'ZHAE2_HEX3';...
    'OC3L2';...
    'ZHAE2_OC3L2'};

% Cosmetics
SET.Colors = [...
    200,200,200;...N2
    102,166,030;...ZHAE2
    217,095,002;...Col
    231,041,138;...ZHAE2_COL
    117,112,179;...HEX3
    27,158,119;...ZHAE2_HEX3
    230,171,002;...OC3L2
    166,118,29;...ZHAE2_OC3L2
    ]/255;

% Colormap for heatmaps
SET.HeatmapColors.gregarious = Widefield_SubFcn.ColMapPlasma(1000);
SET.HeatmapColors.solitarious = [...
    SET.HeatmapColors.gregarious(:,2),...
    SET.HeatmapColors.gregarious(:,1),...
    SET.HeatmapColors.gregarious(:,3)];
SET.HeatmapColors.other = gray(1000);


%% Pool data

% Iterate over both phases
for iPhase = 1:length(SET.phases)

    % Contruct the path
    basepath = [SET.main_path, SET.phases{iPhase},'\'];
    all_animals = dir(basepath);

    % Prepare data structures
    PoolData.(SET.phases{iPhase}).stim_dist_val = [];
    PoolData.(SET.phases{iPhase}).stim_dist_act = [];
    PoolData.(SET.phases{iPhase}).n_resp = [];
    PoolData.(SET.phases{iPhase}).n_unique = [];
    PoolData.(SET.phases{iPhase}).avg_RespVec = [];
    PoolData.(SET.phases{iPhase}).avg_RespProfile = [];
    PoolData.(SET.phases{iPhase}).avg_value = [];
    PoolData.(SET.phases{iPhase}).active_mask = [];
    PoolData.(SET.phases{iPhase}).dist2center = [];
    PoolData.(SET.phases{iPhase}).VennData = [];
    PoolData.(SET.phases{iPhase}).info = [];
    PoolData.(SET.phases{iPhase}).info_RespProfile = [];
    PoolData.(SET.phases{iPhase}).info_roi_mask = {};
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
                    PoolData.(SET.phases{iPhase}).info_RespProfile = [PoolData.(SET.phases{iPhase}).info_RespProfile; iPhase, iAni, iStim, iP];
                end%iP
                curr_tc = [curr_tc, avg_RespVec];
                PoolData.(SET.phases{iPhase}).avg_RespProfile = [PoolData.(SET.phases{iPhase}).avg_RespProfile; avg_RespVec];
            end%iStim

            % Get distances between stimuli
            % --- Based on  avg value
            dist = tril(squareform(pdist(curr_avg')));
            dist(dist==0)=NaN;
            dist = dist/nanmean(dist(:));
            PoolData.(SET.phases{iPhase}).stim_dist_val = cat(3,...
                PoolData.(SET.phases{iPhase}).stim_dist_val,...
                dist);
            % --- Based on active mask
            dist = triu(squareform(pdist(curr_active')));
            dist(dist==0)=NaN;
            dist = dist/nanmean(dist(:));
            PoolData.(SET.phases{iPhase}).stim_dist_act = cat(3,...
                PoolData.(SET.phases{iPhase}).stim_dist_act,...
                dist);

            % Extract info on how uniquely ROIs respond: Number of stimuli
            % a ROI is responding to
            curr_n_resp = hist(sum(curr_active,2),1:length(SET.StimList));
            curr_n_resp = curr_n_resp/sum(curr_n_resp);
            PoolData.(SET.phases{iPhase}).n_resp = [...
                PoolData.(SET.phases{iPhase}).n_resp;...
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
            PoolData.(SET.phases{iPhase}).n_unique = [...
                PoolData.(SET.phases{iPhase}).n_unique;...
                curr_unique];
            % --- info on the number of unique stimuli a ROI is responding to
            PoolData.(SET.phases{iPhase}).n_unique_info = curr_unique_info;
            % --- avg tc
            PoolData.(SET.phases{iPhase}).avg_RespVec = [...
                PoolData.(SET.phases{iPhase}).avg_RespVec;...
                curr_tc];
            % --- avg activity
            PoolData.(SET.phases{iPhase}).avg_value = [...
                PoolData.(SET.phases{iPhase}).avg_value;...
                curr_avg];
            % --- activity mask
            PoolData.(SET.phases{iPhase}).active_mask = [...
                PoolData.(SET.phases{iPhase}).active_mask;...
                curr_active];
            % --- information
            PoolData.(SET.phases{iPhase}).info = [...
                PoolData.(SET.phases{iPhase}).info; ...
                ones(length(curr_avg),1)*iAni, pocket_list(:)];
            % --- pocket distance to center
            PoolData.(SET.phases{iPhase}).dist2center = [...
                PoolData.(SET.phases{iPhase}).dist2center; curr_centroid];
            % --- Venn diagram proportions
            PoolData.(SET.phases{iPhase}).VennData = [...
                PoolData.(SET.phases{iPhase}).VennData; curr_Venn];
            % --- info on the mask of ROIs
            PoolData.(SET.phases{iPhase}).info_roi_mask{cnt,1} = iAni;
            PoolData.(SET.phases{iPhase}).info_roi_mask{cnt,2} = all_animals(iAni).name;
            PoolData.(SET.phases{iPhase}).info_roi_mask{cnt,3} = pocket_list;
            PoolData.(SET.phases{iPhase}).info_roi_mask{cnt,4} = Segmentation.pockets_labeled;
            cnt = cnt+1;

            % Clean
            clearvars -except PoolData iStim iAni iPhase SET basepath all_animals cnt

        end%if animal
    end%iAni
end%iPhase

%% Pool phases and get response motive IDs
% -------------------------------------------------------------------------
active_mask = [];
dat_RespProfile = [];
dat_AvgVal = [];
dat_TC = [];
anis_RespProfile = [];
anis_AvgVal = [];
phase_info = [];
respprofile_info = [];
dist2center = [];
for iPhase = 1:length(SET.phases)
    active_mask =      [active_mask;      PoolData.(SET.phases{iPhase}).active_mask];
    dat_RespProfile =  [dat_RespProfile;  PoolData.(SET.phases{iPhase}).avg_RespProfile];
    dat_AvgVal =       [dat_AvgVal;       PoolData.(SET.phases{iPhase}).avg_value];
    dat_TC =           [dat_TC;           PoolData.(SET.phases{iPhase}).avg_RespVec];
    anis_RespProfile = [anis_RespProfile; PoolData.(SET.phases{iPhase}).info_RespProfile(:,2)*(2*((iPhase==1)-0.5))];
    anis_AvgVal =      [anis_AvgVal;      PoolData.(SET.phases{iPhase}).info(:,1)*(2*((iPhase==1)-0.5))];
    phase_info =       [phase_info;       ones(size(PoolData.(SET.phases{iPhase}).avg_value,1),1)*iPhase];
    respprofile_info = [respprofile_info; PoolData.(SET.phases{iPhase}).info_RespProfile];
    dist2center =      [dist2center;      PoolData.(SET.phases{iPhase}).dist2center(:,3)];
end%iPhase
% Keep original version (not normalized) for later
dat_TC_original = dat_TC;
dat_AvgVal_original = dat_AvgVal;
dat_RespProfile_original = dat_RespProfile;



%% Normalize data from each animal
% --- response motive
ani_list = unique(anis_RespProfile);
for iAni = 1:length(ani_list)
    idx_k = find(anis_RespProfile == ani_list(iAni));
    curr_dat = dat_RespProfile(idx_k,:);
    curr_dat = curr_dat - min(curr_dat(:));
    curr_dat = sqrt(curr_dat);
    curr_dat = curr_dat - nanmean(curr_dat(:));
    curr_dat = curr_dat / nanstd(curr_dat(:));
    dat_RespProfile(idx_k,:) = curr_dat;
end%iAni

% --- Average value
ani_list = unique(anis_AvgVal);
for iAni = 1:length(ani_list)
    idx_k = find(anis_AvgVal == ani_list(iAni));
    curr_dat = dat_AvgVal(idx_k,:);
    curr_dat = curr_dat - min(curr_dat(:));
    curr_dat = sqrt(curr_dat);
    curr_dat = curr_dat - nanmean(curr_dat(:));
    curr_dat = curr_dat / nanstd(curr_dat(:));
    dat_AvgVal(idx_k,:) = curr_dat;
    curr_dat = dat_TC(idx_k,:);
    curr_dat = curr_dat - min(curr_dat(:));
    curr_dat = sqrt(curr_dat);
    curr_dat = curr_dat - nanmean(curr_dat(:));
    curr_dat = curr_dat / nanstd(curr_dat(:));
    dat_TC(idx_k,:) = curr_dat;
end%iAni


%% Cluster response motives
% Here we cluster all time courses from all ROIs responsing to all stimuli
% to check whether there are unique response motives (e.g., slow rise,
% long peak, inhibition, etc)

% Normalize each row
dat = normalize(dat_RespProfile,2);

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
    VRC(k) = Widefield_SubFcn.VarianceRatioCriterion(dat, k, idx, C);
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
elbow_x_value = Widefield_SubFcn.find_elbow_point(xData_interp, VRC_fit);
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
            plot(mean(dat_RespProfile_original(idx,:)), 'LineWidth',3)
            title(iK)
            figure(hFig2); hold on
            plot(mean(dat_RespProfile_original(idx,:)), 'LineWidth',3)
        end
        figure(hFig2); legend;
        % Ask user
        clc
        disp('Give response motive cluster order as, e.g., [1 3 2]')
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
end%switch

% Get the corresponding "real" centroids
for iPhase = 1:length(SET.phases)
    % Get list of animal and IDs
    ani_list = unique(respprofile_info(respprofile_info(:,1) == iPhase,2));
    id_list = unique(ID(respprofile_info(:,1) == iPhase));
    avg_val = nan(kmeans_result.best_k,1);
    for iK = 1:kmeans_result.best_k
        kmeans_result.centroid{iPhase}{iK} = nan(length(ani_list), size(dat_RespProfile_original,2));
        for iAni = 1:length(ani_list)
            idx_ani = find(respprofile_info(:,2) == ani_list(iAni) & ID == id_list(iK));
            kmeans_result.centroid{iPhase}{iK}(iAni,:) = nanmean(dat_RespProfile_original(idx_ani,:),1);
        end%iAni
        avg_val(iK) = nanmean(nanmean(kmeans_result.centroid{iPhase}{iK}(:,SET.active_region(1):SET.active_region(2))));
    end%iK
end%iPhase

% Save
PoolData.kmeans_result.ResponseProfile.best_k = kmeans_result.best_k;
PoolData.kmeans_result.ResponseProfile.ID = ID;
PoolData.kmeans_result.ResponseProfile.dat_RespProfile = dat_RespProfile_original;
PoolData.kmeans_result.ResponseProfile.centroid = kmeans_result.centroid;
PoolData.kmeans_result.ResponseProfile.info = respprofile_info;

% Clean up
clearvars -except PoolData SET dat_* anis_* phase_info respprofile_info active_mask dist2center


%% Cluster Response Vectors
% Here, we cluster ROIs based on their response motives to each stimulus.
% We do not take average values only since the profile offers considerable
% information here.

% Normalize each row
dat = normalize(dat_AvgVal,2);
dat_avg = normalize(dat_AvgVal,2);

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
    VRC(k) = Widefield_SubFcn.VarianceRatioCriterion(dat, k, idx, C);
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
elbow_x_value = Widefield_SubFcn.find_elbow_point(xData_interp, VRC_fit);
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
    ani_list = unique(anis_AvgVal(phase_info==iPhase));
    id_list = unique(ID(phase_info==iPhase));
    % Iterate over all clusters
    for iK = 1:kmeans_result.best_k
        % Preallocation
        kmeans_result.centroid{iPhase}{iK} = nan(length(ani_list), size(dat_AvgVal_original,2));
        kmeans_result.centroid_pca{iPhase}{iK} = nan(length(ani_list), 2);
        % Iterate over all animal
        for iAni = 1:length(ani_list)
            idx_ani = find(anis_AvgVal == ani_list(iAni) & ID == id_list(iK));
            kmeans_result.centroid{iPhase}{iK}(iAni,:) =       nanmean(dat_AvgVal_original(idx_ani,:),1);
            kmeans_result.centroid_norm{iPhase}{iK}(iAni,:) =  nanmean(dat_avg(idx_ani,:),1);
            kmeans_result.centroid_pca{iPhase}{iK}(iAni,:) =   nanmean(dat(idx_ani,1:2),1);
            dist2center_sorted{iPhase}{iK}(iAni,1) =           nanmean(dist2center(idx_ani));
        end%iAni
    end%iK
end%iPhase

% Save
% --- The extracted data
PoolData.kmeans_result.AverageValue.dat_active_mask = active_mask;      % whether the ROI was actvive or not
PoolData.kmeans_result.AverageValue.dat_TC = dat_TC_original;           % time traces (not normalized)
PoolData.kmeans_result.AverageValue.dat_AvgVal = dat_AvgVal_original;   % average values (not normalized)
PoolData.kmeans_result.AverageValue.dat_pca = dat;                      % the PCA pre-proceed data used for clustering
PoolData.kmeans_result.AverageValue.dat_pca_explained = explained;      % the PCA pre-proceed data used for clustering
% --- Each ROIs distance to the respective AL's center
PoolData.kmeans_result.AverageValue.dat_dist2center = dist2center;
PoolData.kmeans_result.AverageValue.dist2center_sorted = dist2center_sorted;
% --- The clustring results
PoolData.kmeans_result.AverageValue.kmeans_best_k = kmeans_result.best_k;                % optimal number of clusters
PoolData.kmeans_result.AverageValue.kmeans_ID = ID;                                      % each ROI's assigned ID
PoolData.kmeans_result.AverageValue.kmeans_centroid = kmeans_result.centroid;            % each animal's avg value centroid (not normalized)
PoolData.kmeans_result.AverageValue.kmeans_centroid_norm = kmeans_result.centroid_norm;  % each animal's avg value centroid (not normalized)
PoolData.kmeans_result.AverageValue.kmeans_centroid_pca = kmeans_result.centroid_pca;    % each animal's centroid in PCA space (pca-1 & pca-2)
% --- Information on each ROI's phase (col-1) and animal (col-2)
PoolData.kmeans_result.AverageValue.info = [phase_info, anis_AvgVal];

clearvars -except PoolData SET dat_* anis_* phase_info respprofile_info active_mask dist2center

% Save and clean up
save('PureOdorAnalysis_OdorSpace\OdorSpace_kmeans_PoolData.mat', 'PoolData', 'SET')
clearvars -except PoolData SET


%% Visualize the results.
% response motive vectors

% Define the range of the y-axis
ylim_range = [-1 2.5];

hFig = figure('units', 'normalized', 'Position', [0 0 0.5 0.5]);
% Iterate over both phases and plot the average motives for each phase
% in a subplot
for iPhase = 1:length(SET.phases)
    % One subplot per phase
    subplot(1,3,iPhase) ; hold on
    % Get colors
    cols = SET.HeatmapColors.(SET.phases{iPhase});
    cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), PoolData.kmeans_result.ResponseProfile.best_k+2)),:);
    cols = cols(2:end-1,:);
    % Get time vector
    xvec = linspace(0, size(PoolData.kmeans_result.ResponseProfile.dat_RespProfile,2)/SET.fps, size(PoolData.kmeans_result.ResponseProfile.dat_RespProfile,2));
    % Get position of active window
    pos = [xvec(SET.active_region(1)), -1, xvec(SET.active_region(2))-xvec(SET.active_region(1)),3];
    % Plot active window as gray box
    rectangle('position', pos, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
    % Iterate over all motives
    for iK = 1:PoolData.kmeans_result.ResponseProfile.best_k
        % Get the average
        avg = nanmean(PoolData.kmeans_result.ResponseProfile.centroid{iPhase}{iK},1);
        % Get the corresponding confidence interval
        ci = bootci(5000, {@nanmean, PoolData.kmeans_result.ResponseProfile.centroid{iPhase}{iK}});
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
%Also, plot the motives for both phases combined in the third subplot
subplot(1,3,3) ; hold on
% Use a different colormap
cols = Widefield_SubFcn.ColMapInferno(1000);
cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), PoolData.kmeans_result.ResponseProfile.best_k+2)),:);
cols = cols(2:end-1,:);
% Get the time vector
xvec = linspace(0, size(PoolData.kmeans_result.ResponseProfile.dat_RespProfile,2)/SET.fps, size(PoolData.kmeans_result.ResponseProfile.dat_RespProfile,2));
% Get the position of the active region
pos = [xvec(SET.active_region(1)), -1, xvec(SET.active_region(2))-xvec(SET.active_region(1)),3];
% Plot the active region as a gray box
rectangle('position', pos, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
% Iterate over all motives
for iK = 1:PoolData.kmeans_result.ResponseProfile.best_k
    % Pool gregarious and solitarious data
    currDat = [];
    for iPhase = 1:length(SET.phases)
        currDat = [currDat; PoolData.kmeans_result.ResponseProfile.centroid{iPhase}{iK}];
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
ylim(ylim_range)
ylabel('\DeltaF/F (%)')
xlabel('time (s)')
title('both')
% Save figure and close it
export_fig('PureOdorAnalysis_OdorSpace\kmeans_ResponseProfile', '-pdf', '-painters')
close(hFig)
clearvars -except PoolData SET


%% Visualize the results.
% Number stimuli a ROI responds to

hFig = figure('Name', 'n_resp'); hold on
% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Get colors
    col = SET.HeatmapColors.(SET.phases{iPhase}); col = col(floor(size(col,1)/2),:);
    % Get average and CI across animals
    avg = nanmean(bootstrp(5000, @nanmean, PoolData.(SET.phases{iPhase}).n_resp));
    ci = bootci(5000, @nanmean, PoolData.(SET.phases{iPhase}).n_resp);
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
export_fig('PureOdorAnalysis_OdorSpace\kmeans_n_resp', '-pdf', '-painters')
close(hFig)
clearvars -except PoolData SET


%% Visualize the results.
% Proportion of ROIS responding to unique stimuli

hFig = figure('Name', 'n_unique'); hold on
% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Get color
    col = SET.HeatmapColors.(SET.phases{iPhase}); col = col(floor(size(col,1)/2),:);
    % Get average and CI across animals
    avg = nanmean(bootstrp(5000, @nanmean, PoolData.(SET.phases{iPhase}).n_unique));
    ci = bootci(5000, @nanmean, PoolData.(SET.phases{iPhase}).n_unique);
    % Iterate over everything
    for iStim = 1:size(PoolData.(SET.phases{iPhase}).n_unique,2)
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
        text(iStim, 0.175, PoolData.(SET.phases{iPhase}).n_unique_info{1,iStim}, 'horizontalAlignment', 'center', 'FontSize', 4, 'Color',[0 1 0],"Interpreter","none")
        text(iStim, 0.17,   PoolData.(SET.phases{iPhase}).n_unique_info{2,iStim}, 'horizontalAlignment', 'center', 'FontSize', 4, 'Color',[1 0 0],"Interpreter","none")
    end%iStim
end%iPhase
% Cosmetics
xlim([0.5 size(PoolData.(SET.phases{iPhase}).n_unique,2)+0.5])
ylim([0 0.5])
set(gca, 'xtick',1:size(PoolData.(SET.phases{iPhase}).n_unique,2) )
ylabel('prop. ROIs')
title('Proportion of ROIs uniquely responding to given stimulus')
% Save
export_fig('PureOdorAnalysis_OdorSpace\kmeans_n_unique', '-pdf', '-painters')
close(hFig)

% Statistics. Note that this can take considerable time.
PoolData.Stats.seed=1234;
PoolData.Stats.n = 5e7;
PoolData.Stats.nComp = 1;
% Preallocation
PoolData.Stats.kmeans_n_unique.p = nan(size(PoolData.(SET.phases{1}).n_unique,2),1);
PoolData.Stats.kmeans_n_unique.s = nan(size(PoolData.(SET.phases{1}).n_unique,2),1);
PoolData.Stats.kmeans_n_unique.c = nan(size(PoolData.(SET.phases{1}).n_unique,2),1);
% Iterate over all stimuli
for iStim = 1:size(PoolData.(SET.phases{1}).n_unique,2)
    [...
        PoolData.Stats.kmeans_n_unique.p(iStim),...
        PoolData.Stats.kmeans_n_unique.s(iStim),...
        ~,...
        PoolData.Stats.kmeans_n_unique.c(iStim)] = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample', PoolData.(SET.phases{1}).n_unique(:,iStim), PoolData.(SET.phases{2}).n_unique(:,iStim), PoolData.Stats.n, PoolData.Stats.seed, PoolData.Stats.nComp);
end%iStim

%Save
save('PureOdorAnalysis_OdorSpace\OdorSpace_kmeans_PoolData.mat', 'PoolData', 'SET')

clearvars -except PoolData SET


%% Visualize the results.
% Venn diagram for proportion of ROIS responding to combinations of stimuli

% Give names for the different regions (hard coded)
venn_names = {'Laa only', 'Lct only', 'LaaLct only', 'Laa and Lct', 'Laa and LaaLct', 'Lct and LaaLct', 'all'};

hFig = figure('Name', 'venn'); hold on
% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Plot diagram in one subplot. Provide info further down (subplot
    % below)
    subplot(2,2,iPhase); hold on
    % Get the average
    avg = mean(PoolData.(SET.phases{iPhase}).VennData);
    % And the corresponding confidence interval
    ci = bootci(5000, @mean, PoolData.(SET.phases{iPhase}).VennData);
    % Make sure everything adds up to one
    avg_norm = avg/sum(avg);
    % Plot the Venn diagram. Here, provide the proportions for every
    % region (only_A, only_B, only_C, A_and_B, A_and_C, B_and_C, A_B_C)
    Widefield_SubFcn.draw_venn_diagram(avg_norm(1), avg_norm(2), avg_norm(3), avg_norm(4), avg_norm(5), avg_norm(6), avg_norm(7),...
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
export_fig('PureOdorAnalysis_OdorSpace\venn_diagrams', '-pdf', '-painters')
close(hFig)

% Do statistics. Note that this can take some time
PoolData.Stats.seed=1234;
PoolData.Stats.n = 5e7;
PoolData.Stats.nComp = 1;
for iVenn = 1:length(venn_names)
    [...
        PoolData.Stats.VennDiagram.p(iVenn),...
        PoolData.Stats.VennDiagram.s(iVenn),...
        ~,...
        PoolData.Stats.VennDiagram.c(iVenn)] = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample', PoolData.(SET.phases{1}).VennData(:,iVenn), PoolData.(SET.phases{2}).VennData(:,iVenn), PoolData.Stats.n, PoolData.Stats.seed, PoolData.Stats.nComp);
end%iVenn

%Save
save('PureOdorAnalysis_OdorSpace\OdorSpace_kmeans_PoolData.mat', 'PoolData', 'SET')

clearvars -except PoolData SET


%% Visualize the results.
% Distance between stimuli based on amplitude (top left) or on active
% regions (bottom right)

% Iterate over both phases
for iPhase = 1:length(SET.phases)
    hFig = figure('Name', 'stim_dist', 'units', 'normalized', 'Position', [0 0 0.5 0.5], 'Color', 'w');

    % Distance between stimuli based on their amplitudes
    subplot(1,2,1); hold on
    % Get the average
    avg1 = nanmean(PoolData.(SET.phases{iPhase}).stim_dist_val,3);
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
    avg2 = nanmean(PoolData.(SET.phases{iPhase}).stim_dist_act,3);
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
    export_fig(['PureOdorAnalysis_OdorSpace\kmeans_stim_dist_',SET.phases{iPhase}], '-pdf', '-painters')
    close(hFig)

end%iPhase
clearvars -except PoolData SET


%% Visualize the results.
% Animals per response motive

hFig = figure('Units','normalized','Position',[0 0 1 1]);
% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Get list of respective cluster IDs for this phase
    id_list = unique( PoolData.kmeans_result.ResponseProfile.ID(PoolData.kmeans_result.ResponseProfile.info(:,1) == iPhase));
    % Get list of all possible animals for this phase
    ani_list = unique(PoolData.kmeans_result.ResponseProfile.info(PoolData.kmeans_result.ResponseProfile.info(:,1) == iPhase, 2));
    % Iterate over all clusters
    for iK = 1:length(id_list)
        % Get index positions and observations
        idx = find(PoolData.kmeans_result.ResponseProfile.ID == id_list(iK));
        anis = PoolData.kmeans_result.ResponseProfile.info(idx,2);
        % Count how many times each animal was present
        cnts = histcounts(anis, [ani_list(:); max(ani_list)+1]);
        % Normalize by number of ROIs per animal
        cnts = cnts./histcounts(PoolData.kmeans_result.ResponseProfile.info(:,2), [ani_list(:); max(ani_list)+1]);
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
        plot(nanmean(PoolData.kmeans_result.ResponseProfile.centroid{iPhase}{iK}), 'color', cols(iK,:), 'LineWidth', 2)
        title(['cluster ',num2str(iK)])
        axis tight off
    end%iC
end%iPhase
% Save
export_fig('PureOdorAnalysis_OdorSpace\pie_RespProfile', '-pdf', '-painters')
close(hFig)
clearvars -except PoolData SET


%% Visualize the results.
% Animals per response vector cluster

hFig = figure('Units','normalized','Position',[0 0 1 1]);
% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Get list of respective cluster IDs
    id_list = unique(PoolData.kmeans_result.AverageValue.kmeans_ID(PoolData.kmeans_result.AverageValue.info(:,1) == iPhase));
    % Get list of all possible animals
    ani_list = unique(PoolData.kmeans_result.AverageValue.info(PoolData.kmeans_result.AverageValue.info(:,1) == iPhase,2));
    % Iterate over all clusters
    for iK = 1:length(id_list)
        % Get index positions and observations
        idx = find(PoolData.kmeans_result.AverageValue.kmeans_ID == id_list(iK));
        anis = PoolData.kmeans_result.AverageValue.info(idx,2);
        % Count how many times each animal was present
        cnts = histcounts(anis, [ani_list(:); max(ani_list)+1]);
        % Normalize by number of ROIs per animal
        cnts = cnts./histcounts(PoolData.kmeans_result.AverageValue.info(:,2), [ani_list(:); max(ani_list)+1]);
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
export_fig('PureOdorAnalysis_OdorSpace\pie_RespVec', '-pdf', '-painters')
close(hFig)
clearvars -except PoolData SET


%% Visualize the results.
% Show PCA olfaction space and clustering

hFig = figure('Units','normalized','Position',[0 0 1 1]);

% Bootstrap PCA data to generate p-confidence ellipse around the
% resulting distribution
% --- Keep track of the data's range
xaxis_lim = []; yaxis_lim = [];
dist2center_vec = [];
% --- Iterate over all clusters
for iK = 1:PoolData.kmeans_result.AverageValue.kmeans_best_k
    % Pool data from both phases
    curr_dat = [];
    for iPhase = 1:length(SET.phases)
        curr_dat = [curr_dat; PoolData.kmeans_result.AverageValue.kmeans_centroid_pca{iPhase}{iK}];
    end%iPhase
    %Kick out NaNs
    curr_dat = curr_dat(~all(isnan(curr_dat), 2),:);
    % Bootstrap data for the ellipse
    dat_pca_boot = bootstrp(5000, @nanmean, curr_dat);
    % Get a 99%-confidence ellipse around the bootstrapped data
    r_ellipse{iK} = Widefield_SubFcn.dataEllipse(dat_pca_boot, 0.99);
    % Get the range of the data
    xaxis_lim = [xaxis_lim; min(curr_dat(:,1)), max(curr_dat(:,1))];
    yaxis_lim = [yaxis_lim; min(curr_dat(:,2)), max(curr_dat(:,2))];
    dist2center_vec = [dist2center_vec; min(PoolData.kmeans_result.AverageValue.dist2center_sorted{iPhase}{iK}(:)), max(PoolData.kmeans_result.AverageValue.dist2center_sorted{iPhase}{iK}(:))];
end%iK
% Cosmetics
xaxis_lim = [min(xaxis_lim(:)), max(xaxis_lim(:))]*1.15;
yaxis_lim = [min(yaxis_lim(:)), max(yaxis_lim(:))]*1.15;
dist2center_vec = [min(dist2center_vec(:)), max(dist2center_vec(:))];

% Raw PCA space (both phases combined) --------------------------------
nexttile; hold on
% Get nice colors
cols = Widefield_SubFcn.ColMapInferno(size(SET.HeatmapColors.(SET.phases{iPhase}),1));
cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1),sum(PoolData.kmeans_result.AverageValue.kmeans_best_k)+2)),:);
cols = cols(2:end-1,:);
% Iterate over all clusters
for iK = 1:PoolData.kmeans_result.AverageValue.kmeans_best_k
    % Pool data from both phases
    curr_dat = [];
    for iPhase = 1:length(SET.phases)
        curr_dat = [curr_dat; PoolData.kmeans_result.AverageValue.kmeans_centroid_pca{iPhase}{iK}];
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
cols = Widefield_SubFcn.ColMapInferno(size(SET.HeatmapColors.(SET.phases{iPhase}),1));
cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1),sum(PoolData.kmeans_result.AverageValue.kmeans_best_k)+2)),:);
cols = cols(2:end-1,:);
% Iterate over all clusters
for iK = 1:PoolData.kmeans_result.AverageValue.kmeans_best_k
    % Iterate over both phases and get data
    curr_dat = []; curr_phase = []; col = [];
    for iPhase = 1:length(SET.phases)
        % Get data
        curr_dat = [curr_dat; PoolData.kmeans_result.AverageValue.kmeans_centroid_pca{iPhase}{iK}];
        curr_phase = [curr_phase; ones(size(PoolData.kmeans_result.AverageValue.kmeans_centroid_pca{iPhase}{iK},1),1)*iPhase];
        % Get nice colors
        col = [col;...
            ones(size(PoolData.kmeans_result.AverageValue.kmeans_centroid_pca{iPhase}{iK},1),1)*0.75*(iPhase==1),...
            zeros(size(PoolData.kmeans_result.AverageValue.kmeans_centroid_pca{iPhase}{iK},1),1),...
            ones(size(PoolData.kmeans_result.AverageValue.kmeans_centroid_pca{iPhase}{iK},1),1)*0.75*(iPhase==2)];
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
cols = Widefield_SubFcn.ColMapInferno(size(SET.HeatmapColors.(SET.phases{iPhase}),1));
dist2center_vec = linspace(min(dist2center_vec), max(dist2center_vec), size(cols,1));
% Iterate over all cluster
for iK = 1:PoolData.kmeans_result.AverageValue.kmeans_best_k
    % Pool phases
    curr_dat = [];
    curr_dist = [];
    for iPhase = 1:length(SET.phases)
        curr_dat = [curr_dat; PoolData.kmeans_result.AverageValue.kmeans_centroid_pca{iPhase}{iK}];
        curr_dist = [curr_dist; PoolData.kmeans_result.AverageValue.dist2center_sorted{iPhase}{iK}];
    end%iPhase
    % Plot elipse
    [~,d_index] = min(abs(dist2center_vec - nanmean(PoolData.kmeans_result.AverageValue.dist2center_sorted{iPhase}{iK})));
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
    curr_dat = [curr_dat; cell2mat(PoolData.kmeans_result.AverageValue.dist2center_sorted{iPhase})];
end%iPhase
% Do stats
% --- ANOVA
[p,~,aov] = anova1(curr_dat, 1:PoolData.kmeans_result.AverageValue.kmeans_best_k,'off');
% --- Multiple comparison
comp_results = multcompare(aov, 'CriticalValueType', 'hsd', 'Display', 'off');
comp_results = sortrows(comp_results);
comp_results = array2table(comp_results, 'VariableNames', {'Group1', 'Group2', 'MeanDifferenceLower', 'MeanDifference', 'MeanDifferenceUpper', 'pValue'});
% --- Extract the comparison group indices and significant differences
CLD = Widefield_SubFcn.multcompLetters(comp_results);
PoolData.Stats.(SET.phases{iPhase}).dist2center_vector.anova1_p = p;
PoolData.Stats.(SET.phases{iPhase}).dist2center_vector.comp_results = comp_results;
PoolData.Stats.(SET.phases{iPhase}).dist2center_vector.CLD = CLD;

% Properties for the violin plot
properties.MinVal =             [];                 % Smallest possible value (e.g. errors = 0)
properties.MaxVal =             [];                 % Biggest possible value
properties.AvgType =            'mean';             %(Set which measure should be plotted: 'median', 'mean' or 'both')
properties.MeanWidth =          2;                  %(Set width of the median's line and/or symbol)
properties.SeparateOutliers =   0;                  %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
properties.MeanCol =            [0.251 0.251 0.251];
% Iterate over all clusters
for iK = 1:PoolData.kmeans_result.AverageValue.kmeans_best_k
    Widefield_SubFcn.violinplot_advanced(curr_dat(:,iK), iK, 0.75, properties)
    text(iK, 0.1, CLD{iK}, 'color', [1 0 0], 'HorizontalAlignment', 'center')
end%iK
% Get the distance of an uniformly distributed sample
dat = rand(1e6,2)*2;
dat = dat-(mean(dat));
% Get each point's distance
dat = sqrt(sum(dat'.*dat'))';
dat(dat>1) = [];
plot([0, PoolData.kmeans_result.AverageValue.kmeans_best_k+1], [1 1]*mean(dat))
clear dat
% Cosmetics
title(['avg dist2center | ',SET.phases{iPhase}])
set(gca, 'XTick', 1:PoolData.kmeans_result.AverageValue.kmeans_best_k)
ylabel('norm. distance to center')
xlabel('response vector ID')
ylim([0 1])
xlim([0, PoolData.kmeans_result.AverageValue.kmeans_best_k+1])

% Export & Save -----------------------------------------------------------
export_fig('PureOdorAnalysis_OdorSpace\kmeans_PCA_space', '-pdf', '-painters')
save('PureOdorAnalysis_OdorSpace\OdorSpace_kmeans_PoolData.mat', 'PoolData', 'SET')
close(hFig)


%% Visualize the results.
% Show how the clusters are different from each other (response vectors)

% Set the ylim (zscored)
ylim_range = [-2 2];


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
for iK = 1:PoolData.kmeans_result.AverageValue.kmeans_best_k

    % Pool data across phases
    resp_vec = [];
    temp_PCA_centroid = [];
    temp_PCA_centroid_pca = [];
    for iPhase = 1:length(SET.phases)
        % resp_vec = [resp_vec; PoolData.kmeans_result.AverageValue.kmeans_centroid{iPhase}{iK}];
        resp_vec =              [resp_vec;              PoolData.kmeans_result.AverageValue.kmeans_centroid_norm{iPhase}{iK}];
        temp_PCA_centroid =     [temp_PCA_centroid;     nanmean(PoolData.kmeans_result.AverageValue.kmeans_centroid{iPhase}{iK})];
        temp_PCA_centroid_pca = [temp_PCA_centroid_pca; nanmean(PoolData.kmeans_result.AverageValue.kmeans_centroid_pca{iPhase}{iK})];
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
    CLD = Widefield_SubFcn.multcompLetters(comp_results);

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
        prop = length(find(PoolData.kmeans_result.AverageValue.kmeans_ID == iK & PoolData.kmeans_result.AverageValue.info(:,1) == iPhase)) / sum(PoolData.kmeans_result.AverageValue.info(:,1) == iPhase);
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
export_fig('PureOdorAnalysis_OdorSpace\resp_vec_porp_pie', '-pdf', '-painters')
close(hFig_pi)
% --- cluster location
figure(hFig_cluster)
axis equal
xlim(xaxis_lim)
ylim(yaxis_lim)
xticks([])
yticks([])
box on
export_fig('PureOdorAnalysis_OdorSpace\kmeans_cluster_locations', '-pdf', '-painters')
close(hFig_cluster)
% --- response vectors
figure(hFig_vec)
nexttile
D = squareform(pdist(PCA_centroid, SET.dist_pca));
tree = linkage(D,'single', SET.dist_pca);
dendrogram(tree, 'Labels', PCA_label)
export_fig('PureOdorAnalysis_OdorSpace\kmeans_cluster_responses', '-pdf', '-painters')
close(hFig_vec)

% Clean
clearvars -except PoolData SET


%% Visualize the results.
% Show how response motive assignments change between stimuli

% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Get avg values for each ROI (row) and each stimulus (column)
    data = PoolData.kmeans_result.AverageValue.dat_AvgVal(PoolData.kmeans_result.AverageValue.info(:,1) == iPhase, :);
    % Get corresponding grouping factors (response motives)
    motives = PoolData.kmeans_result.ResponseProfile.ID(PoolData.kmeans_result.ResponseProfile.info(:,1) == iPhase,:);
    stims = PoolData.kmeans_result.ResponseProfile.info(PoolData.kmeans_result.ResponseProfile.info(:,1) == iPhase,3);
    data_grouping = [];
    % Iterate ove all stimuli and populate the grouping variable
    for iStim = 1:length(SET.StimList)
        idx = find(stims == iStim);
        data_grouping = [data_grouping, motives(idx)];
    end
    % Plot response hive (circular stacked bar plot)
    hFig = figure;
    Widefield_SubFcn.HivePlot(data(:,1:3), data_grouping(:,1:3), [0 0 0], [], 0.5, length(unique(PoolData.kmeans_result.ResponseProfile.ID(PoolData.kmeans_result.ResponseProfile.info(:,1) == iPhase,:)))*0.5)
    % Add a colorbar for the different levels
    cols = SET.HeatmapColors.(SET.phases{iPhase});
    cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), PoolData.kmeans_result.ResponseProfile.best_k+2)),:);
    cols = cols(2:end-1,:);
    colormap(cols); colorbar; clim([1 PoolData.kmeans_result.ResponseProfile.best_k])
    % Add labels
    set(gca,'ThetaTickLabel', SET.StimList, 'TickLabelInterpreter', 'none')
    % Save
    export_fig(['PureOdorAnalysis_OdorSpace\response_hive_',SET.phases{iPhase}], '-pdf', '-painters')
    close(hFig)

    % Also, pool proportions
    % Get list of animals
    ani_list = unique(PoolData.(SET.phases{iPhase}).info(:,1));
    % Get list of motives
    motive_list = unique(data_grouping);
    % Iterate over all but the first motive. The first motive mostly
    % contains noise
    for iMotive = 1:length(motive_list)
        % Prepare variables
        PropSun = {};
        cnt = 1;
        % Iterate over all stimuli in a nested for-loop in order to get
        % all possbile combinations of stimulus changes
        for iStim = 1:length(SET.StimList)
            for jStim = 1:length(SET.StimList)
                if iStim ~= jStim
                    % Iterate all motives again to get all possible
                    % motive changes. Here, also inlcude whether a ROI
                    % kept the same motive despite the stimulus
                    % changing
                    for jMotive = 1:length(motive_list)
                        % Get the proportion of ROIs that was in motive
                        % iMotive during stimulus iStim and in motive
                        % jMotive during stimulus jStim
                        cnts_pool = [];
                        p_pool = [];
                        for iAni = 1:length(ani_list)
                            cnts = sum(data_grouping(:,iStim)==motive_list(iMotive) & data_grouping(:,jStim)==motive_list(jMotive) & PoolData.(SET.phases{iPhase}).info(:,1) == ani_list(iAni));
                            % Check whether this was the case for this
                            % animal
                            if sum(data_grouping(:,iStim)==motive_list(iMotive) & PoolData.(SET.phases{iPhase}).info(:,1) == ani_list(iAni)) == 0
                                p = NaN;
                            else
                                p = cnts / sum(data_grouping(:,iStim)==motive_list(iMotive) & PoolData.(SET.phases{iPhase}).info(:,1) == ani_list(iAni));
                            end
                            cnts_pool = [cnts_pool; cnts];
                            p_pool = [p_pool; p];
                        end%iAni
                        % Put everything together
                        PropSun{cnt,1} = iStim;
                        PropSun{cnt,2} = jStim;
                        PropSun{cnt,3} = motive_list(jMotive);
                        PropSun{cnt,4} = nanmean(p_pool);
                        PropSun{cnt,5} = p_pool;
                        PropSun{cnt,6} = mean(cnts_pool);
                        PropSun{cnt,7} = cnts_pool;
                        cnt = cnt+1;
                    end%iProfile
                end%if
            end%jStim
        end%iStim
        % Normalize for each animal
        t = cell2mat(PropSun(:,1:3));
        % Again, iterate over all stimuli in a nested fashion
        for iStim = 1:length(SET.StimList)
            for jStim = 1:length(SET.StimList)
                if iStim ~= jStim
                    % Get the correct indices
                    idx = find(t(:,1) == iStim & t(:,2) == jStim);
                    % Iterate over all animals
                    for iAni = 1:length(ani_list)
                        p = [];
                        for iM = 1:length(idx)
                            p = [p; PropSun{idx(iM),5}(iAni)];
                        end%iM
                        p = sum(p);
                        for iM = 1:length(idx)
                            PropSun{idx(iM),5}(iAni) = PropSun{idx(iM),5}(iAni)/p;
                        end%iM
                    end%iAni
                    for iM = 1:length(idx)
                        PropSun{idx(iM),4} = nanmean(PropSun{idx(iM),5});
                    end%iM
                end%if
            end%jStim
        end%iStim
        % Pool results
        PoolData.(SET.phases{iPhase}).PropSun.(['motive_',num2str(motive_list(iMotive))]) = PropSun;
        clear PropSun
    end%iMotive
end%iPhase

% Do stats on the PropSun, comparing greg and solitarious
PoolData.PropSun.seed=1234;
PoolData.PropSun.n = 5e7;
PoolData.PropSun.nComp = 1;
% Also pool everything to depict as heatmap
PoolData.gregarious.PropSun.Matrix = [];
PoolData.solitarious.PropSun.Matrix = [];
% Iterate over all motives
for iMotive = 1:length(motive_list)
    % Iterate over all entries for this motive
    for iSun = 1:size(PoolData.gregarious.PropSun.(['motive_',num2str(motive_list(iMotive))]), 1)
        % [Hard coded] get data for gregarious and solitarious animals
        greg = PoolData.gregarious.PropSun.(['motive_',num2str(motive_list(iMotive))]){iSun,5};
        soli = PoolData.solitarious.PropSun.(['motive_',num2str(motive_list(iMotive))]){iSun,5};
        greg(isnan(greg)) = [];
        soli(isnan(soli)) = [];
        % Get p-value (only do the bootstrap hypothesis testing if a
        % normal ranksum indicates a trend
        if ranksum(greg, soli)<0.5
            p = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample', greg, soli, PoolData.PropSun.n, PoolData.PropSun.seed, PoolData.PropSun.nComp);
        else
            p=1;
        end
        PoolData.gregarious.PropSun.(['motive_',num2str(motive_list(iMotive))]){iSun,8} = p;
        PoolData.solitarious.PropSun.(['motive_',num2str(motive_list(iMotive))]){iSun,8} = p;
        % Check whether it is significant
        if p<0.05
            h = 1;
        elseif p<0.01
            h = 2;
        elseif p<0.001
            h = 3;
        else
            h = 0;
        end%if
        % Add the results to the big overview
        PoolData.gregarious.PropSun.(['motive_',num2str(motive_list(iMotive))]){iSun,9} = h;
        PoolData.solitarious.PropSun.(['motive_',num2str(motive_list(iMotive))]){iSun,9} = h;
        % Pool
        PoolData.gregarious.PropSun.Matrix = [PoolData.gregarious.PropSun.Matrix;...
            motive_list(iMotive), [PoolData.gregarious.PropSun.(['motive_',num2str(motive_list(iMotive))]){iSun,[1:4,6,8:9]}]];
        PoolData.solitarious.PropSun.Matrix = [PoolData.solitarious.PropSun.Matrix;...
            motive_list(iMotive), [PoolData.solitarious.PropSun.(['motive_',num2str(motive_list(iMotive))]){iSun,[1:4,6,8:9]}]];
    end%iSun
end%iMotive

% Get a list of all stimulus changes
stim_change_list = fliplr(unique(PoolData.gregarious.PropSun.Matrix(:,[2 3]),'rows'));
% Get list of all motive changed
motive_change_list = fliplr(unique(PoolData.gregarious.PropSun.Matrix(:,[1 4]),'rows'));

% Create heatmaps that illustrate the proportions of ROIs that change
% from one motive to another given a change in stimulus.
% First preallocate some space
map = zeros(size(stim_change_list,1), size(motive_change_list,1), 3);
motive_labels = cell(1, size(motive_change_list,1));
stim_labels = cell(1, size(stim_change_list,1));
% Iterate over all stimulus changes
for iStim = 1:size(stim_change_list,1)
    % Iterate over all motive changes
    for iMotive = 1:size(motive_change_list,1)
        % Get the correct index
        idx = find(...
            PoolData.gregarious.PropSun.Matrix(:,2) == stim_change_list(iStim,1) &...
            PoolData.gregarious.PropSun.Matrix(:,3) == stim_change_list(iStim,2) &...
            PoolData.gregarious.PropSun.Matrix(:,1) == motive_change_list(iMotive,1) &...
            PoolData.gregarious.PropSun.Matrix(:,4) == motive_change_list(iMotive,2));
        % Extract data for gregarious and solitarious
        map(iStim,iMotive,1) = PoolData.gregarious.PropSun.Matrix(idx,5);
        map(iStim,iMotive,2) = PoolData.solitarious.PropSun.Matrix(idx,5);
        % For plottnig, create corresponding labels
        motive_labels{iMotive} = [num2str(motive_change_list(iMotive,1)),'>',num2str(motive_change_list(iMotive,2))];
    end%iMotive
    stim_labels{iStim} = [num2str(stim_change_list(iStim,1)),'>',num2str(stim_change_list(iStim,2))];
end%iStim
% Also get the difference between gregarious and solitarious animals
map(:,:,3) = map(:,:,1)-map(:,:,2);

% Get a nice diverging colormap for the difference between solitarious
% (blue, negative) and gregariuos (red, positive) animals
temp_col = map(:, :, 3);
cols = [...
    interp1([min(temp_col(:)) 0 max(temp_col(:))], [0, 1, 0.75], linspace(min(temp_col(:)),max(temp_col(:)),1000))',...
    interp1([min(temp_col(:)) 0 max(temp_col(:))], [0, 1, 0], linspace(min(temp_col(:)),max(temp_col(:)),1000))',...
    interp1([min(temp_col(:)) 0 max(temp_col(:))], [0.75, 1, 0], linspace(min(temp_col(:)),max(temp_col(:)),1000))'];

% One figure to show all three (gregarious, solitarious, and the
% difference between the phases)
hFig = figure('color', 'w', 'units', 'normalized', 'Position', [0 0 1 1 ]);
% --- Gregarious ---
% Show heatmap containing all proportions of changing motives for all
% stimulus changes
subplot(3,3,1)
imagesc(map(:, :, 1))
% Cosmetcs
axis equal tight
clim([0, 0.45]); colorbar
colormap(gca, SET.HeatmapColors.gregarious)
xticks(1:size(motive_change_list,1))
yticks(1:size(stim_change_list,1))
set(gca, 'XTickLabel', motive_labels, 'YTickLabel', stim_labels)
xtickangle(90)
xlabel('motive change')
ylabel('stim. change')
title('gregarious')

% Show the 5 most prominent motive changes when changing
% from Laa to LaaLct
subplot(3,3,2); hold on
% Sort proportions
[val, idx] = sort(map(find(stim_change_list(:,1)==1 & stim_change_list(:,2)==3), :, 1), 'descend');
% Keep the 5 highest only
val = val(1:5);
idx = idx(1:5);
b = bar(1:length(idx), val); b.FaceColor = [0.75 0 0]; b.EdgeColor = 'none';
set(gca, 'XTick', 1:length(idx), 'XTickLabel', motive_labels(idx))
xtickangle(90)
xlabel('motive change')
ylabel('proportion')
title('S1 -> S3')
xlim([0 length(idx)+1])
ylim([0 0.6])

% Show the 5 most prominent motive changes when changing
% from Lct to LaaLct
subplot(3,3,3); hold on
% Sort proportions
[val, idx] = sort(map(find(stim_change_list(:,1)==2 & stim_change_list(:,2)==3), :, 1), 'descend');
% Keep the 5 highest only
val = val(1:5);
idx = idx(1:5);
b = bar(1:length(idx), val); b.FaceColor = [0.75 0 0]; b.EdgeColor = 'none';
set(gca, 'XTick', 1:length(idx), 'XTickLabel', motive_labels(idx))
xtickangle(90)
xlim([0 length(idx)+1])
xlabel('motive change')
ylabel('proportion')
title('S2 -> S3')
xlim([0 length(idx)+1])
ylim([0 0.6])

% --- Solitarious ---
subplot(3,3,4)
% Show heatmap containing all proportions of changing motives for all
% stimulus changes
imagesc(map(:, :, 2))
% Cosmetics
axis equal tight
clim([0, 0.45]); colorbar
colormap(gca, SET.HeatmapColors.solitarious)
xticks(1:size(motive_change_list,1))
yticks(1:size(stim_change_list,1))
set(gca, 'XTickLabel', motive_labels, 'YTickLabel', stim_labels)
xtickangle(90)
xlabel('motive change')
ylabel('stim. change')
title('solitarious')

% Show the 5 most prominent motive changes when changing
% from Laa to LaaLct
subplot(3,3,5); hold on
% Sort proportions
[val, idx] = sort(map(find(stim_change_list(:,1)==1 & stim_change_list(:,2)==3), :, 2), 'descend');
% Keep the 5 highest only
val = val(1:5);
idx = idx(1:5);
b = bar(1:length(idx), val); b.FaceColor = [0 0 0.75]; b.EdgeColor = 'none';
set(gca, 'XTick', 1:length(idx), 'XTickLabel', motive_labels(idx))
xtickangle(90)
xlabel('motive change')
ylabel('proportion')
title('S1 -> S3')
xlim([0 length(idx)+1])
ylim([0 0.6])

% Show the 5 most prominent motive changes when changing
% from Lct to LaaLct
subplot(3,3,6); hold on
% Sort proportions
[val, idx] = sort(map(find(stim_change_list(:,1)==2 & stim_change_list(:,2)==3), :, 2), 'descend');
% Keep the 5 highest only
val = val(1:5);
idx = idx(1:5);
b = bar(1:length(idx), val); b.FaceColor = [0 0 0.75]; b.EdgeColor = 'none';
set(gca, 'XTick', 1:length(idx), 'XTickLabel', motive_labels(idx))
xtickangle(90)
xlim([0 length(idx)+1])
xlabel('motive change')
ylabel('proportion')
title('S2 -> S3')
xlim([0 length(idx)+1])
ylim([0 0.6])

% --- Difference greg-soli ---
subplot(3,3,7)
% Show heatmap containing all proportions of changing motives for all
% stimulus changes
imagesc(map(:, :, 3))
% Cosmetics
axis equal tight
clim([min(temp_col(:)) max(temp_col(:))]); colorbar
colormap(gca, cols)
xticks(1:size(motive_change_list,1))
yticks(1:size(stim_change_list,1))
set(gca, 'XTickLabel', motive_labels, 'YTickLabel', stim_labels)
xtickangle(90)
xlabel('motive change')
ylabel('stim. change')
title('difference (greg-soli)')

% Show the 5 most prominent differences when changing
% from Laa to LaaLct
subplot(3,3,8); hold on
% Sort proportions
[val, idx] = sort(abs(map(find(stim_change_list(:,1)==1 & stim_change_list(:,2)==3), :, 3)), 'descend');
% Keep the 5 highest only
val = val(1:5);
idx = idx(1:5);
val = val.*sign(map(find(stim_change_list(:,1)==1 & stim_change_list(:,2)==3), idx, 3));
b = bar(1:length(idx), val); b.FaceColor = [0.251 0.251 0.251]; b.EdgeColor = 'none';
set(gca, 'XTick', 1:length(idx), 'XTickLabel', motive_labels(idx))
xtickangle(90)
xlabel('motive change')
ylabel('d proportion')
title('S1 -> S3')
xlim([0 length(idx)+1])
ylim([-0.2 0.2])

% Show the 5 most prominent differences when changing
% from Lct to LaaLct
subplot(3,3,9); hold on
% Sort proportions
[val, idx] = sort(abs(map(find(stim_change_list(:,1)==2 & stim_change_list(:,2)==3), :, 3)), 'descend');
% Keep the 5 highest only
val = val(1:5);
idx = idx(1:5);
val = val.*sign(map(find(stim_change_list(:,1)==1 & stim_change_list(:,2)==3), idx, 3));
b = bar(1:length(idx), val); b.FaceColor = [0.251 0.251 0.251]; b.EdgeColor = 'none';
set(gca, 'XTick', 1:length(idx), 'XTickLabel', motive_labels(idx))
xtickangle(90)
xlabel('motive change')
ylabel('d proportion')
title('S2 -> S3')
xlim([0 length(idx)+1])
ylim([-0.2 0.2])

% Save and close figure
export_fig('PureOdorAnalysis_OdorSpace\prop_matrix', '-pdf', '-painters')
close(hFig)

save('PureOdorAnalysis_OdorSpace\OdorSpace_kmeans_PoolData.mat', 'PoolData', 'SET')
clearvars -except PoolData SET


%% Visualize the results.
% Show how how the proportion of response motives changes between stimuli

% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Get avg values for each ROI (row) and each stimulus (column)
    data = PoolData.kmeans_result.AverageValue.dat_AvgVal(PoolData.kmeans_result.AverageValue.info(:,1) == iPhase, :);
    % Get each ROI's response motive ID for each stimulus and the
    % corresponding grouping factors (response motives)
    motives = PoolData.kmeans_result.ResponseProfile.ID(PoolData.kmeans_result.ResponseProfile.info(:,1) == iPhase,:);
    stims = PoolData.kmeans_result.ResponseProfile.info(PoolData.kmeans_result.ResponseProfile.info(:,1) == iPhase,3);
    data_grouping = [];
    % Iterate ove all stimuli and populate the grouping variable
    for iStim = 1:length(SET.StimList)
        idx = find(stims == iStim);
        data_grouping = [data_grouping, motives(idx)];
    end
    % Get info on animals
    info = PoolData.(SET.phases{iPhase}).info(:,1);
    % Kick out the first motive as this mostly contains noise
    % info(find(any(data_grouping==1,2)), :) = [];
    % data_grouping(find(any(data_grouping==1,2)), :) = [];
    % Get list of motives
    id_list = unique(data_grouping(:));
    % Get the average proportion of motive IDs per stimulus across all
    % animals
    prop_data{iPhase} = zeros(PoolData.kmeans_result.ResponseProfile.best_k-1, size(data,2));
    % Get list of animals
    ani_list = unique(PoolData.(SET.phases{iPhase}).info(:,1));
    % Iterate over all stimui
    for iStim = 1:length(SET.StimList)
        % Iterate over all motives
        for iK = 1:length(id_list)
            % Iterate over all animals and get the proportion of
            % motives for the current stimulus
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
    Widefield_SubFcn.HivePropPlot(prop_data{iPhase}(:,1:3), 'connected', cols(2:end,:), [], 0.5, 0.05)
    % Add a colorbar for the different levels
    colormap(cols(2:end,:)); colorbar; clim([1 PoolData.kmeans_result.ResponseProfile.best_k])
    % Add labels
    set(gca, 'ThetaTickLabel', SET.StimList, 'TickLabelInterpreter', 'none')
    % Save
    export_fig(['PureOdorAnalysis_OdorSpace\response_profile_prop_hive_',SET.phases{iPhase}], '-pdf', '-painters')
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
Widefield_SubFcn.HiveCircPropPlot(currData(:,1:3), cols, [], 0.5, 0.25)
% Save
export_fig('PureOdorAnalysis_OdorSpace\response_profile_prop_hive_diff', '-pdf', '-painters')
close(hFig)

clearvars -except PoolData SET


%% Visualize the results.
% Show the response vector locations within the AL

% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % One figure per phase
    hFig = figure('units','normalized','Position',[0 0 1 1], 'color', 'w');
    % Iterte over all animals
    ani_list = [PoolData.(SET.phases{iPhase}).info_roi_mask{:,1}];
    for iAni = 1:length(ani_list)
        % Get all you need
        idx_ani = find(PoolData.kmeans_result.AverageValue.info(:,1)==iPhase & abs(PoolData.kmeans_result.AverageValue.info(:,2))==ani_list(iAni));
        pocket_list = PoolData.(SET.phases{iPhase}).info_roi_mask{iAni,3};
        [~,~,id_list] = unique(PoolData.kmeans_result.AverageValue.kmeans_ID(idx_ani));
        pockets_labeled = PoolData.(SET.phases{iPhase}).info_roi_mask{iAni,4};
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
        clim([1, PoolData.kmeans_result.AverageValue.kmeans_best_k])
        colorbar
        cols = SET.HeatmapColors.(SET.phases{iPhase});
        cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), PoolData.kmeans_result.AverageValue.kmeans_best_k+2)),:);
        cols = cols(2:end-1,:);
        colormap(cols)
        axis equal tight
        xticks([])
        yticks([])
        box on
        title(PoolData.(SET.phases{iPhase}).info_roi_mask{iAni,2}, 'interpreter', 'none')
    end%iAni
    print(['PureOdorAnalysis_OdorSpace\resp_vec_loc_',SET.phases{iPhase}], '-dpdf')
    close(hFig)
end%iPhase
clearvars -except PoolData SET
