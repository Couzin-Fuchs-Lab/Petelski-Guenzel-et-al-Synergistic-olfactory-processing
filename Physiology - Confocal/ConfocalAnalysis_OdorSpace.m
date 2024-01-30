%% ConfocalAnalysis_OdorSpace
% This script applies a deeper analysis of the population data extracted by
% ConfocalAnalysis_Population. The three major aspects are (i) the
% clustering of odor response vectors, (ii) the clustering of odor response
% motifs, (iii) and an analysis of how informative the dynamics of motif
% transitions are for predicting the animal's social phenotype. These
% analyses are accompanied by their visualization.
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


%% Settings

% Set paths
SET.path2data = '...\PooledData.mat';

% Give list of phases
SET.phases = {'gregarious'; 'solitarious'};

% Give layer names
SET.layers = {'layer01_top', 'layer02_middle'};

% Give field names for the data required for the response vectors and
% motifs, respectively
SET.FieldName_RespVecData = 'StimCombiList6';
SET.FieldName_MotifData = 'StimCombiList1';

% Provide stimulus names for reponse vectors and motifs
SET.StimList_RespVec =  {'COL', 'ZHAE', 'ZHAECOL', 'Z3HL', 'Z3HLCOL', 'OCT', 'OCTCOL', 'Z3HLZHAE', 'OCTZHAE'};
SET.StimList_Motif = {'COL', 'ZHAE', 'ZHAECOL'};

% Set the reference stimulus that should be excluded
SET.RefStim = 'MOL';

% Set framerate
SET.fps = 1;
% Set the number of frames the final trace should have
SET.totalFrameNumber = 17;
SET.time_vec = 0:(1/SET.fps):((SET.totalFrameNumber/SET.fps)-(1/SET.fps));
% Set the analysis window. Note this is hard-coded previous knowledge for
% data that has been cut.
SET.AnalysisWindow = [...
    10 11;... layer 01
    8 9;...   layer 02
    ];

% Clusterng & Dimensionality reduction
SET.kmeans.MaxIter = 1e6;
SET.kmeans.OnlinePhase = 'on';
SET.kmeans.Replicates = 5e2;
SET.maxk_eval = 20;
SET.kmeans.dist_RespProfile = 'sqeuclidean';
SET.kmeans.dist_AvgValue = 'sqeuclidean';
SET.pca_explained = 99;

% Set how many of the most frequent motif-based Bliss score time courses
% should be show
SET.Num_most_freq_BS = 12;


% Cosmetics
SET.Colors = [...
    191,191,191;... MOL

    217,095,002;... COL
    102,166,030;... ZHAE
    231,041,138;... ZHAECOL

    166,168,016;... Z3HL
    231,041,138;... Z3HLCOL

    230,171,002;... OCT
    117,011,179;... OCTCOL

    027,158,119;... Z3HLZHAE
    167,103,090;... OCTZHAE
    ]/255;

SET.ColorNames = {...
    'MOL';...
    'COL';...
    'ZHAE';...
    'ZHAECOL';...
    'Z3HL';...
    'Z3HLCOL';...
    'OCT';...
    'OCTCOL';...
    'Z3HLZHAE';...
    'OCTZHAE'};

% Colormap for heatmaps
SET.HeatmapColors.gregarious = Confocal_SubFcn.ColMapPlasma(1000);
SET.HeatmapColors.solitarious = [...
    SET.HeatmapColors.gregarious(:,2),...
    SET.HeatmapColors.gregarious(:,1),...
    SET.HeatmapColors.gregarious(:,3)];
SET.HeatmapColors.other = gray(1000);

% Load data
PooledData = load(SET.path2data);
PooledData = PooledData.PooledData;

%% Pool phases and get response motif data

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % Prepare folder structure
    mkdir(['ConfocalAnalysis_OdorSpace\', (SET.layers{iLayer})])
    % Prepare variables
    % --- Response motifs based on avg time courses
    OdorSpace.(SET.layers{iLayer}).Motif.data = [];
    OdorSpace.(SET.layers{iLayer}).Motif.info = [];
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % Motif data (exculde the reference stimulus)
        ind = find(~strcmp(PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.FieldName_MotifData).poolInfo(:,4), SET.RefStim));
        OdorSpace.(SET.layers{iLayer}).Motif.data = [OdorSpace.(SET.layers{iLayer}).Motif.data;...
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.FieldName_MotifData).poolTC(ind,:)];
        % Info on motif data [Phase, Animal, Pocket, Stimulus]
        OdorSpace.(SET.layers{iLayer}).Motif.info = [OdorSpace.(SET.layers{iLayer}).Motif.info;...
            ones(size(cell2mat(PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.FieldName_MotifData).poolInfo(ind,1:2)),1),1)*iPhase,...
            cell2mat(PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.FieldName_MotifData).poolInfo(ind,1:3))];
        OdorSpace.(SET.layers{iLayer}).Motif.info(:,2) = OdorSpace.(SET.layers{iLayer}).Motif.info(:,2)*((iPhase-1.5)*2);
    end%iPhase
    % Keep original versions (not normalized) for later
    OdorSpace.(SET.layers{iLayer}).Motif.data_original = OdorSpace.(SET.layers{iLayer}).Motif.data;
end%iLayer
% Save data
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')


%% Pool phases and get response vector data

% Do this for the top layer only
iLayer = 1;

% Prepare variables
% --- Response vector based on avg amplitude during active window
OdorSpace.(SET.layers{iLayer}).RespVec.data = [];
OdorSpace.(SET.layers{iLayer}).RespVec.info = [];
% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Reshape response vector data (exculde the reference stimulus)
    ind = find(~strcmp(PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.FieldName_RespVecData).poolInfo(:,4), SET.RefStim));
    temp_resp_vec = PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.FieldName_RespVecData).poolAvg(ind,:);
    temp_info  = cell2mat(PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.FieldName_RespVecData).poolInfo(ind,1:3));
    temp_ani_list = unique(temp_info(:,1));
    temp_stim_list = unique(temp_info(:,3));
    temp_data_shaped = [];
    temp_info_shaped = [];
    % Iterate over all animals
    for iAni = 1:length(temp_ani_list)
        v = nan( length(unique(temp_info(temp_info(:,1)==temp_ani_list(iAni),2)) ), length(temp_stim_list));
        for iStim = 1:length(temp_stim_list)
            ind = find(temp_info(:,1)==temp_ani_list(iAni) & temp_info(:,3)==temp_stim_list(iStim));
            v(:,iStim) = temp_resp_vec(ind);
        end%iStim
        temp_data_shaped = [temp_data_shaped; v];
        temp_info_shaped = [temp_info_shaped; ones(length(ind),1)*iPhase, ones(length(ind),1)*temp_ani_list(iAni)*((iPhase-1.5)*2), temp_info(ind,2)];
    end%iAni
    % Response vector data
    OdorSpace.(SET.layers{iLayer}).RespVec.data = [OdorSpace.(SET.layers{iLayer}).RespVec.data;...
        temp_data_shaped];
    % Info on response vectors [Phase, Animal, Pocket]
    OdorSpace.(SET.layers{iLayer}).RespVec.info = [OdorSpace.(SET.layers{iLayer}).RespVec.info;...
        temp_info_shaped];
    clear temp*
end%iPhase
% Keep original versions (not normalized) for later
OdorSpace.(SET.layers{iLayer}).RespVec.data_original = OdorSpace.(SET.layers{iLayer}).RespVec.data;
% Save data and clear workspace a bit
clear PooledData
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')


%% Normalize data from each animal

% Iterate over both layers
for iLayer = 1:length(SET.layers)
    % --- response motif
    ani_list = unique(OdorSpace.(SET.layers{iLayer}).Motif.info( OdorSpace.(SET.layers{iLayer}).Motif.info(:,1)==iPhase, 2 ));
    for iAni = 1:length(ani_list)
        idx_k = find(OdorSpace.(SET.layers{iLayer}).Motif.info(:, 2) == ani_list(iAni));
        curr_dat = OdorSpace.(SET.layers{iLayer}).Motif.data(idx_k,:);
        curr_dat = curr_dat - min(curr_dat(:));
        curr_dat = sqrt(curr_dat);
        curr_dat = curr_dat - nanmean(curr_dat(:));
        curr_dat = curr_dat / nanstd(curr_dat(:));
        OdorSpace.(SET.layers{iLayer}).Motif.data(idx_k,:) = curr_dat;
    end%iAni
    % --- response vector
    if iLayer == 1
        ani_list = unique(OdorSpace.(SET.layers{iLayer}).RespVec.info(OdorSpace.(SET.layers{iLayer}).RespVec.info(:,1)==iPhase, 2));
        for iAni = 1:length(ani_list)
            idx_k = find(OdorSpace.(SET.layers{iLayer}).RespVec.info(:, 2) == ani_list(iAni));
            curr_dat = OdorSpace.(SET.layers{iLayer}).RespVec.data(idx_k,:);
            curr_dat = curr_dat - min(curr_dat(:));
            curr_dat = sqrt(curr_dat);
            curr_dat = curr_dat - nanmean(curr_dat(:));
            curr_dat = curr_dat / nanstd(curr_dat(:));
            OdorSpace.(SET.layers{iLayer}).RespVec.data(idx_k,:) = curr_dat;
        end%iAni
    end
end%iLayer
% Save
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')


%% Cluster response motifs
% Here we cluster all time courses from all ROIs responsing to all stimuli
% to check whether there are unique response motifs (e.g., slow rise,
% long peak, inhibition, etc)

% Iterate over all layers
for iLayer = 1:length(SET.layers)

    % Normalize each row
    dat = normalize(OdorSpace.(SET.layers{iLayer}).Motif.data, 2);

    % PCA preprocessing
    opts = statset('pca'); opts.MaxIter = 1e7; opts.TolFun = 1e-8; opts.TolFunX = 1e-8;
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
            'MaxIter', ceil(1+SET.kmeans.MaxIter/10),...
            'OnlinePhase', SET.kmeans.OnlinePhase,...
            'Replicates', ceil(1+SET.kmeans.Replicates/10),...
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

    % Sort centroids manually
    % Depict all clusters
    hFig1 = figure('units', 'normalized', 'Position', [0 0 0.5 1], 'Color', 'w');
    hFig2 = figure('units', 'normalized', 'Position', [0.5 0 0.5 1], 'Color', 'w');
    for iC = 1:kmeans_result.best_k
        idx = find(ID==iC);
        figure(hFig1)
        nexttile
        plot(mean(OdorSpace.(SET.layers{iLayer}).Motif.data_original(idx,:)), 'LineWidth',3)
        title(iC)
        figure(hFig2); hold on
        plot(mean(OdorSpace.(SET.layers{iLayer}).Motif.data_original(idx,:)), 'LineWidth',3)
    end
    figure(hFig2); legend;
    % Ask user
    clc
    disp('Give response motif cluster order as, e.g., [1 3 2]')
    order_idx = input('Order: ');
    % --- Sort centroids
    temp_ID = nan(size(ID));
    for iC = 1:kmeans_result.best_k
        temp_ID(ID==order_idx(iC)) = iC;
    end%iK
    ID = temp_ID;
    close(hFig1); close(hFig2)

    % Run a LDA on the input data for later visualization
    dat = normalize(OdorSpace.(SET.layers{iLayer}).Motif.data, 2);
    Mdl = fitcdiscr(dat,ID);
    [W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma); %Must be in the right order!
    lambda = diag(LAMBDA);
    [lambda, SortOrder] = sort(lambda, 'descend');
    W = W(:, SortOrder);
    dat_lda = dat*W;

    % Get list of IDs
    ID_list = unique(ID);
    % Iterate over all phases
    for iPhase = 1:length(SET.phases)
        % Get list of all animals for the current phase
        ani_list = unique(OdorSpace.(SET.layers{iLayer}).Motif.info( OdorSpace.(SET.layers{iLayer}).Motif.info(:,1)==iPhase, 2 ));
        % Preallocation
        kmeans_result.centroid_prop{iPhase} = zeros(length(ani_list), length(ID_list));
        % Iterate over all clusters
        for iC = 1:length(ID_list)
            % Preallocation
            kmeans_result.centroid{iPhase}{iC} = nan(length(ani_list), size(OdorSpace.(SET.layers{iLayer}).Motif.data_original,2));
            % Iterate over all animals, get their centroid and store it
            for iAni = 1:length(ani_list)
                % Get indices for the current combination of phase, animal,
                % and cluster
                ind = find(...
                    OdorSpace.(SET.layers{iLayer}).Motif.info(:,1)==iPhase & ...
                    OdorSpace.(SET.layers{iLayer}).Motif.info(:,2)==ani_list(iAni) & ...
                    ID == ID_list(iC));
                kmeans_result.centroid{iPhase}{iC}(iAni,:) = nanmean(OdorSpace.(SET.layers{iLayer}).Motif.data_original(ind,:),1);
                kmeans_result.centroid_lda{iPhase}{iC}(iAni,:) = nanmean(dat_lda(ind,:),1);
                kmeans_result.centroid_prop{iPhase}(iAni,iC) = length(ind);
            end%iAni
        end%iC
        kmeans_result.centroid_prop{iPhase} = kmeans_result.centroid_prop{iPhase}./sum(kmeans_result.centroid_prop{iPhase},2);
    end%iPhase

    % Save
    OdorSpace.(SET.layers{iLayer}).Motif.best_k = kmeans_result.best_k;
    OdorSpace.(SET.layers{iLayer}).Motif.ID = ID;
    OdorSpace.(SET.layers{iLayer}).Motif.centroid = kmeans_result.centroid;
    OdorSpace.(SET.layers{iLayer}).Motif.centroid_lda = kmeans_result.centroid_lda;
    OdorSpace.(SET.layers{iLayer}).Motif.centroid_prop = kmeans_result.centroid_prop;

end%iLayer

% Save
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')
% Clean up
clearvars -except OdorSpace SET


%% Cluster Response Vectors
% Here, we cluster ROIs based on their response magnitude to each stimulus.

% Do this for Layer 1 only
iLayer = 1;

% Normalize each row
dat = normalize(OdorSpace.(SET.layers{iLayer}).RespVec.data, 2, 'range');
dat = (dat-0.5)*2;

% PCA preprocessing
opts = statset('pca'); opts.MaxIter = 1e6; opts.TolFun = 1e-8; opts.TolFunX = 1e-8;
[~,dat_pca,~,~,explained,~] = pca(dat, 'Economy', false, 'Options', opts);
ind = find(cumsum(explained)>SET.pca_explained, 1);
dat_pca = dat_pca(:,1:ind);

% Run multiple times to identify the optimal number of clusters
VRC = nan(SET.maxk_eval,1);
hWait = waitbar(1, 'Evaluating clusters ...');
% Test several cluster sizes
for k = 2:SET.maxk_eval
    waitbar(k/SET.maxk_eval, hWait)
    rng(1234)
    [idx,C,sumd] = kmeans(dat_pca, k,...
        "Distance", SET.kmeans.dist_AvgValue,...
        'MaxIter', ceil(1+SET.kmeans.MaxIter/10),...
        'OnlinePhase', SET.kmeans.OnlinePhase,...
        'Replicates', ceil(1+SET.kmeans.Replicates/10),...
        'Options', statset('UseParallel',1));
    % Get the variance ratio criterion (VRC)
    VRC(k) = Confocal_SubFcn.VarianceRatioCriterion(dat_pca, k, idx, C);
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
[ID, C] = kmeans(dat_pca, kmeans_result.best_k,...
    "Distance", SET.kmeans.dist_AvgValue,...
    'MaxIter', SET.kmeans.MaxIter,...
    'OnlinePhase', SET.kmeans.OnlinePhase,...
    'Replicates', SET.kmeans.Replicates,...
    'Options', statset('UseParallel',1));
close(hWait)

% Run a LDA on the input data for later visualization
dat = normalize(OdorSpace.(SET.layers{iLayer}).RespVec.data, 2, 'range');
dat = (dat-0.5)*2;
Mdl = fitcdiscr(dat,ID);
[W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma); %Must be in the right order!
lambda = diag(LAMBDA);
[lambda, SortOrder] = sort(lambda, 'descend');
W = W(:, SortOrder);
dat_lda = dat*W;

% Get list of IDs
ID_list = unique(ID);
% Iterate over all phases
for iPhase = 1:length(SET.phases)
    % Get list of all animals for the current phase
    ani_list = unique(OdorSpace.(SET.layers{iLayer}).RespVec.info( OdorSpace.(SET.layers{iLayer}).RespVec.info(:,1)==iPhase, 2 ));
    % Preallocation
    kmeans_result.centroid_prop{iPhase} = zeros(length(ani_list), length(ID_list));
    % Iterate over all clusters
    for iC = 1:length(ID_list)
        % Preallocation
        kmeans_result.centroid_orig{iPhase}{iC} = nan(length(ani_list), size(OdorSpace.(SET.layers{iLayer}).RespVec.data_original,2));
        kmeans_result.centroid_norm{iPhase}{iC} = nan(length(ani_list), size(OdorSpace.(SET.layers{iLayer}).RespVec.data_original,2));
        kmeans_result.centroid_norm_pca{iPhase}{iC} = nan(length(ani_list), size(OdorSpace.(SET.layers{iLayer}).RespVec.data_original,2));
        kmeans_result.centroid_lda{iPhase}{iC} = nan(length(ani_list), size(OdorSpace.(SET.layers{iLayer}).RespVec.data_original,2));
        % Iterate over all animals, get their centroid and store it
        for iAni = 1:length(ani_list)
            % Get indices for the current combination of phase, animal,
            % and cluster
            ind = find(...
                OdorSpace.(SET.layers{iLayer}).RespVec.info(:,1)==iPhase & ...
                OdorSpace.(SET.layers{iLayer}).RespVec.info(:,2)==ani_list(iAni) & ...
                ID == ID_list(iC));
            kmeans_result.centroid_orig{iPhase}{iC}(iAni,:) = nanmean(OdorSpace.(SET.layers{iLayer}).RespVec.data_original(ind,:),1);
            kmeans_result.centroid_norm{iPhase}{iC}(iAni,:) = nanmean(dat(ind,:),1);
            kmeans_result.centroid_norm_pca{iPhase}{iC}(iAni,:) = nanmean(dat_pca(ind,:),1);
            kmeans_result.centroid_lda{iPhase}{iC}(iAni,:) = nanmean(dat_lda(ind,:),1);
            kmeans_result.centroid_prop{iPhase}(iAni,iC) = length(ind);
        end%iAni
    end%iC
    kmeans_result.centroid_prop{iPhase} = kmeans_result.centroid_prop{iPhase}./sum(kmeans_result.centroid_prop{iPhase},2);
end%iPhase

% Save
OdorSpace.(SET.layers{iLayer}).RespVec.best_k = kmeans_result.best_k;
OdorSpace.(SET.layers{iLayer}).RespVec.ID = ID;
OdorSpace.(SET.layers{iLayer}).RespVec.centroid_orig = kmeans_result.centroid_orig;
OdorSpace.(SET.layers{iLayer}).RespVec.centroid_norm = kmeans_result.centroid_norm;
OdorSpace.(SET.layers{iLayer}).RespVec.centroid_norm_pca = kmeans_result.centroid_norm_pca;
OdorSpace.(SET.layers{iLayer}).RespVec.centroid_lda = kmeans_result.centroid_lda;
OdorSpace.(SET.layers{iLayer}).RespVec.centroid_prop = kmeans_result.centroid_prop;

% Save
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')
% Clean up
clearvars -except OdorSpace SET


%% Visualize the results.
% response motif time courses

% Define the range of the y-axis
ylim_range = [...
    -5 10;... layer 01
    -2 8.0... layer 02
    ];

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
        cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), OdorSpace.(SET.layers{iLayer}).Motif.best_k+2)),:);
        cols = cols(2:end-1,:);
        % Get position of active window
        pos = [SET.time_vec(SET.AnalysisWindow(iLayer,1)), ylim_range(iLayer,1), SET.time_vec(SET.AnalysisWindow(iLayer,end))-SET.time_vec(SET.AnalysisWindow(iLayer,1)), diff(ylim_range(iLayer,:))];
        % Plot active window as gray box
        rectangle('position', pos, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
        % Iterate over all motifs
        for iC = 1:OdorSpace.(SET.layers{iLayer}).Motif.best_k
            % Get the average
            avg = nanmean(OdorSpace.(SET.layers{iLayer}).Motif.centroid{iPhase}{iC},1);
            % Get the corresponding confidence interval
            ci = bootci(5000, {@nanmean, OdorSpace.(SET.layers{iLayer}).Motif.centroid{iPhase}{iC}});
            % Plot both. Turn confidence intervals into shaded areas during
            % figure post-processing
            plot(SET.time_vec, avg, 'LineWidth', 2, 'Color', cols(iC,:))
            plot(SET.time_vec, ci', 'LineWidth', 1, 'Color', cols(iC,:))
            drawnow
        end%iK
        % Cosmetics
        ylabel('\DeltaF/F (%)')
        xlabel('time (s)')
        title(SET.phases{iPhase})
        xlim([SET.time_vec(1), SET.time_vec(end)])
        xticks(SET.time_vec)
        ylim(ylim_range(iLayer,:))
        yticks(ylim_range(iLayer,1):ylim_range(iLayer,2))
        axis square
    end%iPhase
    %Also, plot the motifs for both phases combined in the third subplot
    subplot(1,3,3) ; hold on
    % Use a different colormap
    cols = Confocal_SubFcn.ColMapInferno(1000);
    cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), OdorSpace.(SET.layers{iLayer}).Motif.best_k+2)),:);
    cols = cols(2:end-1,:);
    % Get the position of the active region
    pos = [SET.time_vec(SET.AnalysisWindow(iLayer,1)), ylim_range(iLayer,1), SET.time_vec(SET.AnalysisWindow(iLayer,end))-SET.time_vec(SET.AnalysisWindow(iLayer,1)), diff(ylim_range(iLayer,:))];
    % Plot the active region as a gray box
    rectangle('position', pos, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
    % Iterate over all motifs
    for iC = 1:OdorSpace.(SET.layers{iLayer}).Motif.best_k
        % Pool gregarious and solitarious data
        currDat = [];
        for iPhase = 1:length(SET.phases)
            currDat = [currDat; OdorSpace.(SET.layers{iLayer}).Motif.centroid{iPhase}{iC}];
        end%iPhase
        % Get the average
        avg = nanmean(currDat,1);
        % And the corresponding confidence interval
        ci = bootci(5000, {@nanmean, currDat});
        % Plot both. Turn confidence intervals into shaded areas during
        % figure post-processing
        plot(SET.time_vec, avg, 'LineWidth', 2, 'Color', cols(iC,:))
        plot(SET.time_vec, ci', 'LineWidth', 1, 'Color', cols(iC,:))
        drawnow
    end%iK
    % Cosmetics
    ylabel('\DeltaF/F (%)')
    xlabel('time (s)')
    title(SET.phases{iPhase})
    xlim([SET.time_vec(1), SET.time_vec(end)])
    xticks(SET.time_vec)
    ylim(ylim_range(iLayer,:))
    yticks(ylim_range(iLayer,1):ylim_range(iLayer,2))
    title('both')
    axis square
    % Save figure and close it
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\kmeans_ResponseProfile'], '-pdf', '-painters')
    close(hFig)
end%iLayer

% Save
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')
% Clean up
clearvars -except OdorSpace SET


%% Visualize the results.
% Motif combinations and how they relate to synergy

curr_ylim_avg{1} = [-0.5 1.5];
curr_ylim_avg{2} = [-1 2];
curr_ylim_single{1} = [-12 12];
curr_ylim_single{2} = [-12 12];
% Iterate over all layers
for iLayer = 1:length(SET.layers)

    % Iterate over both phases and get all motif combinations
    for iPhase = 1:length(SET.phases)
        % Get corresponding grouping factors (response motifs)
        motifs = OdorSpace.(SET.layers{iLayer}).Motif.ID(OdorSpace.(SET.layers{iLayer}).Motif.info(:,1) == iPhase, :);
        stims = OdorSpace.(SET.layers{iLayer}).Motif.info(OdorSpace.(SET.layers{iLayer}).Motif.info(:,1) == iPhase, 4);
        anis = OdorSpace.(SET.layers{iLayer}).Motif.info(OdorSpace.(SET.layers{iLayer}).Motif.info(:,1) == iPhase, 2);
        data_grouping.(SET.phases{iPhase}) = [];
        data_grouping_info.(SET.phases{iPhase}) = [];
        stim_list = unique(stims);
        % Iterate ove all stimuli and populate the grouping variable
        for iStim = 1:length(stim_list)
            idx = find(stims == stim_list(iStim));
            data_grouping.(SET.phases{iPhase}) = [data_grouping.(SET.phases{iPhase}), motifs(idx)];
            data_grouping_info.(SET.phases{iPhase}) = [data_grouping_info.(SET.phases{iPhase}), anis(idx)];
        end%iStim
        data_grouping_info.(SET.phases{iPhase}) = data_grouping_info.(SET.phases{iPhase})(:,1);
        % Save for later
        OdorSpace.(SET.layers{iLayer}).Motif.grouping.(SET.phases{iPhase}) = data_grouping.(SET.phases{iPhase});
        OdorSpace.(SET.layers{iLayer}).Motif.grouping_info.(SET.phases{iPhase}) = data_grouping_info.(SET.phases{iPhase});
        % For now, kick out regions that always had the same motif
        data_grouping_info.(SET.phases{iPhase})(std(data_grouping.(SET.phases{iPhase}), [], 2) == 0, :) = [];
        data_grouping.(SET.phases{iPhase})(std(data_grouping.(SET.phases{iPhase}), [], 2) == 0, :) = [];
    end%iPhase

    % Get all motifs as grand mean over all animals of both phases
    for iC = 1:OdorSpace.(SET.layers{iLayer}).Motif.best_k
        avg = [];
        for iPhase = 1:length(SET.phases)
            avg = [avg; OdorSpace.layer01_top.Motif.centroid{iPhase}{iC}];
        end%iPhase
        motif.(['m', num2str(iC)]) = nanmean(avg,1)/100;
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
        % Average by animal
        ani_list = unique(data_grouping_info.(SET.phases{iPhase}));
        syn_avg.(SET.phases{iPhase}) = nan(length(ani_list), length(SET.time_vec));
        for iAni = 1:length(ani_list)
            ind = find(data_grouping_info.(SET.phases{iPhase}) == ani_list(iAni));
            syn_avg.(SET.phases{iPhase})(iAni,:) = nanmean(syn.(SET.phases{iPhase})(ind,:), 1);
        end%iAni
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
        grouping_count(iGrp,3) = (grouping_count(iGrp,1) - grouping_count(iGrp,2));
    end%iGrp

    % Get the most substantial differences and their corresponding
    [sort_grouping_diff, sort_grouping_ind] = sort(grouping_count(:,3), 'descend', 'ComparisonMethod', 'abs');

    % Depict results - overall avg Bliss score
    hFig = figure('units', 'normalized', 'Position', [0 0 1 1], 'Color', 'w');
    % Overall Bliss score
    nexttile; hold on
    cols = [191 0 0; 0 0 191]/255;
    for iPhase = 1:length(SET.phases)
        avg = nanmean(syn.(SET.phases{iPhase}),1);
        ci = bootci(5000, {@nanmean, syn.(SET.phases{iPhase})});
        plot(SET.time_vec, avg, 'LineWidth',3, 'Color', cols(iPhase,:))
        plot(SET.time_vec, ci(1,:), 'LineWidth',1, 'Color', cols(iPhase,:))
        plot(SET.time_vec, ci(2,:), 'LineWidth',1, 'Color', cols(iPhase,:))
    end%iPhase
    xlim([SET.time_vec(1), SET.time_vec(end)])
    ylim(curr_ylim_avg{iLayer})
    xlabel('time (s)')
    ylabel('Bliss score')
    axis square
    % Depict those with the largest difference between soli and greg
    motif_groups = grouping_unique(sort_grouping_ind,:);
    for iM = 1:SET.Num_most_freq_BS
        % Calculate synergy
        a = motif.(['m', num2str(motif_groups(iM,1))]);
        b = motif.(['m', num2str(motif_groups(iM,2))]);
        c = motif.(['m', num2str(motif_groups(iM,3))]);
        % Get expected and observed
        f_expected = (a/2 + b/2) - ((a/2).*(b/2));
        f_observed = c;
        % Get BS
        curr_syn = 100*(f_observed - f_expected);
        % Plot
        nexttile; hold on
        plot(SET.time_vec, curr_syn, 'Color', [sort_grouping_diff(iM)>=0, 0, sort_grouping_diff(iM)<=0], 'LineWidth', 2)
        xlim([SET.time_vec(1), SET.time_vec(end)])
        title([num2str(iM), '|',...
            ' [',...
            num2str(motif_groups(iM,1)),',',...
            num2str(motif_groups(iM,2)),',',...
            num2str(motif_groups(iM,3)),']'])
        ylim(curr_ylim_single{iLayer})
        xlim([SET.time_vec(1), SET.time_vec(end)])
        axis square
    end%iM
    % Save
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\motif_combinations_BS'], '-pdf', '-painters')
    close(hFig)

end%iLayer

% Save
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')
% Clean up
clearvars -except OdorSpace SET


%% Visualize the results.
% Animals per response motif

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % One figure per layer
    hFig = figure('Units','normalized','Position',[0 0 1 1]);
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % Get list of respective cluster IDs for this phase
        id_list = unique( OdorSpace.(SET.layers{iLayer}).Motif.ID(OdorSpace.(SET.layers{iLayer}).Motif.info(:,1) == iPhase));
        % Get list of all possible animals for this phase
        ani_list = unique(OdorSpace.(SET.layers{iLayer}).Motif.info(OdorSpace.(SET.layers{iLayer}).Motif.info(:,1) == iPhase, 2));
        % Iterate over all clusters
        for iC = 1:length(id_list)
            % Get index positions and observations
            idx = find(OdorSpace.(SET.layers{iLayer}).Motif.ID == id_list(iC));
            anis = OdorSpace.(SET.layers{iLayer}).Motif.info(idx,2);
            % Count how many times each animal was present
            cnts = histcounts(anis, [ani_list(:); max(ani_list)+1]);
            % Normalize by number of ROIs per animal
            cnts = cnts./histcounts(OdorSpace.(SET.layers{iLayer}).Motif.info(:,2), [ani_list(:); max(ani_list)+1]);
            % Get percentage
            cnts = cnts*100;
            % Depict as pie chart
            nexttile
            ax = gca;
            pie(ax, cnts);
            colormap(ax, SET.HeatmapColors.(SET.phases{iPhase}))
            title([SET.phases{iPhase}(1), num2str(iC)])
        end%iC
        % Depict centroids of the cluster
        % First get colors
        cols = SET.HeatmapColors.(SET.phases{iPhase});
        cols = cols(floor(linspace(1,size(SET.HeatmapColors.(SET.phases{iPhase}),1), length(id_list)+2)),:);
        cols = cols(2:end-1,:);
        % Then, iterate over all clusters
        for iC = 1:length(id_list)
            % Depict centroid
            nexttile
            plot(nanmean(OdorSpace.(SET.layers{iLayer}).Motif.centroid{iPhase}{iC}), 'color', cols(iC,:), 'LineWidth', 2)
            title(['cluster ',num2str(iC)])
            axis tight off
        end%iC
    end%iPhase
    % Save
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\pie_RespProfile'], '-pdf', '-painters')
    close(hFig)
end%iLayer

% Save
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')
% Clean up
clearvars -except OdorSpace SET

%% Visualize the results.
% Animals per response vector cluster

% Only for layer 1
iLayer = 1;
% One figure per layer
hFig = figure('Units','normalized','Position',[0 0 1 1]);
% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Get list of respective cluster IDs
    id_list = unique(OdorSpace.(SET.layers{iLayer}).RespVec.ID(OdorSpace.(SET.layers{iLayer}).RespVec.info(:,1) == iPhase));
    % Get list of all possible animals
    ani_list = unique(OdorSpace.(SET.layers{iLayer}).RespVec.info(OdorSpace.(SET.layers{iLayer}).RespVec.info(:,1) == iPhase,2));
    % Iterate over all clusters
    for iC = 1:length(id_list)
        % Get index positions and observations
        idx = find(OdorSpace.(SET.layers{iLayer}).RespVec.ID == id_list(iC));
        anis = OdorSpace.(SET.layers{iLayer}).RespVec.info(idx,2);
        % Count how many times each animal was present
        cnts = histcounts(anis, [ani_list(:); max(ani_list)+1]);
        % Normalize by number of ROIs per animal
        cnts = cnts./histcounts(OdorSpace.(SET.layers{iLayer}).RespVec.info(:,2), [ani_list(:); max(ani_list)+1]);
        % Get percentage
        cnts = cnts*100;
        % Depict as pie chart
        nexttile
        ax = gca;
        pie(ax, cnts);
        colormap(ax, SET.HeatmapColors.(SET.phases{iPhase}))
        title([SET.phases{iPhase}(1), num2str(iC)])
    end%iC
end%iPhase
export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\pie_RespVec'], '-pdf', '-painters')
close(hFig)

% Save
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')
% Clean up
clearvars -except OdorSpace SET

%% Visualize the results.
% Show LDA olfaction space and clustering

% Only layer 1
iLayer = 1;

hFig = figure('Units','normalized','Position',[0 0 1 1]);
% Response motif space
subplot(1,2,1); hold on
% --- Get list of IDs
ID_list = unique(OdorSpace.(SET.layers{iLayer}).Motif.ID);
% Get nice colors
% Use a different colormap
cols = Confocal_SubFcn.ColMapInferno(1000);
cols = cols(floor(linspace(1,size(cols,1), OdorSpace.(SET.layers{iLayer}).Motif.best_k+2)),:);
cols = cols(2:end-1,:);
% Iterate over all clusters
for iC = 1:OdorSpace.(SET.layers{iLayer}).Motif.best_k
    % Draw a convex hull
    X = [OdorSpace.(SET.layers{iLayer}).Motif.centroid_lda{1}{iC}(:,1); OdorSpace.(SET.layers{iLayer}).Motif.centroid_lda{2}{iC}(:,1)];
    Y = [OdorSpace.(SET.layers{iLayer}).Motif.centroid_lda{1}{iC}(:,2); OdorSpace.(SET.layers{iLayer}).Motif.centroid_lda{2}{iC}(:,2)];
    ind = find(any(~isnan([X,Y]),2));
    k = convhull(X(ind), Y(ind));
    plot(X(ind(k)), Y(ind(k)), 'Color', cols(iC,:));
    % Draw each animal's centroid
    plot(OdorSpace.(SET.layers{iLayer}).Motif.centroid_lda{1}{iC}(:,1), OdorSpace.(SET.layers{iLayer}).Motif.centroid_lda{1}{iC}(:,2), 'o', 'MarkerFaceColor',cols(iC,:), 'MarkerEdgeColor','k')
    plot(OdorSpace.(SET.layers{iLayer}).Motif.centroid_lda{2}{iC}(:,1), OdorSpace.(SET.layers{iLayer}).Motif.centroid_lda{2}{iC}(:,2), 'd', 'MarkerFaceColor',cols(iC,:), 'MarkerEdgeColor','k')
end%iV
axis equal tight
title('motif space')
xlabel('LD 1'); ylabel('LD 2'); zlabel('LD 3')
xticks([]); yticks([]); zticks([])

% Response vector space
subplot(1,2,2); hold on
% --- Get list of IDs
ID_list = unique(OdorSpace.(SET.layers{iLayer}).RespVec.ID);
% Get nice colors
% Use a different colormap
cols = Confocal_SubFcn.ColMapInferno(1000);
cols = cols(floor(linspace(1,size(cols,1), OdorSpace.(SET.layers{iLayer}).RespVec.best_k+2)),:);
cols = cols(2:end-1,:);
% Iterate over all clusters
for iC = 1:OdorSpace.(SET.layers{iLayer}).RespVec.best_k
    % Draw a convex hull
    X = [OdorSpace.(SET.layers{iLayer}).RespVec.centroid_lda{1}{iC}(:,1); OdorSpace.(SET.layers{iLayer}).RespVec.centroid_lda{2}{iC}(:,1)];
    Y = [OdorSpace.(SET.layers{iLayer}).RespVec.centroid_lda{1}{iC}(:,2); OdorSpace.(SET.layers{iLayer}).RespVec.centroid_lda{2}{iC}(:,2)];
    ind = find(any(~isnan([X,Y]),2));
    k = convhull(X(ind), Y(ind));
    plot(X(ind(k)), Y(ind(k)), 'Color', cols(iC,:));
    % Draw each animal's centroid
    plot(OdorSpace.(SET.layers{iLayer}).RespVec.centroid_lda{1}{iC}(:,1), OdorSpace.(SET.layers{iLayer}).RespVec.centroid_lda{1}{iC}(:,2), 'o', 'MarkerFaceColor',cols(iC,:), 'MarkerEdgeColor','k')
    plot(OdorSpace.(SET.layers{iLayer}).RespVec.centroid_lda{2}{iC}(:,1), OdorSpace.(SET.layers{iLayer}).RespVec.centroid_lda{2}{iC}(:,2), 'd', 'MarkerFaceColor',cols(iC,:), 'MarkerEdgeColor','k')
end%iV
axis equal tight
title('vector space')
xlabel('LD 1'); ylabel('LD 2'); zlabel('LD 3')
xticks([]); yticks([]); zticks([])

% Export and save
export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\kmeans_LDA_space'], '-pdf', '-painters')
close(hFig)
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')
% Clean up
clearvars -except OdorSpace SET


%% Visualize the results.
% Show how the clusters are different from each other (response vectors)

% Set the ylim (zscored)
ylim_range = [-5 15];

% Only layer 1
iLayer = 1;

% Prepare a figure for the response vectors
hFig_vec = figure('Units','normalized','Position',[0 0 0.5 1]);
% Prepare a figure for the allocation of vectors per phase
hFig_pi = figure('Units','normalized','Position',[0 0 0.5 1]);

% Create some variables
PCA_centroid = [];
PCA_centroid_pca = [];
PCA_label = {}; cnt = 1;

% Iterate over all clusters
for iC = 1:OdorSpace.(SET.layers{iLayer}).RespVec.best_k

    % Pool data across phases
    resp_vec = [];
    for iPhase = 1:length(SET.phases)
        resp_vec = [resp_vec; OdorSpace.(SET.layers{iLayer}).RespVec.centroid_orig{iPhase}{iC}];
    end%iPhase
    resp_vec = resp_vec(~all(isnan(resp_vec), 2),:);

    % Get average response
    avg = nanmean(resp_vec,1);
    if size(resp_vec,1)==1
        ci = nan(2,length(resp_vec));
    elseif size(resp_vec,1)==2
        ci = nan(2,length(resp_vec));
    else
        ci = bootci(5000, @nanmean, resp_vec);
    end

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
        % text(iStim, avg(iStim), CLD{iStim}, 'HorizontalAlignment', 'center', 'Color', [1 0 0])
    end%iStim
    % Cosmetics
    xlim([0.5, length(avg)+0.55])
    ylim(ylim_range)
    title(iC)
    box on
    axis square
    set(gca, 'XTick', 1:length(SET.StimList_RespVec), 'XTickLabel', SET.StimList_RespVec, 'TickLabelInterpreter', 'none')

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
        prop = length(find(OdorSpace.(SET.layers{iLayer}).RespVec.ID == iC & OdorSpace.(SET.layers{iLayer}).RespVec.info(:,1) == iPhase)) / sum(OdorSpace.(SET.layers{iLayer}).RespVec.info(:,1) == iPhase);
        prop = [1- prop; prop];
        % Plot pie chart
        nexttile; hold on
        ax = gca;
        pie(ax, prop);
        colormap(ax, col)
        axis square equal tight off
        title([SET.phases{iPhase}(1), num2str(iC)])
    end%iPhase

end%iC

% Save figures
% --- allocation
figure(hFig_pi)
export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\resp_vec_porp_pie'], '-pdf', '-painters')
close(hFig_pi)
% --- response vectors
figure(hFig_vec)
export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\kmeans_cluster_responses'], '-pdf', '-painters')
close(hFig_vec)

% Save
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')
% Clean up
clearvars -except OdorSpace SET


%% Visualize the results.
% Show how response motif assignments change between stimuli

% Set how mmany of the PCs should be kept to predict social phenoype. Note,
% this is hard-coded to increase performance.
curr_pca_expl = [83, 87];

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)

        % Get response motifs as corresponding grouping factors
        data_grouping = OdorSpace.(SET.layers{iLayer}).Motif.grouping.(SET.phases{iPhase});
        data_grouping_info = OdorSpace.(SET.layers{iLayer}).Motif.grouping_info.(SET.phases{iPhase});

        % Iterate over all stimulus and motif-change combinations to get an
        % overview for the dynamics of motif transitions.
        % --- Get list of animals
        ani_list = unique(data_grouping_info);
        % --- Get list of motifs
        motif_list = unique(data_grouping);
        % --- Get all possible combinations of motifs
        motif_change_pairs = repmat(motif_list, 1, 2);
        motif_change_list = perms(motif_list);  motif_change_list = motif_change_list(:,1:2);
        motif_change_list = [motif_change_list; motif_change_pairs];
        motif_change_list = fliplr(unique(motif_change_list, 'rows'));
        clear motif_change_pairs
        % --- Get all possible combinations of stimuli
        stim_change_list = perms(1:length(SET.StimList_Motif));  stim_change_list = stim_change_list(:,1:2);
        stim_change_list = fliplr(unique(stim_change_list, 'rows'));
        clear stim_change_pairs

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
                        data_grouping_info == ani_list(iAni));
                    % Normalize by the number of origin granules
                    origin_cnt = sum(...
                        data_grouping(:,stim_change_list(iStim,1))==motif_change_list(iMotif,1) & ...
                        data_grouping_info == ani_list(iAni));
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
        OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.(SET.phases{iPhase}).stim_change_list = stim_change_list;
        OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.(SET.phases{iPhase}).motif_change_list = motif_change_list;
        OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.(SET.phases{iPhase}).motif_dynamics = motif_dynamics;
        OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.(SET.phases{iPhase}).motif_dynamics_features = motif_dynamics_features;
        OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.(SET.phases{iPhase}).motif_change_list_labels = motif_change_list_labels;
        OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.(SET.phases{iPhase}).stim_change_list_labels = stim_change_list_labels;
        OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.(SET.phases{iPhase}).motif_dynamics_features_labels = motif_dynamics_features_labels;
    end%iPhase

    % Compare gregarious and solitarious motif dynamics
    PoolData.Stats.motif_dynamics.seed=1234;
    PoolData.Stats.motif_dynamics.n = 5e3;
    PoolData.Stats.motif_dynamics.nComp = 1;
    PoolData.Stats.motif_dynamics.pValues = nan(length(stim_change_list), length(motif_change_list));
    for iMotif = 1:size(motif_change_list,1)
        % Iterate over all stimulus changes
        for iStim = 1:size(stim_change_list,1)
            greg = squeeze(OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.gregarious.motif_dynamics(iStim, iMotif, :));
            soli = squeeze(OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.solitarious.motif_dynamics(iStim, iMotif, :));
            greg = greg(~isnan(greg)); soli = soli(~isnan(soli));
            PoolData.Stats.(SET.phases{iPhase}).motif_dynamics.pValues(iStim, iMotif) = ranksum(greg, soli);
        end%iStim
    end%iMotif

    % Show heat maps for both phases and the resulting difference
    hFig = figure('color', 'w', 'units', 'normalized', 'Position', [0 0 1 1 ]);
    % Get data
    img_greg = nanmean(OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.gregarious.motif_dynamics, 3);
    img_soli = nanmean(OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.solitarious.motif_dynamics, 3);
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
        OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.gregarious.motif_dynamics_features;...
        OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.solitarious.motif_dynamics_features];
    all_labels = [...
        zeros(size(OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.gregarious.motif_dynamics_features,1),1);...
        ones(size(OdorSpace.(SET.layers{iLayer}).Motif.Dynamics.solitarious.motif_dynamics_features,1),1)];
    % Kick out NaN
    ind = find(~sum(isnan(all_motifs),2)>0);
    all_motifs = all_motifs(ind,:);
    all_labels = all_labels(ind,:);
    label_list = unique(all_labels);

    opt = statset('pca');
    opt.MaxIter = 1e6;
    opt.TolFun = 1e-8;
    opt.TolFun = 1e-8;
    [~, all_motifs, ~, ~, expl] = pca(all_motifs, 'Options', opt);
    ind = find(cumsum(expl) >= curr_pca_expl(iLayer));
    all_motifs = all_motifs(:,1:ind);

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
        rng(iRep)
        mdl = fitcdiscr(all_motifs(idx_train,:),all_labels(idx_train,:), 'DiscrimType', 'pseudolinear', 'Weights', w);
        % Predict label for held-out data
        model_classification(iRep) =  predict(mdl, all_motifs(idx_test,:));
    end%iRep

    hFig = figure;
    plotconfusion(all_labels',model_classification')
    set(gca, 'XTickLabel', {'greg', 'soli', ' '})
    set(gca, 'YTickLabel', {'greg', 'soli', ' '})
    export_fig(['ConfocalAnalysis_OdorSpace\',SET.layers{iLayer},'\motif_transition_prediction'], '-pdf', '-painters')
    close(hFig)

end%iLayer

% Save
save('ConfocalAnalysis_OdorSpace\OdorSpace.mat', 'OdorSpace', 'SET')
% Clean up
clearvars -except OdorSpace SET


%% Visualize the results.
% Show how how the proportion of response motifs changes between stimuli

% Iterate over all layers
for iLayer = 1:length(SET.layers)
    % Iterate over both phases
    for iPhase = 1:length(SET.phases)
        % Get each ROI's response motif ID for each stimulus
        motifs = OdorSpace.(SET.layers{iLayer}).Motif.ID(OdorSpace.(SET.layers{iLayer}).Motif.info(:,1) == iPhase,:);
        stims = OdorSpace.(SET.layers{iLayer}).Motif.info(OdorSpace.(SET.layers{iLayer}).Motif.info(:,1) == iPhase,3);
        % Get info on animals
        info = OdorSpace.(SET.layers{iLayer}).Motif.grouping_info.(SET.phases{iPhase});
        % Get list of animals
        ani_list = unique(info);
        % Get list of motifs
        id_list = unique(OdorSpace.(SET.layers{iLayer}).Motif.grouping.(SET.phases{iPhase})(:));
        % Get the average proportion of motif IDs per stimulus across all
        % animals
        prop_data{iPhase} = zeros(OdorSpace.(SET.layers{iLayer}).Motif.best_k, length(SET.StimList_Motif));
        % Iterate over all stimui
        for iStim = 1:length(SET.StimList_Motif)
            % Iterate over all motifs
            for iC = 1:length(id_list)
                % Iterate over all animals and get the proportion of
                % motifs for the current stimulus
                for iAni = 1:length(ani_list)
                    ind = find(info == ani_list(iAni));
                    p = sum(OdorSpace.(SET.layers{iLayer}).Motif.grouping.(SET.phases{iPhase})(ind,iStim) == id_list(iC));
                    p = p / length(ind);
                    prop_data{iPhase}(iC, iStim) = prop_data{iPhase}(iC, iStim) + p;
                end%iAni
                prop_data{iPhase}(iC, iStim) = prop_data{iPhase}(iC, iStim) / length(ani_list);
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
        colormap(cols(2:end,:)); colorbar; clim([1 OdorSpace.(SET.layers{iLayer}).Motif.best_k])
        % Add labels
        set(gca, 'ThetaTickLabel', SET.StimList_Motif, 'TickLabelInterpreter', 'none')
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

% Clean up
clearvars -except OdorSpace SET