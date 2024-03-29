%% Widefield_PureOdorAnalysis_Population
% This script pools the data across all animals for both social phenotypes.
% It plots the respective activity time courses with corresponding 
% swarm/violin plots.
%
% Version:
% 05-Jan-2024 (R2023a)

% Prepare
clc; clear all; close all
warning('off')

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig)
addpath(genpath('...\export_fig'))
mkdir('PureOdorAnalysis_Population')

%% Settings

% Set paths
SET.main_path = '...';

% Seth the two phases
SET.phases = {'gregarious', 'solitarious'};

% Set framerate
SET.fps = 10;

% Give list of stimuli (We have this in the meta data, but good to double
% check)
SET.StimList = {...
    'N2';...
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

% Ylim for TC and rasterplots
SET.TC_lim = [-0.25 3];
SET.TC_lim_swarm = [0 5];

%% Get data
% Iterate over all animals and pool data from all pockets. Later use this
% to greate three main plots: a raster plot showing all TCs as an image,
% the resulting average TCs, and violin/swarm plots for the avg. activity
% in the active window


% Iterate over both phases
for iPhase = 1:length(SET.phases)

    % Get overview of animal folders
    path2animal = [SET.main_path, SET.phases{iPhase},'\'];
    curr.dir.all = dir(path2animal);    

    % Prepare variables
    poolTC = [];
    poolAvg = [];
    poolInfo = [];

    % Iterate over all animals
    for iAni = 1:size(curr.dir.all,1)
        % Check whether this is what we are looking for, i.e. whether the
        % string "Animal" is present
        if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir

            % Get the path to the current animal
            curr.dir.animal = [path2animal,curr.dir.all(iAni).name,'\'];

            % Get meta data and segmentation
            load([curr.dir.animal,'03_Data_Processed\meta_info.mat'])
            load([curr.dir.animal,'03_Data_Processed\img_segmentation.mat'])

            % Get list of valid pockets. Take all pockets that were active
            % at least once.
            mask_list = zeros(length(meta.unique_stim),1);
            for iStim = 1:length(SET.StimList)
                mask_list = mask_list+ (meta.unique_stim == SET.StimList{iStim});
            end
            mask = sum(Segmentation.pocket_active_img(:,:,find(mask_list)),3);
            mask = mask>0;
            pocket_list_valid = unique(Segmentation.pockets_labeled.*double(mask));
            pocket_list_valid(pocket_list_valid==0) = [];
            clear mask mask_list iStim

            % Iterate over all valid pockets of all stimuli
            curr_poolTC = nan(length(pocket_list_valid)*length(SET.StimList), meta.duration);
            curr_poolAvg = nan(length(pocket_list_valid)*length(SET.StimList), 1); 
            curr_poolInfo = nan(length(pocket_list_valid)*length(SET.StimList), 2);            
            cnt=1;
            for iStim = 1:length(SET.StimList)
                % Load data
                currData = load([curr.dir.animal,'03_Data_Processed\',SET.StimList{iStim},'.mat']);
                currData = reshape(currData.ImageStream.(SET.StimList{iStim}), [size(currData.ImageStream.(SET.StimList{iStim}),1)*size(currData.ImageStream.(SET.StimList{iStim}),2), size(currData.ImageStream.(SET.StimList{iStim}),3)]);
                % Iterate over all pockets
                for iP = 1:length(pocket_list_valid)
                    % Get index positions of the current pocker
                    idx = find(Segmentation.pockets_labeled == pocket_list_valid(iP));
                    % Get timecourse, activity, and info
                    curr_poolTC(cnt,:) = nanmean(currData(idx,:));
                    curr_poolAvg(cnt,1) = nanmean(curr_poolTC(cnt,meta.activeRegion));
                    curr_poolInfo(cnt,:) = [iAni, iStim];
                    % Counter
                    cnt=cnt+1;
                end%ip
            end%iStim        

        % Pool different animals
        poolTC = [poolTC; curr_poolTC];
        poolAvg = [poolAvg; curr_poolAvg];
        poolInfo = [poolInfo; curr_poolInfo];

        end%if animal
    end%iAni


    % ---------------------------------------------------------------------
    % Sort pool by activity. For this, get the indices in order and apply
    % this to the data.
    sort_idx = [];
    for iStim = 1:length(SET.StimList)
        idx = find(poolInfo(:,2)==iStim);
        [~,I] = sort(poolAvg(idx), 'descend');
        sort_idx = [sort_idx; idx(I)];
    end
    poolTC = poolTC(sort_idx,:);
    poolAvg = poolAvg(sort_idx,:);
    poolInfo = poolInfo(sort_idx,:);

    % Get grand means for TCs and AVGs
    animal_list = unique(poolInfo(:,1));
    for iStim = 1:length(SET.StimList)
        grandMean_TC.(SET.StimList{iStim}) = nan(length(animal_list), meta.duration);
        grandMean_AVG.(SET.phases{iPhase}).(SET.StimList{iStim}) = nan(length(animal_list), 1);
        for iAni = 1:length(animal_list)
            idx = find(poolInfo(:,1)==animal_list(iAni) & poolInfo(:,2)==iStim);
            grandMean_TC.(SET.StimList{iStim})(iAni,:) = nanmean(poolTC(idx,:));
            grandMean_AVG.(SET.phases{iPhase}).(SET.StimList{iStim})(iAni,1) = nanmean(grandMean_TC.(SET.StimList{iStim})(iAni,meta.activeRegion));
        end%iAni
    end%iStim


    % ---------------------------------------------------------------------
    % Plot time courses
    hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 1 0.5]);
    % Get data for air control
    avg_air = mean(bootstrp(5000, @mean, grandMean_TC.(SET.StimList{1})));
    ci_air = bootci(5000, @mean, grandMean_TC.(SET.StimList{1}));
    time_vec = ([1:size(poolTC,2)] - meta.baseline(1)) / SET.fps;
    % Iterate over all stimuli
    for iStim = 2:length(SET.StimList)
        subplot(1,length(SET.StimList)-1,iStim-1); hold on
        rectangle(...
            'Position',[time_vec(meta.activeRegion(1)), -10, time_vec(meta.activeRegion(end))-time_vec(meta.activeRegion(1)), 20],...
            'FaceColor',[0.5 0.5 0.5],...
            'EdgeColor','none')
        % Plot air control
        plot(time_vec, avg_air, 'color', SET.Colors(1,:))
        plot(time_vec, ci_air, 'color', SET.Colors(1,:))
        % Plot response to stimulus
        avg = mean(bootstrp(5000, @mean, grandMean_TC.(SET.StimList{iStim})));
        ci = bootci(5000, @mean, grandMean_TC.(SET.StimList{iStim}));
        plot(time_vec, avg, 'color', SET.Colors(iStim,:))
        plot(time_vec, ci, 'color', SET.Colors(iStim,:))
    end%iStim
    % Iterate over all stimuli again to apply cosmetics
    for iStim = 2:length(SET.StimList)
        subplot(1,length(SET.StimList)-1,iStim-1); hold on
        ylim(SET.TC_lim)
        xlim([0, time_vec(meta.stimulus_onset(2))])
        ylabel('\DeltaF/F (%)')
        xlabel('time (s)')
        title(SET.StimList{iStim}, 'Interpreter','none')
    end%iStim
    export_fig(['PureOdorAnalysis_Population', '\', SET.phases{iPhase}, '_TCs'], '-pdf')
    close(hFig)


    % ---------------------------------------------------------------------
    % Plot raster plot
    hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 0.5 1]);
    imagesc(poolTC)
    hold on
    for iStim = 1:length(SET.StimList)-1
        plot([1, meta.duration], 0.5 + [iStim iStim]*size(poolTC,1)/length(SET.StimList), 'w')
    end
    axis off
    xlim([meta.baseline(1), meta.stimulus_onset(2)-1])
    caxis(SET.TC_lim)
    colormap(SET.HeatmapColors.(SET.phases{iPhase}))
    colorbar
    export_fig(['PureOdorAnalysis_Population', '\', SET.phases{iPhase}, '_TRaster'], '-pdf')
    close(hFig)


    % ---------------------------------------------------------------------
    % Plot violin/swarm plots
    hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 0.5 1]);
    % Iterate over all stimuli
    for iStim = 2:length(SET.StimList)
        subplot(1,length(SET.StimList)-1,iStim-1); hold on

        plot([-0.5 0.75],[1 1]*mean(bootstrp(5000, @mean, grandMean_AVG.(SET.phases{iPhase}).(SET.StimList{1}))), 'Color', SET.Colors(1,:))

        % Violin
        vioData = grandMean_AVG.(SET.phases{iPhase}).(SET.StimList{iStim});
        properties.MinVal =             min(vioData);                 
        properties.MaxVal =             max(vioData);                 
        properties.AvgType =            'mean';             
        properties.MeanCol =            'k';                
        properties.MeanSymbol =         'line';             
        properties.MeanWidth =          2;                  
        properties.SeparateOutliers =   0;                  
        Widefield_SubFcn.violinplot_advanced(vioData, 0, 0.5, properties)
        plot([0,0], [min(vioData), max(vioData)], 'k')
        clear vioData properties

        % Swarm
        beeData = grandMean_AVG.(SET.phases{iPhase}).(SET.StimList{iStim});
        properties.MarkerType =         'o';                 
        properties.MarkerFaceColor =    SET.Colors(iStim,:); 
        properties.MarkerSize =         5;
        Widefield_SubFcn.beeswarmplot_advanced(beeData, 0.25, 0.5, properties)
        % Add mean and CI
        avg = mean(bootstrp(5000, @mean, beeData));
        ci = bootci(5000, @mean, beeData);
        plot([0.125 0.375], [avg avg], 'k')
        plot([0.25 0.25], ci, 'k')
        clear beeData properties

        % Cosmetics
        xlim([-0.5 0.75])
        ylim(SET.TC_lim_swarm)
        ylabel('\DeltaF/F (%)')
        xlabel(SET.StimList{iStim}, 'Interpreter','none')

    end%iStim
    export_fig(['PureOdorAnalysis_Population', '\', SET.phases{iPhase}], '-pdf')
    close(hFig)

    % ---------------------------------------------------------------------
    % Lines connecting swarm plots
    hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 0.5 1]);
    % Iterate over all stimuli
    for iStim = 2:length(SET.StimList)
        for jStim = 2:length(SET.StimList)
            if iStim ~= jStim
                nexttile; hold on
                % Get data
                beeData1 = grandMean_AVG.(SET.phases{iPhase}).(SET.StimList{iStim});
                beeData2 = grandMean_AVG.(SET.phases{iPhase}).(SET.StimList{jStim});
                % Connect
                plot([beeData1(:), beeData2(:)]', 'Color', [0.251 0.251 0.251])
                % Cosmetics
                xlim([0 3])
                ylim(SET.TC_lim_swarm)
                ylabel('\DeltaF/F (%)')
                xlabel(SET.StimList{iStim}, 'Interpreter','none')
                set(gca, 'XTick', [1 2], 'XTickLabel', {SET.StimList{[iStim, jStim]}}, 'TickLabelInterpreter', 'none')
                ylim([0 5])
            end%iif
        end%jStim
    end%iStim
    export_fig(['PureOdorAnalysis_Population', '\', SET.phases{iPhase}, '_ConnectingLines'], '-pdf')
    close(hFig)

    % Statistics
    Stats.seed=1234;
    Stats.n = 5e7;
    Stats.nComp = 7;
    % --- food vs mix
    [...
        Stats.(SET.phases{iPhase}).Laa_vs_LaaLct.p,...
        Stats.(SET.phases{iPhase}).Laa_vs_LaaLct.s,...
        ~,...
        Stats.(SET.phases{iPhase}).Laa_vs_LaaLct.c] = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', grandMean_AVG.(SET.phases{iPhase}).ZHAE2, grandMean_AVG.(SET.phases{iPhase}).ZHAE2_COL, Stats.n, Stats.seed, Stats.nComp);

    % --- social vs mix
    [...
        Stats.(SET.phases{iPhase}).Lct_vs_LaaLct.p,...
        Stats.(SET.phases{iPhase}).Lct_vs_LaaLct.s,...
        ~,...
        Stats.(SET.phases{iPhase}).Lct_vs_LaaLct.c] = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', grandMean_AVG.(SET.phases{iPhase}).COL, grandMean_AVG.(SET.phases{iPhase}).ZHAE2_COL, Stats.n, Stats.seed, Stats.nComp);

    % --- food vs social
    [...
        Stats.(SET.phases{iPhase}).Laa_vs_Lct.p,...
        Stats.(SET.phases{iPhase}).Laa_vs_Lct.s,...
        ~,...
        Stats.(SET.phases{iPhase}).Laa_vs_Lct.c] = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', grandMean_AVG.(SET.phases{iPhase}).ZHAE2, grandMean_AVG.(SET.phases{iPhase}).COL, Stats.n, Stats.seed, Stats.nComp);

    % --- food vs. hex
    [...
        Stats.(SET.phases{iPhase}).Laa_vs_Hex.p,...
        Stats.(SET.phases{iPhase}).Laa_vs_Hex.s,...
        ~,...
        Stats.(SET.phases{iPhase}).Laa_vs_Hex.c] = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', grandMean_AVG.(SET.phases{iPhase}).HEX3, grandMean_AVG.(SET.phases{iPhase}).ZHAE2, Stats.n, Stats.seed, Stats.nComp);

    % --- food vs. food_hex
    [...
        Stats.(SET.phases{iPhase}).Laa_vs_LaaHex.p,...
        Stats.(SET.phases{iPhase}).Laa_vs_LaaHex.s,...
        ~,...
        Stats.(SET.phases{iPhase}).Laa_vs_LaaHex.c] = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', grandMean_AVG.(SET.phases{iPhase}).ZHAE2_HEX3, grandMean_AVG.(SET.phases{iPhase}).ZHAE2, Stats.n, Stats.seed, Stats.nComp);

    % --- food vs. ocl
    [...
        Stats.(SET.phases{iPhase}).Laa_vs_Ocl.p,...
        Stats.(SET.phases{iPhase}).Laa_vs_Ocl.s,...
        ~,...
        Stats.(SET.phases{iPhase}).Laa_vs_LaaHex.c] = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', grandMean_AVG.(SET.phases{iPhase}).OC3L2, grandMean_AVG.(SET.phases{iPhase}).ZHAE2, Stats.n, Stats.seed, Stats.nComp);

    % --- food vs. food_ocl
    [...
        Stats.(SET.phases{iPhase}).Laa_vs_LaaOcl.p,...
        Stats.(SET.phases{iPhase}).Laa_vs_LaaOcl.s,...
        ~,...
        Stats.(SET.phases{iPhase}).Laa_vs_LaaOcl.c] = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', grandMean_AVG.(SET.phases{iPhase}).ZHAE2_OC3L2, grandMean_AVG.(SET.phases{iPhase}).ZHAE2, Stats.n, Stats.seed, Stats.nComp);
    

    clearvars -except iPhase SET Stats grandMean_AVG
end%iPhase

% Statistics for gregarious vs. solitarious
Stats.greg_vs_soli.seed=1234;
Stats.greg_vs_soli.n = 5e7;
Stats.greg_vs_soli.nComp = 1;
% Iterate over all stimuli
for iStim = 2:length(SET.StimList)
[...
        Stats.greg_vs_soli.(SET.StimList{iStim}).p,...
        Stats.greg_vs_soli.(SET.StimList{iStim}).s,...
        ~,...
        Stats.greg_vs_soli.(SET.StimList{iStim}).c] = Widefield_SubFcn.BootstrapHypothesisTesting('two-sample', grandMean_AVG.(SET.phases{1}).(SET.StimList{iStim}), grandMean_AVG.(SET.phases{2}).(SET.StimList{iStim}), Stats.greg_vs_soli.n, Stats.greg_vs_soli.seed, Stats.greg_vs_soli.nComp);
end%iStim

% Save
save(['PureOdorAnalysis_Population', '\', 'statistics.mat'], 'Stats')