%% Confocal_ConfocalAnalysis_Population
% ...
%
% Version:
% 2-April-2023 (R2023a) Yannick Günzel

% Prepare
clc; clear all; close all
warning('off')

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig)
addpath(genpath('...\GitHub\export_fig'))
mkdir('ConfocalAnalysis_Population')

%% Settings

% Set paths
SET.main_path = '...\Data\Physiology\Cal520_Confocal\';

% Give list of phases
SET.phases = {'gregarious'; 'solitarious'};

% Give layer names
SET.layers = {'layer01_top', 'layer02_middle'};

% Give list of stimuli (We have this in the meta data, but good to double
% check)
SET.StimList = {...
    'ZHAE';...
    'COL';...
    'ZHAECOL'};

% Set framerate
SET.fps = 1;

% Cosmetics
SET.Colors = [...
    102,166,030;...ZHAE
    217,095,002;...Col
    231,041,138;...ZZHAECOL
    ]/255;
SET.PhaseColors.gregarious = [191 0 0]/255;
SET.PhaseColors.solitarious = [0 0 191]/255;

% Colormap for heatmaps
SET.HeatmapColors.gregarious = Confocal_SubFcn.ColMapPlasma(1000);
SET.HeatmapColors.solitarious = [...
    SET.HeatmapColors.gregarious(:,2),...
    SET.HeatmapColors.gregarious(:,1),...
    SET.HeatmapColors.gregarious(:,3)];
SET.HeatmapColors.other = gray(1000);


% Ylim for TC and rasterplots
SET.TC_lim = [-1 6;-1 9];
SET.TC_lim_swarm = [-1 8; -1 15];
SET.Bliss_lim = [-9 9; -9 9];
SET.Bliss_pdf_lim = [0, 0.35; 0, 0.35];

%% Get data
% Iterate over all animals and pool data from all pockets. Later use this
% to greate three main plots: a raster plot showing all TCs as an image,
% the resulting average TCs, and violin/swarm plots for the avg. activity
% in the active window


% Iterate over both phases
for iPhase = 1:length(SET.phases)
    for iLayer = 1:length(SET.layers)

        % Get overview of animal folders
        path2animal = [SET.main_path, SET.phases{iPhase},'\',SET.layers{iLayer},'\'];
        curr.dir.all = dir(path2animal);

        % Prepare variables
        poolTC = [];
        poolAvg = [];
        poolBinaryBliss = [];
        poolAvgPosBliss = [];
        poolAvgNegBliss = [];
        poolPdfBliss = [];
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

                % Get time vector to cut data
                time_vec = meta.baseline(1):(meta.activeRegion(end)+length(meta.baseline)-1);
                % Get indices of active region
                active_region = [find(time_vec==meta.activeRegion(1)), find(time_vec==meta.activeRegion(end))];

                % Iterate over all valid pockets of all stimuli
                curr_poolTC = nan(length(pocket_list_valid)*length(SET.StimList), length(time_vec));
                curr_poolAvg = nan(length(pocket_list_valid)*length(SET.StimList), 1);
                curr_poolInfo = nan(length(pocket_list_valid)*length(SET.StimList), 2);
                cnt=1;
                for iStim = 1:length(SET.StimList)
                    % Load data
                    currData = load([curr.dir.animal,'03_Data_Processed\',SET.StimList{iStim},'.mat']);
                    currData = reshape(currData.ImageStream.(SET.StimList{iStim}), [size(currData.ImageStream.(SET.StimList{iStim}),1)*size(currData.ImageStream.(SET.StimList{iStim}),2), size(currData.ImageStream.(SET.StimList{iStim}),3)]);
                    % Cute data
                    currData = currData(:,time_vec);
                    % Iterate over all pockets
                    for iP = 1:length(pocket_list_valid)
                        % Get index positions of the current pocker
                        idx = find(Segmentation.pockets_labeled == pocket_list_valid(iP));
                        % Get timecourse, activity, and info
                        curr_poolTC(cnt,:) = nanmean(currData(idx,:));
                        curr_poolAvg(cnt,1) = nanmean(curr_poolTC(cnt,active_region(1):active_region(2)));
                        curr_poolInfo(cnt,:) = [iAni, iStim];
                        % Counter
                        cnt=cnt+1;
                    end%ip
                end%iStim

                % Get the Bliss score to predict the combined effect of
                % Laa and Lct and if they are acting independently of
                % each other.
                % --- Get the correct indices
                idx_1 = find(curr_poolInfo(:,2)==find(strcmp(SET.StimList, 'ZHAE')));
                idx_2 = find(curr_poolInfo(:,2)==find(strcmp(SET.StimList, 'COL')));
                idx_3 = find(curr_poolInfo(:,2)==find(strcmp(SET.StimList, 'ZHAECOL')));
                m1 = curr_poolTC(idx_1,:)/100;
                m2 = curr_poolTC(idx_2,:)/100;
                m3 = curr_poolTC(idx_3,:)/100;
                % --- Calculate the Bliss expected effect
                f_expected = (m1/2 + m2/2) - ((m1/2).*(m2/2));
                % --- Get the observed values
                f_observed = m3;
                % --- Compare the Bliss expected effect with the observed
                bliss_score =  100*(f_observed - f_expected);

                % --- Subtract baseline
                bliss_score = bliss_score - mean(bliss_score(:,1:active_region(1)-1),2);

                % Get Bliss scores in active region
                active_bliss = mean(bliss_score(:,active_region(1):active_region(2)),2);
                % Get the PDF for all granules
                bliss_density = ksdensity(active_bliss(:), linspace(-100,100,1000), 'BoundaryCorrection','reflection');
                bliss_density = 100*(bliss_density/length(active_bliss));

                % Pool different animals
                poolTC = [poolTC; curr_poolTC];
                poolAvg = [poolAvg; curr_poolAvg];
                poolBinaryBliss = [poolBinaryBliss; mean(active_bliss>0)];
                poolAvgPosBliss = [poolAvgPosBliss; mean(active_bliss(active_bliss>0))];
                poolAvgNegBliss = [poolAvgNegBliss; mean(active_bliss(active_bliss<0))];
                poolPdfBliss = [poolPdfBliss; bliss_density];
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
        TC_lim = [];
        for iStim = 1:length(SET.StimList)
            grandMean_TC.(SET.StimList{iStim}) = nan(length(animal_list), length(time_vec));
            grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.StimList{iStim}) = nan(length(animal_list), 1);
            for iAni = 1:length(animal_list)
                idx = find(poolInfo(:,1)==animal_list(iAni) & poolInfo(:,2)==iStim);
                grandMean_TC.(SET.StimList{iStim})(iAni,:) = nanmean(poolTC(idx,:));
                grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.StimList{iStim})(iAni,1) = nanmean(grandMean_TC.(SET.StimList{iStim})(iAni,active_region(1):active_region(2)));
            end%iAni
            TC_lim = [TC_lim; min(grandMean_TC.(SET.StimList{iStim})(:)), max(grandMean_TC.(SET.StimList{iStim})(:))];
        end%iStim


        % ---------------------------------------------------------------------
        % Plot time courses
        hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 1 0.5]);
        hold on
        xvec = time_vec-time_vec(1) / SET.fps;
        rectangle('position', [active_region(1), floor(min(TC_lim(:,1))*2)/2, active_region(2)-active_region(1), ceil(max(TC_lim(:,2))*2)/2-floor(min(TC_lim(:,1))*2)/2],...
            'FaceColor',[0.5 0.5 0.5],...
            'EdgeColor','none')
        % Iterate over all stimuli
        for iStim = 1:length(SET.StimList)
            % Plot response to stimulus
            avg = mean(bootstrp(5000, @mean, grandMean_TC.(SET.StimList{iStim})));
            ci = bootci(5000, @mean, grandMean_TC.(SET.StimList{iStim}));
            plot(xvec, avg, 'color', SET.Colors(iStim,:))
            plot(xvec, ci, 'color', SET.Colors(iStim,:))
        end%iStim
        % Cosmetics
        ylim(SET.TC_lim(iLayer,:))
        xlim([xvec(1), xvec(end)])
        ylabel('\DeltaF/F (%)')
        xlabel('time (s)')
        title([SET.phases{iPhase},' | ', SET.layers{iLayer}], 'Interpreter','none')
        % Save
        export_fig(['ConfocalAnalysis_Population', '\', SET.phases{iPhase}, '_', SET.layers{iLayer}, '_TCs'], '-pdf')
        close(hFig)


        % ---------------------------------------------------------------------
        % Plot raster plot
        hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 0.5 1]);
        imagesc(poolTC)
        hold on
        for iStim = 1:length(SET.StimList)-1
            plot([1, size(poolTC,2)], 0.5 + [iStim iStim]*size(poolTC,1)/length(SET.StimList), 'w')
        end
        axis off
        xlim([1, size(poolTC,2)])
        caxis(SET.TC_lim(iLayer,:))
        colormap(SET.HeatmapColors.(SET.phases{iPhase}))
        colorbar
        export_fig(['ConfocalAnalysis_Population', '\', SET.phases{iPhase}, '_', SET.layers{iLayer}, '_TRaster'], '-pdf')
        close(hFig)


        % ---------------------------------------------------------------------
        % Plot violin/swarm plots
        hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 0.5 1]);
        % Iterate over all stimuli
        for iStim = 1:length(SET.StimList)
            subplot(1,length(SET.StimList),iStim); hold on

            % Violin
            vioData = grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.StimList{iStim});
            properties.MinVal =             min(vioData);
            properties.MaxVal =             max(vioData);
            properties.AvgType =            'mean';
            properties.MeanCol =            'k';
            properties.MeanSymbol =         'line';
            properties.MeanWidth =          2;
            properties.SeparateOutliers =   0;
            Confocal_SubFcn.violinplot_advanced(vioData, 0, 0.5, properties)
            plot([0,0], [min(vioData), max(vioData)], 'k')
            clear vioData properties

            % Swarm
            beeData = grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.StimList{iStim});
            properties.MarkerType =         'o';
            properties.MarkerFaceColor =    SET.Colors(iStim,:);
            properties.MarkerSize =         5;
            Confocal_SubFcn.beeswarmplot_advanced(beeData, 0.25, 0.5, properties)
            % Add mean and CI
            avg = mean(bootstrp(5000, @mean, beeData));
            ci = bootci(5000, @mean, beeData);
            plot([0.125 0.375], [avg avg], 'k')
            plot([0.25 0.25], ci, 'k')
            clear beeData properties

            % Cosmetics
            xlim([-0.5 0.75])
            ylim(SET.TC_lim_swarm(iLayer,:))
            ylabel('\DeltaF/F (%)')
            xlabel(SET.StimList{iStim}, 'Interpreter','none')

        end%iStim
        export_fig(['ConfocalAnalysis_Population', '\', SET.phases{iPhase},'_',SET.layers{iLayer},'_violinswarm'], '-pdf')
        close(hFig)

        % ---------------------------------------------------------------------
        % Lines connecting swarm plots
        hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 0.5 1]);
        % Iterate over all stimuli
        for iStim = 1:length(SET.StimList)
            for jStim = 1:length(SET.StimList)
                if iStim ~= jStim
                    nexttile; hold on
                    % Get data
                    beeData1 = grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.StimList{iStim});
                    beeData2 = grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.StimList{jStim});
                    % Connect
                    plot([beeData1(:), beeData2(:)]', 'Color', [0.251 0.251 0.251])
                    % Cosmetics
                    xlim([0 3])
                    ylim(SET.TC_lim_swarm(iLayer,:))
                    ylabel('\DeltaF/F (%)')
                    xlabel(SET.StimList{iStim}, 'Interpreter','none')
                    set(gca, 'XTick', [1 2], 'XTickLabel', {SET.StimList{[iStim, jStim]}}, 'TickLabelInterpreter', 'none')
                end%iif
            end%jStim
        end%iStim
        export_fig(['ConfocalAnalysis_Population', '\', SET.phases{iPhase},'_',SET.layers{iLayer},'_ConnectingLines'], '-pdf')
        close(hFig)


        % ---------------------------------------------------------------------
        % Plot avg+-CI
        hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 0.5 1]);
        % Iterate over all stimuli
        for iStim = 1:length(SET.StimList)
            subplot(1,length(SET.StimList),iStim); hold on
            crossData = grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).(SET.StimList{iStim});
            avg = mean(bootstrp(5000, @mean, crossData));
            ci = bootci(5000, @mean, crossData);
            plot([-0.125 0.125], [avg avg], 'color', SET.Colors(iStim,:))
            plot([0 0], ci, 'color', SET.Colors(iStim,:))
            clear crossData properties
            % Cosmetics
            xlim([-0.5 0.5])
            ylim(SET.TC_lim(iLayer,:))
            ylabel('\DeltaF/F (%)')
            xlabel(SET.StimList{iStim}, 'Interpreter','none')
            set(gca, 'XTick', 0, 'XTickLabel', SET.StimList{iStim})
        end%iStim
        export_fig(['ConfocalAnalysis_Population', '\', SET.phases{iPhase},'_',SET.layers{iLayer},'_avgCross'], '-pdf')
        close(hFig)


        % Imaging Statistics
        Stats.(SET.phases{iPhase}).seed=1234;
        Stats.(SET.phases{iPhase}).n = 5e4;
        Stats.(SET.phases{iPhase}).nComp = 3;
        % --- food vs mix
        [...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).Laa_vs_LaaLct.p,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).Laa_vs_LaaLct.s,...
            ~,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).Laa_vs_LaaLct.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).ZHAE, grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).ZHAECOL, Stats.(SET.phases{iPhase}).n, Stats.(SET.phases{iPhase}).seed, Stats.(SET.phases{iPhase}).nComp);

        % --- social vs mix
        [...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).Lct_vs_LaaLct.p,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).Lct_vs_LaaLct.s,...
            ~,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).Lct_vs_LaaLct.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).COL, grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).ZHAECOL, Stats.(SET.phases{iPhase}).n, Stats.(SET.phases{iPhase}).seed, Stats.(SET.phases{iPhase}).nComp);

        % --- food vs social
        [...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).Laa_vs_Lct.p,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).Laa_vs_Lct.s,...
            ~,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).Laa_vs_Lct.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).ZHAE, grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).COL, Stats.(SET.phases{iPhase}).n, Stats.(SET.phases{iPhase}).seed, Stats.(SET.phases{iPhase}).nComp);



        % ---------------------------------------------------------------------
        % Plot Bliss score
        hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 1 0.5]);


        subplot(1,5,1); hold on
        % Violin
        properties.MinVal =             min(poolBinaryBliss);
        properties.MaxVal =             max(poolBinaryBliss);
        properties.AvgType =            'mean';
        properties.MeanCol =            'k';
        properties.MeanSymbol =         'line';
        properties.MeanWidth =          2;
        properties.SeparateOutliers =   0;
        Confocal_SubFcn.violinplot_advanced(poolBinaryBliss, 0, 0.5, properties)
        plot([0,0], [min(poolBinaryBliss), max(poolBinaryBliss)], 'k')
        clear properties
        % Swarm
        properties.MarkerType =         'o';
        properties.MarkerFaceColor =    SET.PhaseColors.(SET.phases{iPhase});
        properties.MarkerSize =         5;
        Confocal_SubFcn.beeswarmplot_advanced(poolBinaryBliss, 0.25, 0.5, properties)
        % Add mean and CI
        avg = mean(bootstrp(5000, @mean, poolBinaryBliss));
        ci = bootci(5000, @mean, poolBinaryBliss);
        plot([0.125 0.375], [avg avg], 'k')
        plot([0.25 0.25], ci, 'k')
        clear properties
        % --- Cosmetics
        xlim([-0.5 0.75])
        ylim([0 1])
        ylabel('prop. syn. gran.')
        % Statistics
        Stats.(SET.phases{iPhase}).BinaryBliss.seed=1234;
        Stats.(SET.phases{iPhase}).BinaryBliss.n = 5e4;
        Stats.(SET.phases{iPhase}).BinaryBliss.nComp = 1;
        [...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).BinaryBliss.p,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).BinaryBliss.s] = Confocal_SubFcn.BootstrapHypothesisTesting('one-sample', poolBinaryBliss, 0.5, Stats.(SET.phases{iPhase}).BinaryBliss.n, Stats.(SET.phases{iPhase}).BinaryBliss.seed, Stats.(SET.phases{iPhase}).BinaryBliss.nComp);


        subplot(1,5,2); hold on
        % Violin
        properties.MinVal =             min(poolAvgPosBliss);
        properties.MaxVal =             max(poolAvgPosBliss);
        properties.AvgType =            'mean';
        properties.MeanCol =            'k';
        properties.MeanSymbol =         'line';
        properties.MeanWidth =          2;
        properties.SeparateOutliers =   0;
        Confocal_SubFcn.violinplot_advanced(poolAvgPosBliss, 0, 0.5, properties)
        plot([0,0], [min(poolAvgPosBliss), max(poolAvgPosBliss)], 'k')
        clear properties
        % Swarm
        properties.MarkerType =         'o';
        properties.MarkerFaceColor =    SET.PhaseColors.(SET.phases{iPhase});
        properties.MarkerSize =         5;
        Confocal_SubFcn.beeswarmplot_advanced(poolAvgPosBliss, 0.25, 0.5, properties)
        % Add mean and CI
        avg = mean(bootstrp(5000, @mean, poolAvgPosBliss));
        ci = bootci(5000, @mean, poolAvgPosBliss);
        plot([0.125 0.375], [avg avg], 'k')
        plot([0.25 0.25], ci, 'k')
        clear properties

        subplot(1,5,2); hold on
        % Violin
        properties.MinVal =             min(poolAvgNegBliss);
        properties.MaxVal =             max(poolAvgNegBliss);
        properties.AvgType =            'mean';
        properties.MeanCol =            'k';
        properties.MeanSymbol =         'line';
        properties.MeanWidth =          2;
        properties.SeparateOutliers =   0;
        Confocal_SubFcn.violinplot_advanced(poolAvgNegBliss, 0, 0.5, properties)
        plot([0,0], [min(poolAvgNegBliss), max(poolAvgNegBliss)], 'k')
        clear properties
        % Swarm
        properties.MarkerType =         'o';
        properties.MarkerFaceColor =    SET.PhaseColors.(SET.phases{iPhase});
        properties.MarkerSize =         5;
        Confocal_SubFcn.beeswarmplot_advanced(poolAvgNegBliss, 0.25, 0.5, properties)
        % Add mean and CI
        avg = mean(bootstrp(5000, @mean, poolAvgNegBliss));
        ci = bootci(5000, @mean, poolAvgNegBliss);
        plot([0.125 0.375], [avg avg], 'k')
        plot([0.25 0.25], ci, 'k')
        clear properties
        % --- Cosmetics
        xlim([-0.5 0.75])
        ylim(SET.Bliss_lim(iLayer,:))
        ylabel('avg. Bliss of syn gran')
        % Statistics
        Stats.(SET.phases{iPhase}).AvgPosBliss.seed=1234;
        Stats.(SET.phases{iPhase}).AvgPosBliss.n = 5e4;
        Stats.(SET.phases{iPhase}).AvgPosBliss.nComp = 1;
        [...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgPosBliss.p,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgPosBliss.s] = Confocal_SubFcn.BootstrapHypothesisTesting('one-sample', poolAvgPosBliss, 0, Stats.(SET.phases{iPhase}).AvgPosBliss.n, Stats.(SET.phases{iPhase}).AvgPosBliss.seed, Stats.(SET.phases{iPhase}).AvgPosBliss.nComp);
        
        Stats.(SET.phases{iPhase}).AvgNegBliss.seed=1234;
        Stats.(SET.phases{iPhase}).AvgNegBliss.n = 5e4;
        Stats.(SET.phases{iPhase}).AvgNegBliss.nComp = 1;
        [...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgNegBliss.p,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgNegBliss.s] = Confocal_SubFcn.BootstrapHypothesisTesting('one-sample', poolAvgNegBliss, 0, Stats.(SET.phases{iPhase}).AvgNegBliss.n, Stats.(SET.phases{iPhase}).AvgNegBliss.seed, Stats.(SET.phases{iPhase}).AvgNegBliss.nComp);
        
        Stats.(SET.phases{iPhase}).AvgBliss.seed=1234;
        Stats.(SET.phases{iPhase}).AvgBliss.n = 5e4;
        Stats.(SET.phases{iPhase}).AvgBliss.nComp = 1;
        [...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgBliss.p,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgBliss.s,...
            ~,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgBliss.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', abs(poolAvgPosBliss), abs(poolAvgNegBliss), Stats.(SET.phases{iPhase}).AvgBliss.n, Stats.(SET.phases{iPhase}).AvgBliss.seed, Stats.(SET.phases{iPhase}).AvgBliss.nComp);



        % Average of all animals during active regions
        subplot(1,5,3); hold on
        % --- Get data
        % StandBliss = (sum(poolPdfBliss(:,501:end),2) - sum(poolPdfBliss(:,1:500),2)) ./ (sum(poolPdfBliss(:,501:end),2) + sum(poolPdfBliss(:,1:500),2));
        StandBliss = (poolAvgPosBliss + poolAvgNegBliss) ./ (poolAvgPosBliss + abs(poolAvgNegBliss));
        % Violin
        properties.MinVal =             min(StandBliss);
        properties.MaxVal =             max(StandBliss);
        properties.AvgType =            'mean';
        properties.MeanCol =            'k';
        properties.MeanSymbol =         'line';
        properties.MeanWidth =          2;
        properties.SeparateOutliers =   0;
        Confocal_SubFcn.violinplot_advanced(StandBliss, 0, 0.5, properties)
        plot([0,0], [min(StandBliss), max(StandBliss)], 'k')
        clear properties
        % Swarm
        properties.MarkerType =         'o';
        properties.MarkerFaceColor =    SET.PhaseColors.(SET.phases{iPhase});
        properties.MarkerSize =         5;
        Confocal_SubFcn.beeswarmplot_advanced(StandBliss, 0.25, 0.5, properties)
        % Add mean and CI
        avg = mean(bootstrp(5000, @mean, StandBliss));
        ci = bootci(5000, @mean, StandBliss);
        plot([0.125 0.375], [avg avg], 'k')
        plot([0.25 0.25], ci, 'k')
        clear properties
        % --- Cosmetics
        xlim([-0.5 0.75])
        ylim([-1 1])
        ylabel('standardized Bliss score')
        % Statistics
        Stats.(SET.phases{iPhase}).StandBliss.seed=1234;
        Stats.(SET.phases{iPhase}).StandBliss.n = 5e4;
        Stats.(SET.phases{iPhase}).StandBliss.nComp = 1;
        [...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).StandBliss.p,...
            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).StandBliss.s] = Confocal_SubFcn.BootstrapHypothesisTesting('one-sample', StandBliss, 0, Stats.(SET.phases{iPhase}).StandBliss.n, Stats.(SET.phases{iPhase}).StandBliss.seed, Stats.(SET.phases{iPhase}).StandBliss.nComp);


        % Average PDF
        subplot(1,5,4); hold on
        % --- Get data
        avg = mean(bootstrp(5000, @mean, poolPdfBliss));
        ci = bootci(5000, @mean, poolPdfBliss);
        xvec = linspace(-100,100,1000);
        % --- Plot
        plot(xvec, avg, 'color', SET.PhaseColors.(SET.phases{iPhase}))
        plot(xvec, ci', 'color', SET.PhaseColors.(SET.phases{iPhase}))
        % --- Cosmetics
        xlim([-15 15])
        ylim(SET.Bliss_pdf_lim(iLayer,:))
        xlabel('Bliss score')
        ylabel('PDF')

        export_fig(['ConfocalAnalysis_Population', '\', SET.phases{iPhase},'_',SET.layers{iLayer},'_BSanalysis'], '-pdf')
        close(hFig)

        % Keep Bliss Score data for later to compare phases
        grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).BinaryBliss = poolBinaryBliss;
        grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgPosBliss = poolAvgPosBliss;
        grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgNegBliss = poolAvgNegBliss;
        grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).StandBliss = StandBliss;


        clearvars -except iPhase iLayer SET Stats grandMean_AVG
    end%iLayer
end%iPhase

hFig = figure;
cols{1}{1}=[1 0 1]; cols{1}{2}=[1 0 0];
cols{2}{1}=[0 1 1]; cols{2}{2}=[0 0 1];
for iPhase = 1:length(SET.phases)
    for iLayer = 1:length(SET.layers)

        % Get the data
        phi = rad2deg( bootstrp(50000, @mean, grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).StandBliss) * pi/2);
        r1 = grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgPosBliss;
        r1_boot = bootstrp(50000, @mean, r1);
        r2 = grandMean_AVG.(SET.phases{iPhase}).(SET.layers{iLayer}).AvgNegBliss;
        r2_boot = bootstrp(50000, @mean, r2);
        c1 = r1_boot.*[cosd(phi), sind(phi)];
        c2 = r2_boot.*[cosd(phi), sind(phi)];

        subplot(1,2,1); hold on
        % --- somata
        r_ellipse = Confocal_SubFcn.dataEllipse(c1, 0.95);
        plot(r_ellipse(:,1), r_ellipse(:,2), 'Color', cols{iPhase}{iLayer}, 'LineWidth', 1)
        % --- glomeruli
        r_ellipse = Confocal_SubFcn.dataEllipse(c2, 0.95);
        plot(r_ellipse(:,1), r_ellipse(:,2), 'Color', cols{iPhase}{iLayer}, 'LineWidth', 1)
        axis equal
        xlim([-3 5])
        ylim([-1.5 2.5])
        plot([mean(c1(:,1)), mean(c2(:,1))] ,  [mean(c1(:,2)), mean(c2(:,2))], 'k')
        title('Bliss map')

        subplot(1,2,2); hold on
        plot(r1, r2, 'o', 'MarkerFaceColor', cols{iPhase}{iLayer}, 'MarkerEdgeColor','none')
        % Fit model to data.
        [xData, yData] = prepareCurveData( r1, r2 );
        ft = fittype( 'poly1' );
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft );
        plot(xData, fitresult.p1*xData + fitresult.p2, 'Color', cols{iPhase}{iLayer})
        axis equal
        xlim([0 9])
        ylim([-5 0])
        text(8, fitresult.p1*xData(end) + fitresult.p2,...
            ['m=', num2str(round(fitresult.p1,2)),'|r2=',num2str(round(gof.rsquare,2))],...
            'Color', cols{iPhase}{iLayer})


    end%iLayer
end%iPhase
% Cosmetics
subplot(1,2,1); hold on
title('Bliss map')
subplot(1,2,2); hold on
title('Antagonism vs Synergy')
xlabel('Synergy')
ylabel('Antagonism')
% Save
export_fig('ConfocalAnalysis_Population\BSmap', '-pdf')
close(hFig)




% Stats comparing phases
Stats.gregar_vs_soli.seed=1234;
Stats.gregar_vs_soli.n = 5e4;
Stats.gregar_vs_soli.nComp = 1;
for iLayer = 1:length(SET.layers)
    % --- food
    [...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).Laa.p,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).Laa.s,...
        ~,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).Laa.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample', grandMean_AVG.gregarious.(SET.layers{iLayer}).ZHAE, grandMean_AVG.solitarious.(SET.layers{iLayer}).ZHAE, Stats.gregar_vs_soli.n, Stats.gregar_vs_soli.seed, Stats.gregar_vs_soli.nComp);
    % --- social
    [...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).Lct.p,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).Lct.s,...
        ~,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).Lct.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample', grandMean_AVG.gregarious.(SET.layers{iLayer}).COL, grandMean_AVG.solitarious.(SET.layers{iLayer}).COL, Stats.gregar_vs_soli.n, Stats.gregar_vs_soli.seed, Stats.gregar_vs_soli.nComp);
    % --- mix
    [...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).LaaLct.p,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).LaaLct.s,...
        ~,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).LaaLct.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample', grandMean_AVG.gregarious.(SET.layers{iLayer}).ZHAECOL, grandMean_AVG.solitarious.(SET.layers{iLayer}).ZHAECOL, Stats.gregar_vs_soli.n, Stats.gregar_vs_soli.seed, Stats.gregar_vs_soli.nComp);

    % --- Binary Bliss score
    [...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).BinaryBliss.p,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).BinaryBliss.s,...
        ~,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).BinaryBliss.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample',...
        grandMean_AVG.gregarious.(SET.layers{iLayer}).BinaryBliss,...
        grandMean_AVG.solitarious.(SET.layers{iLayer}).BinaryBliss,...
        Stats.gregar_vs_soli.n, Stats.gregar_vs_soli.seed, Stats.gregar_vs_soli.nComp);

    % --- Average pos. Bliss score
    [...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).AvgPosBliss.p,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).AvgPosBliss.s,...
        ~,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).AvgPosBliss.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample',...
        grandMean_AVG.gregarious.(SET.layers{iLayer}).AvgPosBliss,...
        grandMean_AVG.solitarious.(SET.layers{iLayer}).AvgPosBliss,...
        Stats.gregar_vs_soli.n, Stats.gregar_vs_soli.seed, Stats.gregar_vs_soli.nComp);

    % --- Average neg. Bliss score
    [...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).AvgNegBliss.p,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).AvgNegBliss.s,...
        ~,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).AvgNegBliss.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample',...
        grandMean_AVG.gregarious.(SET.layers{iLayer}).AvgNegBliss,...
        grandMean_AVG.solitarious.(SET.layers{iLayer}).AvgNegBliss,...
        Stats.gregar_vs_soli.n, Stats.gregar_vs_soli.seed, Stats.gregar_vs_soli.nComp);

    % --- Standardized Bliss score
    [...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).StandBliss.p,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).StandBliss.s,...
        ~,...
        Stats.gregar_vs_soli.(SET.layers{iLayer}).StandBliss.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample',...
        grandMean_AVG.gregarious.(SET.layers{iLayer}).StandBliss,...
        grandMean_AVG.solitarious.(SET.layers{iLayer}).StandBliss,...
        Stats.gregar_vs_soli.n, Stats.gregar_vs_soli.seed, Stats.gregar_vs_soli.nComp);
end%iLayer

% Save
save(['ConfocalAnalysis_Population', '\', 'statistics.mat'], 'Stats')





























