%% ConfocalAnalysis_ExampleAnimals
% This script analyzes four exemplary animals, two gregarious and two
% solitarious animals. For each phase, one recording is from the cell body
% layer, and one is from the glomerulus layer. It produces heatmap
% (s.d. projection and mean intensity projection) and time course figures
% for each stimulus.
%
% Version:
% 05-Jan-2024 (R2023a)

% Prepare
clc; clear all; close all
warning('off')

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig)
addpath(genpath('C:\Users\Yannick\Documents\GitHub\export_fig'))
mkdir('ConfocalAnalysis_ExampleAnimals')

%% Settings

% Provide a list of example animals (selected based on BF)
SET.example_animals.gregarious = {...
    '...\gregarious\layer01_top\220411_Animal029_greg_d050\';...       top layer
    '...\gregarious\layer02_middle\220926_Animal055_greg_d160\'};%     middle layer

SET.example_animals.solitarious = {...
    '...\solitarious\layer01_top\221027_Animal065_soli_d095\';...      top layer
    '...\solitarious\layer02_middle\220316_Animal019_soli_d140\'};%    middle layer

% Give list of phases
SET.phases = {'gregarious'; 'solitarious'};

% Give layer names
SET.layers = {'layer01_top', 'layer02_middle'};

% Set framerate
SET.fps = 1;

% Set how to cut the data and which window for analysis should be used
SET.totalFrameNumber = 17;
SET.AnalysisWindow = [...
    10 11;... layer 01
    8 9;...   layer 02
    ];

% Give list of stimuli (We have this in the meta data, but good to double
% check)
SET.StimList = {...
    'COL';...
    'ZHAE';...
    'ZHAECOL'};

% Cosmetics
SET.Colors = [...
    217,095,002;...Col
    102,166,030;...ZHAE
    231,041,138;...ZZHAECOL
    ]/255;

% Colormap for heatmaps
SET.HeatmapColors.gregarious = Confocal_SubFcn.ColMapPlasma(1000);
SET.HeatmapColors.solitarious = [...
    SET.HeatmapColors.gregarious(:,2),...
    SET.HeatmapColors.gregarious(:,1),...
    SET.HeatmapColors.gregarious(:,3)];
SET.HeatmapColors.other = gray(1000);

%% Create figures

% Iterate over both phases
for iPhase = 1:length(SET.phases)

    % Iterate over the three examples from the current phase
    for iAni = 1:length(SET.example_animals.(SET.phases{iPhase}))

        % Get the phase of the current animal
        curr.phase = SET.phases{iPhase};
        % Get meta data and segmentation
        load([SET.example_animals.(SET.phases{iPhase}){iAni},'03_Data_Processed\meta_info.mat'])
        load([SET.example_animals.(SET.phases{iPhase}){iAni},'03_Data_Processed\img_segmentation.mat'])

        % Get list of valid pockets. Take all pockets that were active
        % at least once.
        mask_list = zeros(length(meta.unique_stim),1);
        for iStim = 1:length(SET.StimList)
            mask_list = mask_list+ (meta.unique_stim == SET.StimList{iStim});
        end
        mask = sum(Segmentation.pocket_active_img(:,:,find(mask_list)),3);
        mask = mask>0;
        pocket_list = unique(Segmentation.pockets_labeled);
        pocket_list(pocket_list==0) = [];
        clear mask_list iStim

        % Get the BF image, the avg time courses for each stimulus, and the
        % average heat map for each stimulus.
        BF = [];
        HM_lim = [];
        TC_lim = [];
        % Iterate over all stimuli
        for iStim = 1:length(SET.StimList)
            % BF image
            currData = load([SET.example_animals.(SET.phases{iPhase}){iAni},'03_Data_Processed\',SET.StimList{iStim},'.mat']);
            currData = currData.ImageStream.(SET.StimList{iStim});
            [x,y,t] = size(currData);
            BF = cat(3,BF,currData);
            % Time course with confidence interval
            % --- reshape data
            currData = reshape(currData, [x*y, t]);
            % --- cut data
            time_vec = meta.baseline(1):(meta.baseline(1)+SET.totalFrameNumber-1);
            currData = currData(:, time_vec);
            % --- iterate over all pockets to create the heatmap
            HM.(SET.StimList{iStim}) = nan(size(Segmentation.pockets_labeled));
            for iP = 1:length(pocket_list)
                idx = find(Segmentation.pockets_labeled == pocket_list(iP));
                HM.(SET.StimList{iStim})(idx) = nanmean(nanmean(currData(idx, SET.AnalysisWindow(iAni,1):SET.AnalysisWindow(iAni,end) )));
            end%iP
            HM_lim = [HM_lim; min(HM.(SET.StimList{iStim})(mask==1)), max(HM.(SET.StimList{iStim})(mask==1))];
            % --- iterate over all valid pockets
            TC = nan(length(pocket_list), length(time_vec));
            for iP = 1:length(pocket_list)
                idx = find(Segmentation.pockets_labeled == pocket_list(iP));
                TC(iP,:) = mean(currData(idx,:));
            end%ip
            % --- mean and CI
            avgTC.(SET.StimList{iStim}) = [
                nanmean(bootstrp(5000, @mean, TC));...
                bootci(5000, @mean, TC)];
            TC_lim = [TC_lim; min(avgTC.(SET.StimList{iStim})(:)), max(avgTC.(SET.StimList{iStim})(:))];
        end%iStim

        % Get the BF image as the std-projection across all frames
        hFig_BF = figure('Color','w');
        BF = nanstd(BF,[],3);
        % Normalize
        BF = BF - nanmin(BF(:));
        BF = sqrt(BF);
        BF = BF - nanmean(BF(:));
        BF = BF / nanstd(BF(:));

        % Get range in ROI
        imagesc(BF)
        temp = BF;
        temp(~Segmentation.manualROI) = NaN;
        clim([min(temp(:)), max(temp(:))])
        colormap gray
        axis equal off
        title('BF')
        export_fig(['ConfocalAnalysis_ExampleAnimals', '\', curr.phase,'_',SET.layers{iAni}, '_BF'], '-pdf', '-painters')
        close(hFig_BF)

        % Plot everything
        hFig_TC = figure('Color','w');
        % indicate active region
        xvec = 0 : (1/SET.fps) : (size(TC,2) - (1/SET.fps));
        rectangle('position', [xvec(SET.AnalysisWindow(iAni,1)), floor(min(TC_lim(:,1))*2)/2, xvec(SET.AnalysisWindow(iAni,end))-xvec(SET.AnalysisWindow(iAni,1)), ceil(max(TC_lim(:,2))*2)/2-floor(min(TC_lim(:,1))*2)/2],...
            'FaceColor',[0.5 0.5 0.5],...
            'EdgeColor','none')
        % Iterate over all stimuli and populate the plots
        for iStim = 1:length(SET.StimList)
            % Heat map
            hFig_HM = figure('Color','w');
            % In the background, plot everything in gray scale
            ax1 = axes;
            hB = imagesc(HM.(SET.StimList{iStim})); axis equal off;
            clim([floor(max(HM_lim(:,1))*2)/2, ceil(min(HM_lim(:,2))*2)/2])
            % In the foreground plot everything valid in color
            ax2 = axes;
            hF = imagesc(HM.(SET.StimList{iStim})); axis equal off;
            clim([floor(max(HM_lim(:,1))*2)/2, ceil(min(HM_lim(:,2))*2)/2])
            set(hF,'AlphaData',mask);
            % Cosmetics
            colormap(ax1,SET.HeatmapColors.other)
            colormap(ax2,SET.HeatmapColors.(curr.phase))
            cb = colorbar(ax1); cb.Ticks = [floor(max(HM_lim(:,1))*2)/2, ceil(min(HM_lim(:,2))*2)/2];
            cb = colorbar(ax2); cb.Ticks = [floor(max(HM_lim(:,1))*2)/2, ceil(min(HM_lim(:,2))*2)/2];
            title(SET.StimList{iStim}, 'Interpreter','none')
            print(['ConfocalAnalysis_ExampleAnimals', '\', curr.phase,'_',SET.layers{iAni}, '_heatmap_', SET.StimList{iStim}], '-dpdf')
            close(hFig_HM)
            % Time course
            figure(hFig_TC)
            hold on
            % --- Plot stim response
            xvec = linspace(0, size(TC,2)/SET.fps, size(TC,2));
            plot(xvec, avgTC.(SET.StimList{iStim})', 'color', SET.Colors(iStim,:))
            % --- cosmetics
            xlim([xvec(1), xvec(end)])
            ylim([floor(min(TC_lim(:,1))*2)/2, ceil(max(TC_lim(:,2))*2)/2])
            ylabel('\DeltaF/F (%)')
            xlabel('time (s)')
            clear ax1 hB ax2 hF B L k xvec
        end%iStim
        figure(hFig_TC)
        export_fig(['ConfocalAnalysis_ExampleAnimals', '\', curr.phase,'_',SET.layers{iAni}, '_TC'], '-pdf', '-painters')
        close(hFig_TC)
        clearvars -except SET iAni iPhase
    end%iAni
end%iPhase