%% BerryLeafAnalysis_ExampleAnimals
% ...
%
% Version:
% 24-April-2023 (R2023a) Yannick Günzel

% Prepare
clc; clear all; close all
warning('off')

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig)
addpath(genpath('...\GitHub\export_fig'))
mkdir('BerryLeafAnalysis_ExampleAnimals')

%% Settings

% Provide a list of example animals (selected based on BF)
SET.example_animals.gregarious = {...
    '...\Data\Physiology\Cal520_Widefield\BERRYLEAF_COL_MOL_N2_ZHAE2\gregarious\230215_Animal253_greg_socialmodulation\'};

SET.example_animals.solitarious = {...
    '...\Data\Physiology\Cal520_Widefield\BERRYLEAF_COL_MOL_N2_ZHAE2\solitarious\230418_Animal273_soli_socialmodulation\'};

% Give list of phases
SET.phases = {'gregarious'; 'solitarious'};

% Set framerate
SET.fps = 10;

% Give list of stimuli (We have this in the meta data, but good to double
% check)
SET.StimList = {...
    'N2';...
    'BERRYLEAF';...
    'COL';...
    'BERRYLEAF_COL'};

% Cosmetics
SET.Colors = [...
    200,200,200;...N2
    102,166,030;...BERRYLEAF
    217,095,002;...Col
    231,041,138;...ZBERRYLEAF_COL
    ]/255;

% Colormap for heatmaps
SET.HeatmapColors.gregarious = Widefield_SubFcn.ColMapPlasma(1000);
SET.HeatmapColors.solitarious = [...
    SET.HeatmapColors.gregarious(:,2),...
    SET.HeatmapColors.gregarious(:,1),...
    SET.HeatmapColors.gregarious(:,3)];
SET.HeatmapColors.other = gray(1000);

SET.ylim_range = [-0.5 2];


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
        % Get list of pockets. Take all pockets that were active
        % at least once
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
        for iStim = 1:length(SET.StimList)
            % BF image
            currData = load([SET.example_animals.(SET.phases{iPhase}){iAni},'03_Data_Processed\',SET.StimList{iStim},'.mat']);
            BF = cat(3,BF,currData.ImageStream.(SET.StimList{iStim}));
            % Time course with confidence interval
            % --- reshape data
            currData = reshape(currData.ImageStream.(SET.StimList{iStim}), [size(currData.ImageStream.(SET.StimList{iStim}),1)*size(currData.ImageStream.(SET.StimList{iStim}),2), size(currData.ImageStream.(SET.StimList{iStim}),3)]);
            % --- iterate over all pockets to create the heatmap
            HM.(SET.StimList{iStim}) = nan(size(Segmentation.pockets));
            for iP = 1:length(pocket_list)
                idx = find(Segmentation.pockets_labeled == pocket_list(iP));
                HM.(SET.StimList{iStim})(idx) = mean(mean(currData(idx,meta.activeRegion)));
            end%iP
            HM_lim = [HM_lim; min(HM.(SET.StimList{iStim})(mask==1)), max(HM.(SET.StimList{iStim})(mask==1))];
            % --- get time vector
            time_vec = meta.baseline(1):meta.stimulus_onset(2)-1;
            % --- get indices of active region
            active_region = [find(time_vec==meta.activeRegion(1)), find(time_vec==meta.activeRegion(end))];
            % --- iterate over all valid pockets
            TC = nan(length(pocket_list), length(time_vec));
            for iP = 1:length(pocket_list)
                idx = find(Segmentation.pockets_labeled == pocket_list(iP));
                TC(iP,:) = mean(currData(idx,time_vec));
            end%ip
            % --- mean and CI
            avgTC.(SET.StimList{iStim}) = [
                nanmean(bootstrp(5000, @mean, TC));...
                bootci(5000, @mean, TC)];
            TC_lim = [TC_lim; min(avgTC.(SET.StimList{iStim})(:)), max(avgTC.(SET.StimList{iStim})(:))];
        end%iStim
        % Get the BF image as the max-projection across all frames
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
        export_fig(['BerryLeafAnalysis_ExampleAnimals', '\', curr.phase,'_BF'],'-pdf')
        close(hFig_BF)

        % Plot everything
        hFig_TC = figure('Color','w');
        % indicate active region
        xvec = linspace(0, size(TC,2)/SET.fps, size(TC,2));
        rectangle('position', [xvec(active_region(1)), floor(min(TC_lim(:,1))*2)/2, xvec(active_region(2))-xvec(active_region(1)), ceil(max(TC_lim(:,2))*2)/2-floor(min(TC_lim(:,1))*2)/2],...
            'FaceColor',[0.5 0.5 0.5],...
            'EdgeColor','none')
        % Iterate over all stimuli and populate the plots
        for iStim = 1:length(SET.StimList)
            % Heat map
            hFig_HM = figure('Color','w');
            % In the background, plot everything in gray scale
            ax1 = axes;
            hB = imagesc(HM.(SET.StimList{iStim})); axis equal off;
            clim([floor(min(HM_lim(2:end,1))*2)/2, ceil(max(HM_lim(2:end,2))*2)/2])
            % In the foreground plot everything valid in color
            ax2 = axes;
            hF = imagesc(HM.(SET.StimList{iStim})); axis equal off;
            clim([floor(min(HM_lim(2:end,1))*2)/2, ceil(max(HM_lim(2:end,2))*2)/2])
            set(hF,'AlphaData',mask);
            % Cosmetics
            colormap(ax1,SET.HeatmapColors.other)
            colormap(ax2,SET.HeatmapColors.(curr.phase))
            cb = colorbar(ax1); cb.Ticks = [floor(min(HM_lim(2:end,1))*2)/2, ceil(max(HM_lim(2:end,2))*2)/2];
            cb = colorbar(ax2); cb.Ticks = [floor(min(HM_lim(2:end,1))*2)/2, ceil(max(HM_lim(2:end,2))*2)/2];
            title(SET.StimList{iStim}, 'Interpreter','none')
            print(['BerryLeafAnalysis_ExampleAnimals', '\', curr.phase,'_heatmap_', SET.StimList{iStim}],'-dpdf')
            close(hFig_HM)

            % Time course
            figure(hFig_TC)
            hold on
            % --- Plot stim response
            xvec = linspace(0, size(TC,2)/SET.fps, size(TC,2));
            plot(xvec, avgTC.(SET.StimList{iStim})', 'color', SET.Colors(iStim,:))
            % --- cosmetics
            xlim([xvec(1), xvec(end)])
            ylim(SET.ylim_range)
            ylabel('\DeltaF/F (%)')
            xlabel('time (s)')

            clear ax1 hB ax2 hF B L k xvec
        end%iStim
        figure(hFig_TC)
        export_fig(['BerryLeafAnalysis_ExampleAnimals', '\', curr.phase,'_TC'],'-pdf')
        close(hFig_TC)
        clearvars -except SET iAni iPhase
    end%iAni
end%iPhase












