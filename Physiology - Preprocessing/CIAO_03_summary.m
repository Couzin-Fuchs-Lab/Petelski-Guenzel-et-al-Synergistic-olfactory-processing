%% Calcium Imaging Analysis Operator (summary)
% CIOA_03_summary is the third of a series of scripts aiming for
% an easy-to-use calcium imaging analysis pipeline. Here, get time courses
% for segmented data set and summarize the response to a given stimulus
% based on the average activity in a user-specified window.
% CIOA assumes a strict folder structer and each analysis step is saved in
% a given folder. Thus, each module can be run independently and results
% from different stages can be used without running a long pipeline.
%
% In the following, we will lay out the folder structure needed for CIAO to
% run without errors:
%   > DATA                       : this is the folder you select at the beginning
%   |--> Animal01                : this folder can have any name that contains the string "Animal01"
%   |  |--> 01_Data_raw          : CIAO assumes this folder contains the raw data, separated by trials
%   |  |  |--> Trial01           : CIAO assumes this folder contains the text file "protocol.txt" with meta inforamtion, together with individual tiff images
%   |  |  |--> Trial##           : same locig for all following trials
%   |  |--> 02_Data_Matlab       : folder will be created automatically and will contain pre-processed data
%   |  |  |--> Trial01           : each stimulus will be saved in the order it was given as a separate mat file
%   |  |  |--> Trial##           : same logic for all following trials
%   |  |--> 03_Data_Processed    : folder will be created automatically and will contain segmented data
%   |  |  |--> Trial01           : each stimulus will be saved in the order it was given as a separate mat file
%   |  |  |--> Trial##           : same logic for all following trials
%   |  |--> 04_Data_Summary      : folder will be created automatically and will contain segmented data
%   |  |  |--> Trial01           : each stimulus will be saved in the order it was given as a separate mat file
%   |  |  |--> Trial##           : same logic for all following trials
%   |--> Animal##                : same locic for all following animals
%
% Version: 16-Oct-22 (R2023a)

%% Clean and prepare environment
clc; clear all; close all
disp('----- Ciao, welcome back! -----')

% Add paths:
% --- Current folder and subfolders
addpath(genpath(pwd))

% Add toolboxes
% --- A MATLAB toolbox for exporting publication quality figures
%     (https://github.com/altmany/export_fig)
addpath(genpath('...\export_fig'))
% --- Toolbox for nice colors
%     https://github.com/DrosteEffect/BrewerMap
addpath(genpath('...\BrewerMap'))

%% Settings

% Set the correct fps of the recording
currSET.fps = 1;

% Set whether to run script for all animals, or just for new ones
currSET.overwrite = true;

% Set reference stiulus that should be subtracted
currSET.refStim = [];

% Set to which stimuli regions need to respond. Leave empty to take all
% currSET.StimRegionList = {'BERRYLEAF', 'COL', 'BERRYLEAF_COL'};
% currSET.StimRegionList = {'ZHAE', 'COL', 'ZHAECOL'};
% currSET.StimRegionList = {'COL', 'HEX3', 'HEX3_COL', 'OC3L2', 'OC3L2_COL', 'ZHAE2', 'ZHAE2_COL', 'ZHAE2_HEX3', 'ZHAE2_OC3L2'};
currSET.StimRegionList = [];


%% Iterate over all animals and their corresponding trials
% Get path to raw data
path2animal = uigetdir(pwd,'Select folder listing all animals');

% Set colormap
currSET.colMap = flipud(brewermap(1000,'RdBu'));

% Get overview of animal folders
curr.dir.all = dir(path2animal);

% Put everything in one figure
hFig_overview = figure('Units','normalized','Color','w', 'Position', [0 0 1 1]);
hFig_pool = figure('Units','normalized','Color','w', 'Position', [0 0 0.5 0.5]);

% Pool data
PooledData.StimList = {};

% Iterate over all animals
cnt_ani = 0;
for iAni = 1:size(curr.dir.all,1)

    % Check whether this is what we are looking for, i.e. whether the
    % string "Animal" is present
    if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir

        % Check whether to run preprocessing or not
        if currSET.overwrite || ~isfolder([path2animal,'\',curr.dir.all(iAni).name,'\04_Data_Summary'])
            cnt_ani = cnt_ani + 1;

            % Get information about trials
            curr.path.data = [path2animal,'\',curr.dir.all(iAni).name,'\03_Data_Processed'];

            % Check whether the processed folder is there
            if isfolder(curr.path.data)

                % Get meta info
                load([curr.path.data, '\meta_info.mat']);
                % Sort stimuli
                unique_stimuli = sort(meta.unique_stim);

                % Get info on segmentation
                load([curr.path.data, '\img_segmentation.mat'])

                % Prepare to save results
                curr.path.save = [path2animal,'\',curr.dir.all(iAni).name,'\04_Data_Summary'];
                mkdir(curr.path.save)

                % If not specified otherwise, use all stimuli
                if isempty(currSET.StimRegionList)
                    currSET.StimRegionList = meta.unique_stim;
                    currSET.clear_StimRegionList = 1;
                else
                    currSET.clear_StimRegionList = 0;
                end

                % Get indices of activity maps that should be included
                idx = zeros(1,length(currSET.StimRegionList));
                for iStim = 1:length(currSET.StimRegionList)
                    idx(iStim) = find(strcmp(unique_stimuli, currSET.StimRegionList{iStim}));
                end
                if currSET.clear_StimRegionList
                    currSET.StimRegionList = [];
                end

                % Get mask based on which pockets were active
                if isfield(Segmentation, 'pocket_active_img')
                    Mask = double(sum(Segmentation.pocket_active_img(:,:,idx),3)>0); %==length(idx)
                    pocket_list = unique(Segmentation.pockets_labeled.*Mask);
                else
                    if isfield(Segmentation, 'manualROI')
                        Mask = double(Segmentation.manualROI);
                        pocket_list = unique(Segmentation.pockets_labeled.*double(Segmentation.manualROI));
                    else
                        if isfield(Segmentation, 'autoROI')
                            Mask = double(Segmentation.autoROI);
                            pocket_list = unique(Segmentation.pockets_labeled.*double(Segmentation.autoROI));
                        else
                            Mask = double(Segmentation.pockets_bw);
                            pocket_list = unique(Segmentation.pockets_labeled.*double(Segmentation.pockets_bw));
                        end
                    end
                end
                Segmentation.pockets_bw = Mask;
                save([curr.path.data, '\img_segmentation.mat'], 'Segmentation')
                Mask(~Mask) = NaN;
                pocket_list(pocket_list==0) = [];

                % Get reday to plot everything
                figure(hFig_overview); clf
                figure(hFig_pool); clf
                cols = parula(length(unique_stimuli));

                % Iterate over stimul
                val_range.img_std = [];
                val_range.img_avg = [];
                val_range.tc = [];
                poolStd = zeros([size(meta.Mask),length(unique_stimuli)]);
                for iStim = 1:length(unique_stimuli)

                    % Load data
                    curr.data = load([curr.path.data, '\', unique_stimuli{iStim}, '.mat']);

                    % Load reference
                    if ~isempty(currSET.refStim)
                        curr.ref_data = load([curr.path.data, '\', currSET.refStim, '.mat']);
                        curr.ref_data = curr.ref_data.ImageStream.(currSET.refStim);
                    else
                        curr.ref_data = zeros(size(curr.data.ImageStream.(unique_stimuli{iStim})));
                    end

                    % Subtract reference
                    curr.data.ImageStream.(unique_stimuli{iStim}) = curr.data.ImageStream.(unique_stimuli{iStim})-curr.ref_data;

                    % Apply segmentation to keep active areas only
                    Summary.data.(unique_stimuli{iStim}).ImageStream = curr.data.ImageStream.(unique_stimuli{iStim}) .* Mask;
                    ImageStream_respahed  = reshape(...
                        Summary.data.(unique_stimuli{iStim}).ImageStream,...
                        [size(Summary.data.(unique_stimuli{iStim}).ImageStream,1)*size(Summary.data.(unique_stimuli{iStim}).ImageStream,2), size(Summary.data.(unique_stimuli{iStim}).ImageStream,3)]);

                    % Get the time course per pocket
                    % Iterate over all pockets and get their time courses
                    pocket_tc = zeros(length(pocket_list), size(ImageStream_respahed,2));
                    for iPockte = 1:length(pocket_list)
                        % Get time course
                        idx = find(Segmentation.pockets_labeled == pocket_list(iPockte));
                        pocket_tc(iPockte,:) = nanmean(ImageStream_respahed(idx,:));
                    end%iPocket

                    % Get grand mean
                    Summary.data.(unique_stimuli{iStim}).TC = nanmean(pocket_tc);
                    Summary.data.(unique_stimuli{iStim}).t = linspace(0, length(Summary.data.(unique_stimuli{iStim}).TC)/currSET.fps, length(Summary.data.(unique_stimuli{iStim}).TC));

                    % Get the average activity
                    Summary.data.(unique_stimuli{iStim}).avg_val = mean(Summary.data.(unique_stimuli{iStim}).TC);

                    % Get the mean and std-projection
                    Summary.data.(unique_stimuli{iStim}).img_avg = nanmean(Summary.data.(unique_stimuli{iStim}).ImageStream, 3);
                    Summary.data.(unique_stimuli{iStim}).img_std = nanstd(Summary.data.(unique_stimuli{iStim}).ImageStream, [], 3);

                    % Plot result
                    figure(hFig_overview)
                    % --- avg image
                    ax = subplot(3,length(unique_stimuli), iStim+length(unique_stimuli)*0);
                    imagesc(Summary.data.(unique_stimuli{iStim}).img_avg); axis equal off; colormap(ax, currSET.colMap)
                    title(['avg = ', num2str(round(Summary.data.(unique_stimuli{iStim}).avg_val,2)),'%'])
                    % --- std image
                    ax = subplot(3,length(unique_stimuli), iStim+length(unique_stimuli)*1);
                    imagesc(Summary.data.(unique_stimuli{iStim}).img_std); axis equal off; colormap(ax, gray(1000))
                    title('std')                    
                    % --- TC
                    subplot(3,length(unique_stimuli), iStim+length(unique_stimuli)*2)
                    plot(Summary.data.(unique_stimuli{iStim}).t, Summary.data.(unique_stimuli{iStim}).TC, 'LineWidth', 2, 'Color',	cols(iStim,:))
                    xlim([Summary.data.(unique_stimuli{iStim}).t(1), Summary.data.(unique_stimuli{iStim}).t(end)])
                    title(unique_stimuli{iStim}, 'interpreter', 'none')
                    xlabel('time (s)')
                    ylabel('\DeltaF/F (%)')
                    % --- range of plots
                    val_range.img_std = [val_range.img_std; nanmin(Summary.data.(unique_stimuli{iStim}).img_std(:)), nanmax(Summary.data.(unique_stimuli{iStim}).img_std(:))];
                    val_range.img_avg = [val_range.img_avg; nanmin(Summary.data.(unique_stimuli{iStim}).img_avg(:)), nanmax(Summary.data.(unique_stimuli{iStim}).img_avg(:))];
                    val_range.tc =      [val_range.tc;      nanmin(Summary.data.(unique_stimuli{iStim}).TC),         nanmax(Summary.data.(unique_stimuli{iStim}).TC)];

                    figure(hFig_pool)
                    subplot(2,2,[1 2]); hold on
                    plot(Summary.data.(unique_stimuli{iStim}).t, Summary.data.(unique_stimuli{iStim}).TC, 'LineWidth', 2, 'Color',	cols(iStim,:))
                    xlim([Summary.data.(unique_stimuli{iStim}).t(1), Summary.data.(unique_stimuli{iStim}).t(end)])
                    poolStd(:,:,iStim) = std(Summary.data.(unique_stimuli{iStim}).ImageStream,[],3);

                    % Pool results
                    if ~isfield(PooledData, unique_stimuli{iStim})
                        % Keep track of stimuli
                        PooledData.StimList{length(PooledData.StimList)+1, 1} = unique_stimuli{iStim};
                        PooledData.(unique_stimuli{iStim}) = nan(size(curr.dir.all,1),1);
                    end% if field alreday exists
                    PooledData.(unique_stimuli{iStim})(cnt_ani,1) = Summary.data.(unique_stimuli{iStim}).avg_val;
                end%iStim

                % Correct plots
                figure(hFig_overview)
                for iStim = 1:length(unique_stimuli)

                    % --- avg
                    subplot(3,length(unique_stimuli), iStim+length(unique_stimuli)*0); hold on
                    clim([-max(abs(val_range.img_avg(:))) max(abs(val_range.img_avg(:)))])
                    % Draw AL outlines
                    [B, L] = bwboundaries(Segmentation.manualROI>0);
                    for k = 1:length(B)
                        boundary = B{k};
                        plot(boundary(:, 2), boundary(:, 1), 'w:', 'LineWidth', 1)
                    end%k
                    % Draw area outlines
                    [B, L] = bwboundaries(Mask>0);
                    for k = 1:length(B)
                        boundary = B{k};
                        plot(boundary(:, 2), boundary(:, 1), 'w', 'LineWidth', 1)
                    end%k
                    clear k B L  boundary
                    % --- std
                    subplot(3,length(unique_stimuli), iStim+length(unique_stimuli)*1); hold on
                    caxis([min(val_range.img_std(:,1)), max(val_range.img_std(:,2))])
                    % Draw AL outlines
                    [B, L] = bwboundaries(Segmentation.manualROI>0);
                    for k = 1:length(B)
                        boundary = B{k};
                        plot(boundary(:, 2), boundary(:, 1), 'w:', 'LineWidth', 1)
                    end%k
                    % Draw area outlines
                    [B, L] = bwboundaries(Mask>0);
                    for k = 1:length(B)
                        boundary = B{k};
                        plot(boundary(:, 2), boundary(:, 1), 'w', 'LineWidth', 1)
                    end%k
                    clear k B L  boundary
                    % --- TC
                    subplot(3,length(unique_stimuli), iStim+length(unique_stimuli)*2); hold on
                    ylim([min(floor(val_range.tc(:,1) * 4) / 4), max(ceil(val_range.tc(:,2) * 4) / 4)])
                end%iStim
                % Add colorbars
                % --- avg
                subplot(3,length(unique_stimuli), iStim+length(unique_stimuli)*0);
                pos = get(gca, 'Position');
                colorbar
                set(gca, 'Position',  pos)
                % --- std
                subplot(3,length(unique_stimuli), iStim+length(unique_stimuli)*1);
                pos = get(gca, 'Position');
                colorbar
                set(gca, 'Position',  pos)                
                % Save image
                export_fig([curr.path.save, '\', curr.dir.all(iAni).name, '_overview'], '-pdf', '-painters')

                figure(hFig_pool)
                % --- TC
                subplot(2,2,[1 2]); hold on
                ylim([min(floor(val_range.tc(:,1) * 4) / 4), max(ceil(val_range.tc(:,2) * 4) / 4)])
                xlim([Summary.data.(unique_stimuli{iStim}).t(1), Summary.data.(unique_stimuli{iStim}).t(end)])
                legend(unique_stimuli, 'Location', 'eastoutside', 'Interpreter','none')
                title(curr.dir.all(iAni).name, 'interpreter', 'none')
                xlabel('time (s)')
                ylabel('\DeltaF/F (%)')
                % --- avg std
                subplot(2,2,3); hold on
                imagesc(nanmean(poolStd,3)); axis equal off;
                caxis([nanmin([quantile(poolStd(:),0.01),-1e-6]), nanmax([quantile(poolStd(:),0.99),+1e-6])]); clear poolStd
                title('average std-projection')
                % Draw AL outlines
                [B, L] = bwboundaries(Segmentation.manualROI>0);
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:, 2), boundary(:, 1), 'w:', 'LineWidth', 1)
                end%k
                % Drar region outlines
                [B, L] = bwboundaries(Mask>0);
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:, 2), boundary(:, 1), 'w', 'LineWidth', 1)
                end%k
                set(gca, 'YDir', 'reverse')
                % --- avg in z score
                subplot(2,2,4); hold on
                temp = Segmentation.summary_stats.pocket_Corr_img;
                imagesc(temp); axis equal off; 
                colormap gray
                clim([-1 1]); clear temp
                title('corr')
                % Draw AL outlines
                [B, L] = bwboundaries(Segmentation.manualROI>0);
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:, 2), boundary(:, 1), 'w:', 'LineWidth', 1)
                end%k
                % Drar region outlines
                [B, L] = bwboundaries(Mask>0);
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:, 2), boundary(:, 1), 'w', 'LineWidth', 1)
                end%k
                set(gca, 'YDir', 'reverse')
                colorbar
                export_fig([curr.path.save, '\', curr.dir.all(iAni).name, '_pooled'], '-pdf', '-painters')

                % Put all we need in the summary variable
                Summary.meta = meta;
                Summary.stimuli = unique_stimuli;
                Summary.Segmentation = Segmentation;
                Summary.fps = currSET.fps;

                % Save
                save([curr.path.save, '\', curr.dir.all(iAni).name, '.mat'], 'Summary')
                clear Summary meta Segmentation

            end%if reafy
        end%if overwrite
    end%if animal folder
end%iAni

% % Plot pooled results across all stimuli
% hFig_poolAll = figure('Units','normalized','Color','w', 'Position', [0 0 0.5 0.5]); hold on
% cols = turbo(length(PooledData.StimList));
% Sample = [];
% for iStim = 1:length(PooledData.StimList)
%     % Get data
%     BeeData = PooledData.(PooledData.StimList{iStim});
%     % Properties for swarm plot
%     properties.MarkerType =         'o';                 %(Marker type)
%     properties.MarkerFaceColor =    cols(iStim,:);       %(Marker face color)
%     properties.MarkerSize =         5;                   %(Marker size)
%     % Plot results
%     if length(BeeData)>1
%         % Swarm plot
%         SubFcn.beeswarmplot_advanced(BeeData, iStim, 0.75, properties)
%         % Avg and CI
%         avg = nanmean(bootstrp(5000, @nanmean, BeeData));
%         CI = bootci(5000, {@nanmean, BeeData});
%         plot([iStim iStim]+[-0.4 0.4], [avg avg], 'k', 'LineWidth', 2)
%         plot([iStim iStim], CI, 'k', 'LineWidth', 2)
%     else
%         plot(iStim, BeeData, 'o', 'MarkerFaceColor', cols(iStim,:), 'MarkerEdgeColor', 'none')
%         plot([iStim iStim]+[-0.4 0.4],[BeeData BeeData],'k', 'LineWidth', 2)
%     end
% 
% end%iStim
% % Cosmetics
% set(gca, 'XTick', 1:length(PooledData.StimList), 'XTickLabel', PooledData.StimList, 'TickLabelInterpreter', 'none')
% ylabel('\DeltaF/F (%)')
% xlabel('stimulus')
% title(path2animal)
% export_fig([path2animal, 'PooledData'], '-pdf', '-painters')

% Save
save([path2animal, 'PooledData.mat'], 'PooledData')

%% Clean and prepare environment
clc; clear all; close all
disp('----- Ciao, see you next time! -----')