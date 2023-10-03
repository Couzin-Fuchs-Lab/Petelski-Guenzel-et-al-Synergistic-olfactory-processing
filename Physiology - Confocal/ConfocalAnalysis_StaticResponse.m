%% ConfocalAnalysis_StaticResponse
% ...
%
% Version:
% 28-July-2023 (R2023a) Yannick Günzel

% Prepare
clc; clear all; close all
warning('off')

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig)
addpath(genpath('...\GitHub\export_fig'))
mkdir('ConfocalAnalysis_StaticResponse')

%% Settings

% Give paths to find the data
SET.basepath = '...\Data\Physiology\Cal520_Confocal\static_response\';
SET.layers = {'layer01_top'; 'layer02_middle'};

% Give stimulus names
SET.stim_list = {'OCT', 'ZHAE', 'HEX'};

% PCA settings
SET.exp = 95;

% Set framerate
SET.fps = 1;

% Colormap for heatmaps
SET.HeatmapColors = Confocal_SubFcn.ColMapPlasma(1000);

% Cosmetics
SET.Colors = [...
    117,112,179;...HEX3
    102,166,030;...ZHAE2
    230,171,002;...OC3L2
    ]/255;

%% Get data

% Iterate over all data sets, pool the data, and compute within stimulus
% and across stimuli differences
for iLayer = 1:length(SET.layers)

    % Get overview of animal folders
    path2animal = [SET.basepath,SET.layers{iLayer},'\'];
    curr.dir.all = dir(path2animal);

    % Iterate over all animals
    cnt_ani = 1;
    for iAni = 1:size(curr.dir.all,1)
        % Check whether this is what we are looking for, i.e. whether the
        % string "Animal" is present
        if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir

            % Get the path to the current animal
            curr.dir.animal = [path2animal,curr.dir.all(iAni).name,'\'];

            % Get meta data and segmentation
            load([curr.dir.animal,'03_Data_Processed\meta_info.mat'])
            load([curr.dir.animal,'03_Data_Processed\img_segmentation.mat'])

            % Iterate over all stimuli and get the corresponding std
            % projection and the avg projection during the active window
            for iStim = 1:length(SET.stim_list)
                % Prepare some data structures
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection = [];
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection = [];
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection_mat = [];
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection_mat = [];
                % Make it a little bit more complicated than just using the
                % meta file to ensure a specific stimulus get assigned to a
                % certain color in the RGB overlay.
                ind_stim = find(contains(meta.unique_stim, SET.stim_list{iStim}));
                for jStim = 1:length(ind_stim)
                    % Load the data
                    curr.stim = meta.unique_stim{ind_stim(jStim)};
                    curr.data = load([curr.dir.animal, '03_Data_Processed\', curr.stim, '.mat']);
                    % Get avg and std intensity projection
                    avg_projection = nanmean(curr.data.ImageStream.(curr.stim)(:,:,meta.activeRegion),3);
                    std_projection = nanstd(curr.data.ImageStream.(curr.stim),[],3);
                    % Pool based on segmentation
                    pocket_list = unique(Segmentation.pockets_labeled(:));
                    for iP = 1:length(pocket_list)
                        avg_projection(Segmentation.pockets_labeled == pocket_list(iP)) = nanmean(nanmean(avg_projection(Segmentation.pockets_labeled == pocket_list(iP))));
                        std_projection(Segmentation.pockets_labeled == pocket_list(iP)) = nanmean(nanmean(std_projection(Segmentation.pockets_labeled == pocket_list(iP))));
                    end
                    % Apply manual selection of the AL
                    avg_projection(~Segmentation.manualROI) = NaN;
                    std_projection(~Segmentation.manualROI) = NaN;
                    % Cut data to the smallest area possible
                    linearIndices = find(Segmentation.manualROI);
                    [x, y] = ind2sub(size(Segmentation.manualROI), linearIndices);
                    avg_projection = avg_projection(min(x):max(x), min(y):max(y));
                    std_projection = std_projection(min(x):max(x), min(y):max(y));
                    % Normalize to have values between zero and one
                    avg_projection = avg_projection-nanmin(avg_projection(:));
                    std_projection = std_projection-nanmin(std_projection(:));
                    avg_projection = avg_projection/nanmax(avg_projection(:));
                    std_projection = std_projection/nanmax(std_projection(:));
                    % Concatenate
                    % -- as an image
                    PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection = ...
                        cat(3, PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection, avg_projection);
                    PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection = ...
                        cat(3, PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection, std_projection);
                    % As row vectors
                    avg_projection(isnan(avg_projection)) = [];
                    std_projection(isnan(std_projection)) = [];
                    PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection_mat = ...
                        [PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection_mat; avg_projection(:)'];
                    PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection_mat = ...
                        [PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection_mat; std_projection(:)'];
                end%jStim
                % Get the average across all repetitions
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_avg_projection = ...
                    nanmean(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection, 3);
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_std_projection = ...
                    nanmean(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection, 3);
                % Normalize to have values between zero and one
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_avg_projection = PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_avg_projection-nanmin(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_avg_projection(:));
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_std_projection = PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_std_projection-nanmin(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_std_projection(:));
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_avg_projection = PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_avg_projection/nanmax(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_avg_projection(:));
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_std_projection = PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_std_projection/nanmax(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_std_projection(:));
            end%iStim
            % Get the RGB overlay
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).avg_RGB = cat(3,...
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{1}).mean_avg_projection,...
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{2}).mean_avg_projection,...
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{3}).mean_avg_projection);
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).std_RGB = cat(3,...
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{1}).mean_std_projection,...
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{2}).mean_std_projection,...
                PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{3}).mean_std_projection);

            %% Now, get the within and across difference between stimuli
            % First pool
            all_avg = [];
            all_std = [];
            all_label = [];
            for iStim = 1:length(SET.stim_list)
                all_avg = [all_avg; PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection_mat];
                all_std = [all_std; PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection_mat];
                all_label = [all_label; iStim*ones(size(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection_mat,1),1)];
            end%iStim
            % Normalize. But not the range but zscore to ensure that one
            % class does not have a larger mean
            all_avg = normalize(all_avg,2);
            all_std = normalize(all_std,2);
            % Apply PCA to reduce the dimensions
            % --- avg projection
            [~,all_avg_pca,~,~,avg_explained,~] = pca(all_avg, 'Economy', false);
            ind = find(cumsum(avg_explained)>SET.exp , 1);
            all_avg_pca = all_avg_pca(:,1:ind);
            % --- std projection
            [~,all_std_pca,~,~,std_explained,~] = pca(all_std, 'Economy', false);
            ind = find(cumsum(std_explained)>SET.exp , 1);
            all_std_pca = all_std_pca(:,1:ind);
            % Get the average distances
            D_avg = squareform(pdist(all_avg_pca, "euclidean"));
            D_std = squareform(pdist(all_std_pca, "euclidean"));
            % Iterate over all distances and pool them within/across
            % stimuli
            D_tril_avg = tril(D_avg);
            D_tril_std = tril(D_std);
            dist_within_avg = [];
            dist_within_std = [];
            dist_across_avg = [];
            dist_across_std = [];
            % Iterate over all stimuli
            for iStim = 1:length(all_label)
                for jStim = 1:length(all_label)
                    if iStim~=jStim && D_tril_avg(iStim, jStim)~=0
                        % Check whether we compare within the same stimulus
                        % or across different ones
                        if all_label(iStim) == all_label(jStim)
                            dist_within_avg = [dist_within_avg; D_tril_avg(iStim, jStim)];
                            dist_within_std = [dist_within_std; D_tril_std(iStim, jStim)];
                        else
                            dist_across_avg = [dist_across_avg; D_tril_avg(iStim, jStim)];
                            dist_across_std = [dist_across_std; D_tril_std(iStim, jStim)];
                        end%if within
                    end%if valid
                end%jStim
            end%iStim
            % Pool
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.label = all_label;
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.avg_explained = avg_explained;
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.std_explained = std_explained;
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.all_avg_pca = all_avg_pca;
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.all_std_pca = all_std_pca;
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.dist_within_avg = dist_within_avg;
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.dist_within_std = dist_within_std;
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.dist_across_avg = dist_across_avg;
            PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.dist_across_std = dist_across_std;


            %% Plot heatmaps
            hFig = figure('Units','normalized','Position',[0 0 1 1],'Color','w');
            
            % Average projection
            % ----- Plot the grand mean across everything
            ax = subplot(2,5,1);
            imagesc(nanmean(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).avg_RGB, 3))
            axis equal tight off
            colormap(ax, gray(1000))
            title('mean avg proje.')
            % ----- Plot the projection for each stimulus
            for iStim = 1:length(SET.stim_list)
                ax = subplot(2,5,1+iStim);
                imagesc(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_avg_projection)
                cmap = zeros(1000,3); cmap(:,iStim) = linspace(0,1,1000)';
                axis equal tight off
                colormap(ax, cmap)
                clim([min(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).avg_RGB(:)), max(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).avg_RGB(:))])
                title(['avg proje. ', SET.stim_list{iStim}])
            end
            % ----- Plot RGB overlay
            subplot(2,5,5);
            imagesc(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).avg_RGB)
            axis equal tight off
            clim([min(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).avg_RGB(:)), max(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).avg_RGB(:))])
            title('RGB overlay')
           
            
            % SD projection
            % ----- Plot the grand mean across everything
            ax = subplot(2,5,6);
            imagesc(nanmean(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).std_RGB, 3))
            axis equal tight off
            colormap(ax, gray(1000))
            title('mean std proje.')
            % ----- Plot the projection for each stimulus
            for iStim = 1:length(SET.stim_list)
                ax = subplot(2,5,6+iStim);
                imagesc(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).mean_std_projection)
                cmap = zeros(1000,3); cmap(:,iStim) = linspace(0,1,1000)';
                axis equal tight off
                colormap(ax, cmap)
                clim([min(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).std_RGB(:)), max(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).std_RGB(:))])
                title(['std proje. ', SET.stim_list{iStim}])
            end
            % ----- Plot RGB overlay
            subplot(2,5,10);
            imagesc(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).std_RGB)
            axis equal tight off
            clim([min(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).std_RGB(:)), max(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).std_RGB(:))])
            title('RGB overlay')
            % Save
            export_fig([curr.dir.animal, '04_Data_Summary\static_heatmap'], '-pdf', '-painters')
            close(hFig)


            %% Plot all heatmaps
            hFig = figure('Units','normalized','Position',[0 0 1 1],'Color','w');
            % Get the maximum number of reps
            max_rep = [];
            for iStim = 1:length(SET.stim_list)
                max_rep = [max_rep; size(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection, 3)];
            end%iStim
            % All reps for the the current animal
            cnt_all = 1;
            clim_all = [];
            for iStim = 1:length(SET.stim_list)
                for iRep = 1:size(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection, 3)
                    ax = subplot(length(SET.stim_list), max(max_rep)+1, cnt_all);
                    imagesc(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection(:,:,iRep))
                    axis equal tight off
                    cmap = zeros(1000,3); cmap(:,iStim) = linspace(0,1,1000)';
                    colormap(ax, cmap)
                    title(iRep)
                    cnt_all = cnt_all+1;
                end%iRep
                ax = subplot(length(SET.stim_list), max(max_rep)+1, cnt_all);
                imagesc(mean(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).avg_projection,3))
                axis equal tight off
                cmap = zeros(1000,3); cmap(:,iStim) = linspace(0,1,1000)';
                colormap(ax, cmap)
                title('avg avg')
                cnt_all = cnt_all+1;
            end%iStim
            % Save
            export_fig([curr.dir.animal, '04_Data_Summary\all_avg_reps_heatmaps'], '-pdf', '-painters')
            close(hFig)

            %% Plot all heatmaps
            hFig = figure('Units','normalized','Position',[0 0 1 1],'Color','w');
            % Get the maximum number of reps
            max_rep = [];
            for iStim = 1:length(SET.stim_list)
                max_rep = [max_rep; size(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection, 3)];
            end%iStim
            % All reps for the the current animal
            cnt_all = 1;
            clim_all = [];
            for iStim = 1:length(SET.stim_list)
                for iRep = 1:size(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection, 3)
                    ax = subplot(length(SET.stim_list), max(max_rep)+1, cnt_all);
                    imagesc(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection(:,:,iRep))
                    axis equal tight off
                    cmap = zeros(1000,3); cmap(:,iStim) = linspace(0,1,1000)';
                    colormap(ax, cmap)
                    title(iRep)
                    cnt_all = cnt_all+1;
                end%iRep
                ax = subplot(length(SET.stim_list), max(max_rep)+1, cnt_all);
                imagesc(mean(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).(SET.stim_list{iStim}).std_projection,3))
                axis equal tight off
                cmap = zeros(1000,3); cmap(:,iStim) = linspace(0,1,1000)';
                colormap(ax, cmap)
                title('avg std')
                cnt_all = cnt_all+1;
            end%iStim
            % Save
            export_fig([curr.dir.animal, '04_Data_Summary\all_std_reps_heatmaps'], '-pdf', '-painters')
            close(hFig)
            

            %% PCA space
            hFig = figure('Units','normalized','Position',[0.25 0.25 0.5 0.5],'Color','w');
            for iStim = 1:length(SET.stim_list)
                ind = find(all_label == iStim);
                subplot(1,2,1); hold on
                plot(all_avg_pca(ind,1), all_avg_pca(ind,2), 'o', 'MarkerFaceColor', SET.Colors(iStim,:), 'MarkerEdgeColor','none', 'MarkerSize', 10)
                % Add confidence ellipse
                dat_pca_boot = bootstrp(5000, @nanmean, all_avg_pca(ind,1:2));
                % Get a 99%-confidence ellipse around the bootstrapped data
                r_ellipse = Confocal_SubFcn.dataEllipse(dat_pca_boot, 0.99);
                plot(r_ellipse(:,1), r_ellipse(:,2), 'Color', [0.251 0.251 0.251], 'LineWidth', 1)                
                subplot(1,2,2); hold on
                plot(all_std_pca(ind,1), all_std_pca(ind,2), 'o', 'MarkerFaceColor', SET.Colors(iStim,:), 'MarkerEdgeColor','none', 'MarkerSize', 10)
                % Add confidence ellipse
                dat_pca_boot = bootstrp(5000, @nanmean, all_std_pca(ind,1:2));
                % Get a 99%-confidence ellipse around the bootstrapped data
                r_ellipse = Confocal_SubFcn.dataEllipse(dat_pca_boot, 0.99);
                plot(r_ellipse(:,1), r_ellipse(:,2), 'Color', [0.251 0.251 0.251], 'LineWidth', 1)
            end%iStim
            % Cosmetics
            subplot(1,2,1); hold on
            axis equal tight
            xlabel(['PC1 (', num2str(round(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.avg_explained(1), 2)),'pc)'])
            ylabel(['PC2 (', num2str(round(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.avg_explained(2), 2)),'pc)'])
            title('avg proje.')
            subplot(1,2,2); hold on
            axis equal tight
            xlabel(['PC1 (', num2str(round(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.std_explained(1), 2)),'pc)'])
            ylabel(['PC2 (', num2str(round(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', cnt_ani)]).PCA.std_explained(2), 2)),'pc)'])
            title('std proje.')
            export_fig([curr.dir.animal, '04_Data_Summary\pca_space'], '-pdf', '-painters')
            close(hFig)


            %% Dendrogram
            % --- avg projection
            D = pdist(all_avg_pca, "euclidean");
            tree = linkage(D,'average');
            hFig = figure();
            [~,~,dendro_label] = dendrogram(tree, 0);
            dendro_label_update = cell(1,length(dendro_label));
            % Get custom labels
            for iLabel = 1:length(dendro_label)
                dendro_label_update{iLabel} = SET.stim_list{all_label(dendro_label(iLabel))};
            end%iLabel
            set(gca, 'XTick', 1:length(dendro_label_update), 'XTickLabel', dendro_label_update, 'XTickLabelRotation', 90)
            ylabel('ED')
            title({'all_avg_pca'; ['c=', num2str(round(cophenet(tree,D),2))]}, 'interpreter', 'none')
            export_fig([curr.dir.animal, '04_Data_Summary\dendro_all_avg_pca'], '-pdf', '-painters')
            close(hFig)

            % --- avg projection
            D = pdist(all_std_pca, "euclidean");
            tree = linkage(D,'average');
            hFig = figure();
            [~,~,dendro_label] = dendrogram(tree, 0);
            dendro_label_update = cell(1,length(dendro_label));
            % Get custom labels
            for iLabel = 1:length(dendro_label)
                dendro_label_update{iLabel} = SET.stim_list{all_label(dendro_label(iLabel))};
            end%iLabel
            set(gca, 'XTick', 1:length(dendro_label_update), 'XTickLabel', dendro_label_update, 'XTickLabelRotation', 90)
            ylabel('ED')
            title({'all_std_pca'; ['c=', num2str(round(cophenet(tree,D),2))]}, 'interpreter', 'none')
            export_fig([curr.dir.animal, '04_Data_Summary\dendro_all_std_pca'], '-pdf', '-painters')
            close(hFig)

            cnt_ani = cnt_ani+1;
        end%if animal
    end%iAni

    % Pool each animals within/across distances
    w_avg = []; w_std = [];
    a_avg = []; a_std = [];
    for iAni = 1:length(fieldnames(PoolData.(SET.layers{iLayer})))
        w_avg = [w_avg; nanmean(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', iAni)]).PCA.dist_within_avg)];
        w_std = [w_std; nanmean(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', iAni)]).PCA.dist_within_std)];
        a_avg = [a_avg; nanmean(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', iAni)]).PCA.dist_across_avg)];
        a_std = [a_std; nanmean(PoolData.(SET.layers{iLayer}).(['Ani', sprintf('%02d', iAni)]).PCA.dist_across_std)];
    end%iAni

    hFig = figure('Units','normalized','Position',[0.25 0.25 0.5 0.5],'Color','w');

    % Within - avg projection
    subplot(2,3,1); hold on
    dist_data = w_avg;
    % Violin
    properties.MinVal =             min(dist_data);
    properties.MaxVal =             max(dist_data);
    properties.AvgType =            'mean';
    properties.MeanCol =            'k';
    properties.MeanSymbol =         'line';
    properties.MeanWidth =          2;
    properties.SeparateOutliers =   0;
    Confocal_SubFcn.violinplot_advanced(dist_data, 0, 0.5, properties)
    plot([0,0], [min(dist_data), max(dist_data)], 'k')
    clear properties
    % Swarm
    properties.MarkerType =         'o';
    properties.MarkerFaceColor =    [0.251 0.251 0.251];
    properties.MarkerSize =         5;
    Confocal_SubFcn.beeswarmplot_advanced(dist_data, 0.25, 0.5, properties)
    % Add mean and CI
    avg = mean(bootstrp(5000, @mean, dist_data));
    ci = bootci(5000, @mean, dist_data);
    plot([0.125 0.375], [avg avg], 'k')
    plot([0.25 0.25], ci, 'k')
    clear dist_data properties
    % Cosmetics
    xlim([-0.5 0.75])
    title('within, avg')
    ylim([min([w_avg(:); a_avg(:)])-min([w_avg(:); a_avg(:)])*0.15, max([w_avg(:); a_avg(:)])+max([w_avg(:); a_avg(:)])*0.15])

    % Across - avg projection
    subplot(2,3,2); hold on
    dist_data = a_avg;
    % Violin
    properties.MinVal =             min(dist_data);
    properties.MaxVal =             max(dist_data);
    properties.AvgType =            'mean';
    properties.MeanCol =            'k';
    properties.MeanSymbol =         'line';
    properties.MeanWidth =          2;
    properties.SeparateOutliers =   0;
    Confocal_SubFcn.violinplot_advanced(dist_data, 0, 0.5, properties)
    plot([0,0], [min(dist_data), max(dist_data)], 'k')
    clear properties
    % Swarm
    properties.MarkerType =         'o';
    properties.MarkerFaceColor =    [0.251 0.251 0.251];
    properties.MarkerSize =         5;
    Confocal_SubFcn.beeswarmplot_advanced(dist_data, 0.25, 0.5, properties)
    % Add mean and CI
    avg = mean(bootstrp(5000, @mean, dist_data));
    ci = bootci(5000, @mean, dist_data);
    plot([0.125 0.375], [avg avg], 'k')
    plot([0.25 0.25], ci, 'k')
    clear dist_data properties
    % Cosmetics
    xlim([-0.5 0.75])
    title('across, avg')
    ylim([min([w_avg(:); a_avg(:)])-min([w_avg(:); a_avg(:)])*0.15, max([w_avg(:); a_avg(:)])+max([w_avg(:); a_avg(:)])*0.15])

    subplot(2,3,3); hold on
    plot([1 2], [w_avg(:), a_avg(:)], 'k')
    ylim([min([w_avg(:); a_avg(:)])-min([w_avg(:); a_avg(:)])*0.15, max([w_avg(:); a_avg(:)])+max([w_avg(:); a_avg(:)])*0.15])
    xlim([0 3])
    % Do stats
    Stats.seed=1234;
    Stats.n = 5e7;
    Stats.nComp = 1;
    [...
        Stats.(SET.layers{iLayer}).within_across_avg.p,...
        Stats.(SET.layers{iLayer}).within_across_avg.s,...
        ~,...
        Stats.(SET.layers{iLayer}).within_across_avg.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', w_avg, a_avg, Stats.n, Stats.seed, Stats.nComp);
    title(['p=', num2str(Stats.(SET.layers{iLayer}).within_across_avg.p)])


    % Within - std projection
    subplot(2,3,4); hold on
    dist_data = w_std;
    % Violin
    properties.MinVal =             min(dist_data);
    properties.MaxVal =             max(dist_data);
    properties.stdType =            'mean';
    properties.MeanCol =            'k';
    properties.MeanSymbol =         'line';
    properties.MeanWidth =          2;
    properties.SeparateOutliers =   0;
    Confocal_SubFcn.violinplot_advanced(dist_data, 0, 0.5, properties)
    plot([0,0], [min(dist_data), max(dist_data)], 'k')
    clear properties
    % Swarm
    properties.MarkerType =         'o';
    properties.MarkerFaceColor =    [0.251 0.251 0.251];
    properties.MarkerSize =         5;
    Confocal_SubFcn.beeswarmplot_advanced(dist_data, 0.25, 0.5, properties)
    % Add mean and CI
    std = mean(bootstrp(5000, @mean, dist_data));
    ci = bootci(5000, @mean, dist_data);
    plot([0.125 0.375], [std std], 'k')
    plot([0.25 0.25], ci, 'k')
    clear dist_data properties
    % Cosmetics
    xlim([-0.5 0.75])
    title('within, std')
    ylim([min([w_std(:); a_std(:)])-min([w_std(:); a_std(:)])*0.15, max([w_std(:); a_std(:)])+max([w_std(:); a_std(:)])*0.15])

    % Across - std projection
    subplot(2,3,5); hold on
    dist_data = a_std;
    % Violin
    properties.MinVal =             min(dist_data);
    properties.MaxVal =             max(dist_data);
    properties.stdType =            'mean';
    properties.MeanCol =            'k';
    properties.MeanSymbol =         'line';
    properties.MeanWidth =          2;
    properties.SeparateOutliers =   0;
    Confocal_SubFcn.violinplot_advanced(dist_data, 0, 0.5, properties)
    plot([0,0], [min(dist_data), max(dist_data)], 'k')
    clear properties
    % Swarm
    properties.MarkerType =         'o';
    properties.MarkerFaceColor =    [0.251 0.251 0.251];
    properties.MarkerSize =         5;
    Confocal_SubFcn.beeswarmplot_advanced(dist_data, 0.25, 0.5, properties)
    % Add mean and CI
    std = mean(bootstrp(5000, @mean, dist_data));
    ci = bootci(5000, @mean, dist_data);
    plot([0.125 0.375], [std std], 'k')
    plot([0.25 0.25], ci, 'k')
    clear dist_data properties
    % Cosmetics
    xlim([-0.5 0.75])
    title('across, std')
    ylim([min([w_std(:); a_std(:)])-min([w_std(:); a_std(:)])*0.15, max([w_std(:); a_std(:)])+max([w_std(:); a_std(:)])*0.15])

    subplot(2,3,6); hold on
    plot([1 2], [w_std(:), a_std(:)], 'k')
    ylim([min([w_std(:); a_std(:)])-min([w_std(:); a_std(:)])*0.15, max([w_std(:); a_std(:)])+max([w_std(:); a_std(:)])*0.15])
    xlim([0 3])
    % Do stats
    Stats.seed=1234;
    Stats.n = 5e7;
    Stats.nComp = 1;
    [...
        Stats.(SET.layers{iLayer}).within_across_std.p,...
        Stats.(SET.layers{iLayer}).within_across_std.s,...
        ~,...
        Stats.(SET.layers{iLayer}).within_across_std.c] = Confocal_SubFcn.BootstrapHypothesisTesting('two-sample-pairs', w_std, a_std, Stats.n, Stats.seed, Stats.nComp);
    title(['p=', num2str(Stats.(SET.layers{iLayer}).within_across_std.p)])    
    % Save
    export_fig(['ConfocalAnalysis_StaticResponse\',SET.layers{iLayer},'_population_response'], '-pdf', '-painters')
    close(hFig)
    save('ConfocalAnalysis_StaticResponse\StaticResponse_data.mat', 'PoolData', 'Stats')

end%iLayer



























