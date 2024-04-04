%% EAG_Analysis_Population
% bla
%
% Version: 04-April-24 (R2023a)

% Clean up
clc; clear all; close all;
% Add paths:
% --- Current paths
addpath(genpath(pwd))
% --- Toolbox for exporting publication quality figures
%     https://github.com/altmany/export_fig
addpath(genpath('...\GitHub\export_fig'))
mkdir('FIGs')

%% Settings

% Give list of stimuli (We have this in the meta data, but good to double
% check)
SET.StimList = {...
    'MOL';...
    'COL';...
    'ZHAE';...
    'ZHAECOL'};

% Cosmetics
SET.Colors = [...
    200,200,200;...MOL
    217,095,002;...COL
    102,166,030;...ZHAE
    231,041,138;...ZHAECOL  
    ]/255;

% Phases
SET.phases = {'gregarious'; 'solitarious'};

% Info on data collection
SET.freq = 25e3;
SET.crop_rep = 2;
SET.stim_dur = 2;
SET.stim_window = [2 4]*SET.freq;

% Set whether to pool antennae
SET.pool_ants = 0;

% Paths
% --- where can we find the data
SET.path2data = '...\EAG\';
% --- where can we find the annotation
SET.path2annotation = 'ANNOTATION.csv';

%% Pool each animal's response

% Get annotation
ANNOTATION = readtable(SET.path2annotation);

% Get list of animals
[~,ind_ani] = unique(ANNOTATION(:,1:3), "rows");

% Preallocation
PoolData.gregarious.info = {};
PoolData.solitarious.info = {};
for iStim = 1:length(SET.StimList)
    for iPhase = 1:length(SET.phases)
        PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).TC_original = [];
        PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val_original = [];
        PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).TC = [];
        PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val = [];
    end%iPhase
end%iStim

% Iterate over all animals
for iAni = 1:length(ind_ani)
    % Get rows in annotation
    curr_rows = find(...
        strcmp(ANNOTATION.phase, ANNOTATION.phase{ind_ani(iAni)}) & ...
        ANNOTATION.date == ANNOTATION.date(ind_ani(iAni)) & ...
        ANNOTATION.animal == ANNOTATION.animal(ind_ani(iAni)));

    % Continue if there is something valid
    if any([ANNOTATION.rep1(curr_rows); ANNOTATION.rep2(curr_rows); ANNOTATION.rep3(curr_rows)])
        % Get antennae
        curr_ant = unique(ANNOTATION.ant(curr_rows));
        % Preallocation
        for iStim = 1:length(SET.StimList)
            pool.(SET.StimList{iStim}) = [];
            pool_original.(SET.StimList{iStim}) = [];
        end%iStim
        % Iterate over both antennae
        for iAnt = 1:length(curr_ant)
            curr_trial = find(...
                strcmp(ANNOTATION.phase, ANNOTATION.phase{ind_ani(iAni)}) & ...
                ANNOTATION.date == ANNOTATION.date(ind_ani(iAni)) & ...
                ANNOTATION.animal == ANNOTATION.animal(ind_ani(iAni)) & ...
                ANNOTATION.ant == curr_ant(iAnt));
            % Load each stimulus
            for iStim = 1:length(curr_trial)
                % Construct path to file
                filename = [...
                    num2str(ANNOTATION.date(curr_trial(iStim))),'_Animal',sprintf('%02d',ANNOTATION.animal(curr_trial(iStim))),...
                    '_ant0',num2str(ANNOTATION.ant(curr_trial(iStim))),'_',ANNOTATION.ant_side{curr_trial(iStim)},...
                    '_',ANNOTATION.stim{curr_trial(iStim)},...
                    '_',ANNOTATION.phase{curr_trial(iStim)}(1:4)];
                path2file = [...
                    SET.path2data,...
                    ANNOTATION.phase{curr_trial(iStim)},'\',...
                    num2str(ANNOTATION.date(curr_trial(iStim))),'_Animal',sprintf('%02d',ANNOTATION.animal(curr_trial(iStim))),'\',...
                    'ant0',num2str(ANNOTATION.ant(curr_trial(iStim))),'_',ANNOTATION.ant_side{curr_trial(iStim)},'\',...
                    filename,'.mat'];
                % Load
                curr_data.(ANNOTATION.stim{curr_trial(iStim)}) = load(path2file);
                % Check if there was an invalid rep. If so, set it to NaN
                for iRep = 1:3
                    if ANNOTATION.(['rep',num2str(iRep)])(curr_trial(iStim)) == 0
                        curr_data.(ANNOTATION.stim{curr_trial(iStim)}).data.(['wave_data_original_rep',num2str(iRep)]) = nan(size(curr_data.(ANNOTATION.stim{curr_trial(iStim)}).data.wave_data_original_rep1));
                        curr_data.(ANNOTATION.stim{curr_trial(iStim)}).data.(['wave_data_rep',num2str(iRep)]) = nan(size(curr_data.(ANNOTATION.stim{curr_trial(iStim)}).data.wave_data_rep1));
                    end%if
                end%iRep
            end%iStim

            % Now, since we managed in an overcomplicated way to load some
            % data, normalize each repetition with the respective response
            % to ZHAE. Afterward, pool all 3 reps to one response for the
            % current antenna and normalize again (due to slight
            % differences between Zhae repetitions)
            for iRep = 1:3
                min_val_original = abs(min( curr_data.ZHAE.data.(['wave_data_original_rep',num2str(iRep)])(SET.stim_window(1):SET.stim_window(2)) ));
                min_val = abs(min( curr_data.ZHAE.data.(['wave_data_rep',num2str(iRep)])(SET.stim_window(1):SET.stim_window(2)) ));
                for iStim = 1:length(SET.StimList)
                    curr_data.(SET.StimList{iStim}).data.(['wave_data_original_rep',num2str(iRep)]) = ...
                        curr_data.(SET.StimList{iStim}).data.(['wave_data_original_rep',num2str(iRep)]) / min_val_original;
                    curr_data.(SET.StimList{iStim}).data.(['wave_data_rep',num2str(iRep)]) = ...
                        curr_data.(SET.StimList{iStim}).data.(['wave_data_rep',num2str(iRep)]) / min_val;
                end%iStim
            end%iRep
            for iStim = 1:length(SET.StimList)
                temp_original.(SET.StimList{iStim}) = [];
                temp.(SET.StimList{iStim}) = [];
                for iRep = 1:3
                    temp_original.(SET.StimList{iStim}) = [temp_original.(SET.StimList{iStim}) ; ...
                        curr_data.(SET.StimList{iStim}).data.(['wave_data_original_rep',num2str(iRep)])(:)'];
                    temp.(SET.StimList{iStim}) = [temp.(SET.StimList{iStim}) ; ...
                        curr_data.(SET.StimList{iStim}).data.(['wave_data_rep',num2str(iRep)])(:)'];
                end%iRep
                temp_original.(SET.StimList{iStim}) = nanmean(temp_original.(SET.StimList{iStim}), 1);
                temp.(SET.StimList{iStim}) = nanmean(temp.(SET.StimList{iStim}), 1);
            end%iStim
            min_val_original = abs(min( temp_original.ZHAE(SET.stim_window(1):SET.stim_window(2)) ));
            min_val = abs(min( temp.ZHAE(SET.stim_window(1):SET.stim_window(2)) ));
            for iStim = 1:length(SET.StimList)
                pool_original.(SET.StimList{iStim}) = [pool_original.(SET.StimList{iStim}); temp_original.(SET.StimList{iStim}) / min_val_original];
                pool.(SET.StimList{iStim}) = [pool.(SET.StimList{iStim}); temp.(SET.StimList{iStim}) / min_val];
            end% iStim
        end%iAnt

        % Pool both antennae and normalize again to ensure that the minimum
        % response to Zhae is -1
        if SET.pool_ants
        for iStim = 1:length(SET.StimList)
            pool_original.(SET.StimList{iStim}) = nanmean(pool_original.(SET.StimList{iStim}), 1);
            pool.(SET.StimList{iStim}) = nanmean(pool.(SET.StimList{iStim}), 1);
        end%iStim
        min_val_original = abs(min( pool_original.ZHAE(SET.stim_window(1):SET.stim_window(2)) ));
        min_val = abs(min( pool.ZHAE(SET.stim_window(1):SET.stim_window(2)) ));
        for iStim = 1:length(SET.StimList)
            pool_original.(SET.StimList{iStim}) = pool_original.(SET.StimList{iStim}) / min_val_original;
            pool.(SET.StimList{iStim}) = pool.(SET.StimList{iStim}) / min_val;
        end%iStim
        end%if

        % Pool across all animals
        PoolData.(ANNOTATION.phase{ind_ani(iAni)}).info = [...
            PoolData.(ANNOTATION.phase{ind_ani(iAni)}).info;...
            [num2str(ANNOTATION.date(ind_ani(iAni))), '_Animal', sprintf('%02d',ANNOTATION.animal(ind_ani(iAni)))]];
        for iStim = 1:length(SET.StimList)

            % Whole time course
            PoolData.(ANNOTATION.phase{ind_ani(iAni)}).(SET.StimList{iStim}).TC_original = [...
                PoolData.(ANNOTATION.phase{ind_ani(iAni)}).(SET.StimList{iStim}).TC_original;...
                pool_original.(SET.StimList{iStim})];
            % Only the minimum value
            PoolData.(ANNOTATION.phase{ind_ani(iAni)}).(SET.StimList{iStim}).val_original = [...
                PoolData.(ANNOTATION.phase{ind_ani(iAni)}).(SET.StimList{iStim}).val_original;...
                min(pool_original.(SET.StimList{iStim})(:, SET.stim_window(1):SET.stim_window(2)), [], 2)];

            % Whole time course
            PoolData.(ANNOTATION.phase{ind_ani(iAni)}).(SET.StimList{iStim}).TC = [...
                PoolData.(ANNOTATION.phase{ind_ani(iAni)}).(SET.StimList{iStim}).TC;...
                pool.(SET.StimList{iStim})];
            % Only the minimum value
            PoolData.(ANNOTATION.phase{ind_ani(iAni)}).(SET.StimList{iStim}).val = [...
                PoolData.(ANNOTATION.phase{ind_ani(iAni)}).(SET.StimList{iStim}).val;...
                min(pool.(SET.StimList{iStim})(:, SET.stim_window(1):SET.stim_window(2)), [], 2)];

        end%iStim
    end%if any valid
end%iAni

%% Plot results

hFig = figure('Units','normalized','Position',[0 0 1 1],'Color','w');
time_vec = linspace(-SET.crop_rep, size(PoolData.gregarious.MOL.TC_original,2)/SET.freq, size(PoolData.gregarious.MOL.TC_original,2));
% Plot average traces to each stimulus in its own subplot.
% Add mol to each
% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Get data
    avg_mol = nanmean(PoolData.(SET.phases{iPhase}).MOL.TC_original,1);
    ci_mol = bootci(5000, {@nanmean, PoolData.(SET.phases{iPhase}).MOL.TC_original});
    avg_col = nanmean(PoolData.(SET.phases{iPhase}).COL.TC_original,1);
    ci_col = bootci(5000, {@nanmean, PoolData.(SET.phases{iPhase}).COL.TC_original});
    avg_zhae = nanmean(PoolData.(SET.phases{iPhase}).ZHAE.TC_original,1);
    ci_zhae = bootci(5000, {@nanmean, PoolData.(SET.phases{iPhase}).ZHAE.TC_original});
    avg_mix = nanmean(PoolData.(SET.phases{iPhase}).ZHAECOL.TC_original,1);
    ci_mix = bootci(5000, {@nanmean, PoolData.(SET.phases{iPhase}).ZHAECOL.TC_original});

    subplot(2,4,1+((iPhase-1)*4)); hold on
    % -- MOL
    plot(time_vec, avg_mol, "Color", SET.Colors(1,:), 'LineWidth', 2)
    plot(time_vec, ci_mol(1,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    plot(time_vec, ci_mol(2,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    % --- COL
    plot(time_vec, avg_col, "Color", SET.Colors(2,:), 'LineWidth', 2)
    plot(time_vec, ci_col(1,:), "Color", SET.Colors(2,:), 'LineWidth', 1)
    plot(time_vec, ci_col(2,:), "Color", SET.Colors(2,:), 'LineWidth', 1)
    % Cosmetics
    axis square
    xlim([time_vec(1), time_vec(end)])
    ylim([-1.4 0.2])
    ylabel('relative response magnitude')
    xlabel('time (s)')


    subplot(2,4,2+((iPhase-1)*4)); hold on
    % -- MOL
    plot(time_vec, avg_mol, "Color", SET.Colors(1,:), 'LineWidth', 2)
    plot(time_vec, ci_mol(1,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    plot(time_vec, ci_mol(2,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    % --- ZHAE
    plot(time_vec, avg_zhae, "Color", SET.Colors(3,:), 'LineWidth', 2)
    plot(time_vec, ci_zhae(1,:), "Color", SET.Colors(3,:), 'LineWidth', 1)
    plot(time_vec, ci_zhae(2,:), "Color", SET.Colors(3,:), 'LineWidth', 1)
    % Cosmetics
    axis square
    xlim([time_vec(1), time_vec(end)])
    ylim([-1.4 0.2])
    ylabel('relative response magnitude')
    xlabel('time (s)')


    subplot(2,4,3+((iPhase-1)*4)); hold on
    % -- MOL
    plot(time_vec, avg_mol, "Color", SET.Colors(1,:), 'LineWidth', 2)
    plot(time_vec, ci_mol(1,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    plot(time_vec, ci_mol(2,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    % --- ZHAECOL
    plot(time_vec, avg_mix, "Color", SET.Colors(4,:), 'LineWidth', 2)
    plot(time_vec, ci_mix(1,:), "Color", SET.Colors(4,:), 'LineWidth', 1)
    plot(time_vec, ci_mix(2,:), "Color", SET.Colors(4,:), 'LineWidth', 1)
    % Cosmetics
    axis square
    xlim([time_vec(1), time_vec(end)])
    ylim([-1.4 0.2])
    ylabel('relative response magnitude')
    xlabel('time (s)')


    subplot(2,4,4+((iPhase-1)*4)); hold on
    for iStim = 1:length(SET.StimList)
        % Violin
        vioData = PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val_original;
        clear properties
        properties.MinVal = nanmin(vioData);
        properties.MaxVal = nanmax(vioData);
        properties.AvgType = 'mean';
        properties.MeanCol = 'k';
        properties.MeanSymbol = 'line';
        properties.MeanWidth = 2;
        properties.SeparateOutliers = 0;
        EAG_SubFcn.violinplot_advanced(vioData, iStim, 0.5, properties)
        % Swarm
        beeData = PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val_original;
        properties.MarkerType = 'o';
        properties.MarkerFaceColor = SET.Colors(iStim,:);
        properties.MarkerSize = 5;
        EAG_SubFcn.beeswarmplot_advanced(beeData, iStim, 0.5, properties)
        % Avg and CI
        avg = nanmean(PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val_original);
        ci = bootci(5000, {@nanmean, PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val_original});
        plot([1 1]*iStim, ci,'k')
        plot([-0.25 0.25]+iStim, [avg avg],'k')
        % Statistics
        [...
        STATS.(SET.phases{iPhase}).(SET.StimList{iStim}).p_original,...
        STATS.(SET.phases{iPhase}).(SET.StimList{iStim}).s_original,] = EAG_SubFcn.BootstrapHypothesisTesting(...
        'one-sample',...
        PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val_original,...
        -1, ...
        1e7, 1234, 1);
    end%iStim
    % Cosmetics
    axis square
    set(gca, 'XTick', 1:length(SET.StimList), 'XTickLabel', SET.StimList)
    ylabel('relative response magnitude')
    ylim([-2 0.2])
    

end %iPhase

% Compare phases
for iStim = 1:length(SET.StimList)
    [...
        STATS.greg_vs_soli.(SET.StimList{iStim}).p_original,...
        STATS.greg_vs_soli.(SET.StimList{iStim}).s_original,] = EAG_SubFcn.BootstrapHypothesisTesting(...
        'two-sample',...
        PoolData.gregarious.(SET.StimList{iStim}).val_original(~isnan(PoolData.gregarious.(SET.StimList{iStim}).val_original)),...
        PoolData.solitarious.(SET.StimList{iStim}).val_original(~isnan(PoolData.solitarious.(SET.StimList{iStim}).val_original)),...
        1e7, 1234, 1);
end%iStim

% Save
export_fig('FIGs\EAG_results_original', '-pdf', '-painters')
close(hFig)

%% Plot results

hFig = figure('Units','normalized','Position',[0 0 1 1],'Color','w');
time_vec = linspace(-SET.crop_rep, size(PoolData.gregarious.MOL.TC,2)/SET.freq, size(PoolData.gregarious.MOL.TC,2));
% Plot average traces to each stimulus in its own subplot.
% Add mol to each
% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Get data
    avg_mol = nanmean(PoolData.(SET.phases{iPhase}).MOL.TC,1);
    ci_mol = bootci(5000, {@nanmean, PoolData.(SET.phases{iPhase}).MOL.TC});
    avg_col = nanmean(PoolData.(SET.phases{iPhase}).COL.TC,1);
    ci_col = bootci(5000, {@nanmean, PoolData.(SET.phases{iPhase}).COL.TC});
    avg_zhae = nanmean(PoolData.(SET.phases{iPhase}).ZHAE.TC,1);
    ci_zhae = bootci(5000, {@nanmean, PoolData.(SET.phases{iPhase}).ZHAE.TC});
    avg_mix = nanmean(PoolData.(SET.phases{iPhase}).ZHAECOL.TC,1);
    ci_mix = bootci(5000, {@nanmean, PoolData.(SET.phases{iPhase}).ZHAECOL.TC});

    subplot(2,4,1+((iPhase-1)*4)); hold on
    % -- MOL
    plot(time_vec, avg_mol, "Color", SET.Colors(1,:), 'LineWidth', 2)
    plot(time_vec, ci_mol(1,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    plot(time_vec, ci_mol(2,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    % --- COL
    plot(time_vec, avg_col, "Color", SET.Colors(2,:), 'LineWidth', 2)
    plot(time_vec, ci_col(1,:), "Color", SET.Colors(2,:), 'LineWidth', 1)
    plot(time_vec, ci_col(2,:), "Color", SET.Colors(2,:), 'LineWidth', 1)
    % Cosmetics
    axis square
    xlim([time_vec(1), time_vec(end)])
    ylim([-1.2 0.5])
    ylabel('relative response magnitude')
    xlabel('time (s)')


    subplot(2,4,2+((iPhase-1)*4)); hold on
    % -- MOL
    plot(time_vec, avg_mol, "Color", SET.Colors(1,:), 'LineWidth', 2)
    plot(time_vec, ci_mol(1,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    plot(time_vec, ci_mol(2,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    % --- ZHAE
    plot(time_vec, avg_zhae, "Color", SET.Colors(3,:), 'LineWidth', 2)
    plot(time_vec, ci_zhae(1,:), "Color", SET.Colors(3,:), 'LineWidth', 1)
    plot(time_vec, ci_zhae(2,:), "Color", SET.Colors(3,:), 'LineWidth', 1)
    % Cosmetics
    axis square
    xlim([time_vec(1), time_vec(end)])
    ylim([-1.2 0.5])
    ylabel('relative response magnitude')
    xlabel('time (s)')


    subplot(2,4,3+((iPhase-1)*4)); hold on
    % -- MOL
    plot(time_vec, avg_mol, "Color", SET.Colors(1,:), 'LineWidth', 2)
    plot(time_vec, ci_mol(1,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    plot(time_vec, ci_mol(2,:), "Color", SET.Colors(1,:), 'LineWidth', 1)
    % --- ZHAECOL
    plot(time_vec, avg_mix, "Color", SET.Colors(4,:), 'LineWidth', 2)
    plot(time_vec, ci_mix(1,:), "Color", SET.Colors(4,:), 'LineWidth', 1)
    plot(time_vec, ci_mix(2,:), "Color", SET.Colors(4,:), 'LineWidth', 1)
    % Cosmetics
    axis square
    xlim([time_vec(1), time_vec(end)])
    ylim([-1.2 0.5])
    ylabel('relative response magnitude')
    xlabel('time (s)')


    subplot(2,4,4+((iPhase-1)*4)); hold on
    for iStim = 1:length(SET.StimList)
        % Violin
        vioData = PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val;
        clear properties
        properties.MinVal = nanmin(vioData);
        properties.MaxVal = nanmax(vioData);
        properties.AvgType = 'mean';
        properties.MeanCol = 'k';
        properties.MeanSymbol = 'line';
        properties.MeanWidth = 2;
        properties.SeparateOutliers = 0;
        EAG_SubFcn.violinplot_advanced(vioData, iStim, 0.5, properties)
        % Swarm
        beeData = PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val;
        properties.MarkerType = 'o';
        properties.MarkerFaceColor = SET.Colors(iStim,:);
        properties.MarkerSize = 5;
        EAG_SubFcn.beeswarmplot_advanced(beeData, iStim, 0.5, properties)
        % Avg and CI
        avg = nanmean(PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val);
        ci = bootci(5000, {@nanmean, PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val});
        plot([1 1]*iStim, ci,'k')
        plot([-0.25 0.25]+iStim, [avg avg],'k')
        % Statistics
        [...
        STATS.(SET.phases{iPhase}).(SET.StimList{iStim}).p,...
        STATS.(SET.phases{iPhase}).(SET.StimList{iStim}).s,] = EAG_SubFcn.BootstrapHypothesisTesting(...
        'one-sample',...
        PoolData.(SET.phases{iPhase}).(SET.StimList{iStim}).val,...
        -1, ...
        1e7, 1234, 1);
    end%iStim
    % Cosmetics
    axis square
    set(gca, 'XTick', 1:length(SET.StimList), 'XTickLabel', SET.StimList)
    ylabel('relative response magnitude')
    ylim([-1.5, 0.125])      

end %iPhase

% Compare phases
for iStim = 1:length(SET.StimList)
    [...
        STATS.greg_vs_soli.(SET.StimList{iStim}).p,...
        STATS.greg_vs_soli.(SET.StimList{iStim}).s,] = EAG_SubFcn.BootstrapHypothesisTesting(...
        'two-sample',...
        PoolData.gregarious.(SET.StimList{iStim}).val(~isnan(PoolData.gregarious.(SET.StimList{iStim}).val)),...
        PoolData.solitarious.(SET.StimList{iStim}).val(~isnan(PoolData.solitarious.(SET.StimList{iStim}).val)),...
        1e7, 1234, 1);
end%iStim

% Save
export_fig('FIGs\EAG_results', '-pdf', '-painters')
close(hFig)

%% Save data

save('PoolData.mat', 'SET', 'PoolData', 'STATS')











