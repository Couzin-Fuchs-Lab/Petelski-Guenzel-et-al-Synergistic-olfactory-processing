%% FourCageDecision_Analysis
% Extract, analyze, and pool data from the 4cage experiments.
% For these experiments, a single animal (gregarious or solitarious) was
% places in a circular arena (d=90cm) with four stimulus cages, placed on
% opposing sides. The four cages were (i) an empty control Ctr, (ii) a cage
% with fresh blackberry leaves Lvs, (iii) a cage with eitgh gregarious 5th
% instar hoppers Lct, and (iv) a cage with fresh blackberry leaves and
% eight gregarious 5th instar hoppers LvsLct. These cages were either
% transparent with holes (access to both visual and olfactory information),
% opaque (blocked visual inforamtion), or sealed (blocked olfactory
% inforamtion).
%
% Version:
% 20-April-2023 (R2023a) Yannick Günzel

% Prepare
clc; clear all; close all
mkdir('Output')
warning('off')

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig) 
addpath(genpath('...\GitHub\export_fig'))

%% Settings

% ***** Info ON DATA COLLECTION *****
% Set frame rate of recordings
SET.FrameRate = 25; %[fps]
% Set arena diameter
SET.dArena = 90; %[cm]

% ***** DATA ANALYSIS *****
% Set smoothing parameters
SET.g_filter = 'gauss';
SET.g_hw = 2; %(s)
SET.g_hw = SET.g_hw*SET.FrameRate;
SET.g_sigma = SET.g_hw/3;
% Set additional distance to patch edge as tolerance
SET.ToleranceDistance = 4; %[cm, 4cm=1BL]
% Set bins in which data will be divided for the heatmaps. Resulting number
% of bins will be SET.HeatMapGrid^2.
SET.HeatMapGrid = 1000;
% Interaction range for density filter
SET.InteractionRange_SD = 7; %(cm)
% Std for imgausfilt()
SET.dArena = 90; %[cm]
SET.SmoothValue = (SET.HeatMapGrid/SET.dArena)*SET.InteractionRange_SD;
% Bootsrapping
SET.BootSamples = 5000;
% The minimum inter-bout-interval and bout duration
SET.minInterBoutInterval = 5; %(s, animal-centric)
SET.minBoutDuration = 5;      %(s)
SET.minInterBoutInterval = SET.minInterBoutInterval*SET.FrameRate;
SET.minBoutDuration = SET.minBoutDuration*SET.FrameRate;
% Size of pizza slices for heatmaps
SET.HeatmapSliceSize = 45; %(deg, halfwidth)
% Set vector for binning angular location
SET.AngularLocBinning = 0:10:360;


% ***** PATHS *****
SET.BasePath = '...\Data\Behaviour\FourCageDecision\';
SET.Phase = {'gregarious', 'solitarious'};
SET.Conditions.cages = {...
    'VpOpFp',... visual cues (see-through),   feeding
    'VpOpFm',... visual cues (see-through),   not feeding
    'VmOpFp',... no visual cues(opaque),      feeding
    'VmOpFm',... no visual cues(opaque),      not feeding
    'VpOmFp',... no olfactory cues (glass),   feeding
    'VpOmFm',... no olfactory cues (glass),   not feeding
    };

SET.Cages = {...
    '_Ctr_',...       empty control
    '_Leaves_',...    blackberry leaves
    '_Lo_'...         group of locusts
    '_LoLeaves_',...  group of locusts + blackberry leaves
    };

% ***** COSMETICS *****
% Cosmetics
SET.Colors = [...
    191,191,191;...Ctr
    102,166,030;...Leaves
    217,095,002;...Lo
    221,041,138;...LoLeaves
    ]/255;
SET.phaseColors.gregarious  = [191 0 0]/255;
SET.phaseColors.solitarious = [0 0 191]/255;
SET.caxis_pizzamap = [0, 0.2];
SET.Colors_Cotton = {' ','LoLeavesZ','LoLeavesCo','LoLeavesZCo'};
SET.HeatmapColors.gregarious = FourCageDecision_SubFcn.ColMapPlasma(1000);
SET.HeatmapColors.solitarious = [...
    SET.HeatmapColors.gregarious(:,2),...
    SET.HeatmapColors.gregarious(:,1),...
    SET.HeatmapColors.gregarious(:,3)];
%%

% Iterate over phases
for iPhase = 1:length(SET.Phase)

    %----------------------------------------------------------------------
    %                                            Experiment with four cages
    %----------------------------------------------------------------------

    % Preallocation
    Data4Dendro.(SET.Phase{iPhase}) = nan(length(SET.Conditions.cages), length(SET.Cages));
    Data4Dendro_std.(SET.Phase{iPhase}) = nan(length(SET.Conditions.cages), length(SET.Cages));

    % Iterate over conditions
    for iCond = 1:length(SET.Conditions.cages)
        % Get list of files
        currFileList = dir([SET.BasePath, SET.Phase{iPhase}, '\', SET.Conditions.cages{iCond}]);
        % Check whether there is something in there
        if sum(~[currFileList(:).isdir])>2
            % Preallocation
            for iCage = 1:length(SET.Cages)
                PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).Pizza.(SET.Cages{iCage}(2:end-1)) = [];
            end
            % Iterate over file.Note that the following code assumes that
            % both annotation and trajectories are in the same folder
            % Counter
            cnt = 1;
            % Prespecification
            PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D = [];
            PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D_ID = [];
            for iFile = 1:size(currFileList,1)
                % Check whether it is a file containing trajectory data
                idx = strfind(currFileList(iFile).name, '_tracked');
                if ~isempty(idx)

                    % ***** PREPARE DATA *****
                    % Get base name
                    currFile.Basename = currFileList(iFile).name(1:idx-1);
                    % Get data
                    currFile = FourCageDecision_SubFcn.GetRawData(currFile, SET, [iPhase, iCond, 1]);
                    % Get order of conditions
                    currFile = FourCageDecision_SubFcn.GetConditionOrder(currFile, SET.Cages);
                    % Center and normalize data
                    currFile = FourCageDecision_SubFcn.CenterNormData(currFile);

                    % ***** EXTRACT BASIC METRICS *****
                    % Speed profile
                    currFile = FourCageDecision_SubFcn.SpeedProfile(currFile);
                    % Heatmaps with quaters arranged like pizza slices
                    currFile = FourCageDecision_SubFcn.PizzaHeatMap(currFile, SET);
                    % Distance to cages
                    currFile = FourCageDecision_SubFcn.Dist2Cage(currFile, SET.Cages);
                    % Determine whether animal is at the cage
                    currFile = FourCageDecision_SubFcn.AtCage(currFile, SET, SET.Cages);

                    % ***** STORE DATA FOR EACH TRIAL *****
                    DATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}){cnt,1} = currFile;
                    
                    % ***** STORE DATA FOR EACH CONDITION *****
                    for iCage = 1:length(SET.Cages)
                        % --- median moving speed
                        PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('avgSpeedMoving')(cnt,1) = median(currFile.SpeedProfile(currFile.SpeedProfile>quantile(currFile.SpeedProfile,0.25)));
                        % --- Pizza
                        PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).Pizza.(SET.Cages{iCage}(2:end-1)) = [...
                            PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).Pizza.(SET.Cages{iCage}(2:end-1));...
                            currFile.Pizza.(SET.Cages{iCage}(2:end-1))];
                        % --- Dist2Cage
                        PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('Dist2Cage')(cnt).(SET.Cages{iCage}(2:end-1)) = nanmean(currFile.Dist2Cage(iCage));
                        % --- PizzaTime
                        PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('PizzaTime')(cnt).(SET.Cages{iCage}(2:end-1)) = numel(currFile.Pizza.(SET.Cages{iCage}(2:end-1)))/numel([currFile.Pizza.Ctr; currFile.Pizza.Leaves; currFile.Pizza.Lo; currFile.Pizza.LoLeaves]);
                        % --- PropTimeAtCage
                        PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('PropTimeAtCage')(cnt).(SET.Cages{iCage}(2:end-1)) = sum(currFile.AtCage(:,iCage))/sum(sum(currFile.AtCage));
                        % --- NumberOfVisits
                        PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('NumberOfVisits')(cnt).(SET.Cages{iCage}(2:end-1)) = size(regionprops(currFile.AtCage(:,iCage)),1);
                        % --- Average density in slice, PizzaDensity "Topping density"
                        PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('PizzaDensity')(cnt).(SET.Cages{iCage}(2:end-1)) = nanmean(currFile.PizzaDensity.(SET.Cages{iCage}(2:end-1))(:));
                        % --- Average density at cage
                        PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('CageDensity')(cnt).(SET.Cages{iCage}(2:end-1)) = currFile.CageDensity.(SET.Cages{iCage}(2:end-1));
                        % --- 2D pizza density
                        PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D = cat(3,PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D, currFile.PizzaDensity.(SET.Cages{iCage}(2:end-1)));
                        PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D_ID = cat(3,PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D_ID, ~isnan(currFile.PizzaDensity.(SET.Cages{iCage}(2:end-1)))*iCage );
                    end%iCage
                    cnt = cnt+1;
                    clear currFile
                end% if trajectory data
            end%iFile

            % Get the average map across trials
            PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D = nanmean(PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D, 3);
            PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D_ID = nanmax(PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D_ID, [],3);

            % Crop map
            % Get indices of pixels
            clear idx
            [idx(:,1), idx(:,2)] = find(PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D>=-inf);
            % Get outlying pixels
            helper = idx-[SET.HeatMapGrid/2, SET.HeatMapGrid/2];
            dist = sqrt(sum(helper'.*helper'))';
            idx_dist = find(dist>SET.HeatMapGrid/2);
            % Set all outlying pixels to zero
            PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D(idx_dist) = nan(1,length(idx_dist));


            % ***** EXPORT RAW RESULTS *****
            T = table(PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('avgSpeedMoving'), 'VariableNames', {'avgSpeedMoving'});
            filename = ['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'_avgSpeedMoving.csv'];
            writetable(T,filename)
            clear T filename
            T = struct2table(PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('PropTimeAtCage'));
            filename = ['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'_PropTimeAtCage.csv'];
            writetable(T,filename)
            clear T filename
            T = struct2table(PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('PizzaTime'));
            filename = ['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'PizzaTime.csv'];
            writetable(T,filename)
            clear T filename
            T = struct2table(PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('PizzaDensity'));
            filename = ['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'PizzaDensity.csv'];
            writetable(T,filename)
            clear T filename
            T = struct2table(PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).('CageDensity'));
            filename = ['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'CageDensity.csv'];
            writetable(T,filename)
            clear T filename


            % ***** CREATE Pizza PLOTS *****
            currData = PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).Pizza;
            currData_ALL = [currData.Ctr; currData.Leaves; currData.Lo; currData.LoLeaves];
            % Bin and depict data
            bins = hist3(currData_ALL, 'ctrs',{linspace(-1,1,SET.HeatMapGrid) linspace(-1,1,SET.HeatMapGrid)});
            bins = (bins/sum(sum(bins)))*100;
            hFig = figure('Color', 'w');
            imagesc(bins)
            hold on
            plot([0, SET.HeatMapGrid],[0, SET.HeatMapGrid],'w', 'LineWidth', 2)
            plot([0, SET.HeatMapGrid],[SET.HeatMapGrid, 0],'w', 'LineWidth', 2)
            caxis([0, quantile(bins(:), 1-0.05)])
            colormap(SET.HeatmapColors.(SET.Phase{iPhase}))
            axis equal; axis off
            title(SET.Conditions.cages{iCond})
            set(gca, 'view', [-90 90])
            title(SET.Conditions.cages{iCond})
            export_fig(['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'_PizzaMap'], '-pdf', '-painters')
            close(hFig)
            clear currData* bins hFig

            % ***** CREATE PizzaDesnity PLOTS *****
            hFig = figure('Color', 'w');
            imagesc(PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity_2D)
            caxis(SET.caxis_pizzamap)
            cols = SET.HeatmapColors.(SET.Phase{iPhase});
            cols = [1 1 1; cols];
            colormap(cols)
            hold on
            plot([0, SET.HeatMapGrid],[0, SET.HeatMapGrid],'w', 'LineWidth', 2)
            plot([0, SET.HeatMapGrid],[SET.HeatMapGrid, 0],'w', 'LineWidth', 2)
            FourCageDecision_SubFcn.DrawCircle(SET.HeatMapGrid/2,SET.HeatMapGrid/2,SET.HeatMapGrid/2, 'w', 2)
            axis equal
            xlim([0 SET.HeatMapGrid]); xticks([])
            ylim([0 SET.HeatMapGrid]); yticks([])
            text(SET.HeatMapGrid*0.75, SET.HeatMapGrid/2,     SET.Cages{1},"HorizontalAlignment","center")
            text(SET.HeatMapGrid*0.5,  SET.HeatMapGrid*0.25,  SET.Cages{2},"HorizontalAlignment","center")
            text(SET.HeatMapGrid*0.5,  SET.HeatMapGrid*0.75,  SET.Cages{3},"HorizontalAlignment","center")
            text(SET.HeatMapGrid*0.25, SET.HeatMapGrid*0.5,   SET.Cages{4},"HorizontalAlignment","center")
            title(SET.Conditions.cages{iCond})
            export_fig(['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'_PizzaDensityMap'], '-pdf', '-painters')
            close(hFig)

            % ***** CREATE PropTimeAtPatch PLOTS *****
            % Collect data
            currData = PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PropTimeAtCage;
            AllData = [[currData.Ctr]', [currData.Leaves]', [currData.Lo]', [currData.LoLeaves]'];
            % Bootstrapping
            bootMean = mean(bootstrp(5000,@nanmean,AllData));
            bootCI = bootci(5000,@nanmean,AllData);
            % Figure
            hFig = figure('Color', 'w');
            hold on
            % Line at zero
            plot([0.75,length(bootMean)+0.75+0.25],[0 0],'k')
            % Iterate over ythe 4 cages
            for iStim = 1:length(bootMean)
                % Plot mean+-CI as bar plot
                rectangle('Position', [iStim, 0, 0.75, bootMean(iStim)], 'EdgeColor', 'none', 'FaceColor', [0 0 0])
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(2,iStim)],'k')
                plot([iStim, iStim+0.75],[bootCI(2,iStim), bootCI(2,iStim)],'k')
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(1,iStim)],'w')
                plot([iStim, iStim+0.75],[bootCI(1,iStim), bootCI(1,iStim)],'w')
                % Properties for swarm plot
                properties.Orientation =        'vertical';          %(Set how the box plot should be oriented 'vertical' or 'horizontal')
                properties.MarkerType =         'o';                 %(Marker type)
                properties.MarkerFaceColor =    'k';                 %(Marker face color)
                properties.MarkerEdgeColor =    'none';              %(Marker edge color)
                properties.MarkerSize =         5;                   %(Marker size)
                % Only color is changeing
                properties.MarkerFaceColor = SET.Colors(iStim,:);
                % Plot swarm
                FourCageDecision_SubFcn.beeswarmplot_advanced(AllData(:,iStim), iStim+0.75/2, 0.75/2, properties)
            end%iStim
            % Cosmetics
            xlim([0.5,length(bootMean)+0.75+0.5])
            xticks([1+0.75/2 : length(bootMean)+0.75/2])
            set(gca, 'XTickLabel', {'ctr', 'food', 'colony', 'mixture'})
            ylim([0 1])
            ylabel('time at cage (%)', 'interpreter', 'none')
            title(SET.Conditions.cages{iCond})
            export_fig(['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'_PropTimeAtCage'], '-pdf', '-painters')
            % Clean
            close(hFig)
            clear AllData bootCI bootMean currData hFig iStim properties

            % ***** CREATE NumberOfVisits PLOTS *****
            % Collect data
            currData = PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).NumberOfVisits;
            AllData = [[currData.Ctr]', [currData.Leaves]', [currData.Lo]', [currData.LoLeaves]'];
            % Bootstrapping
            bootMean = mean(bootstrp(5000,@nanmean,AllData));
            bootCI = bootci(5000,@nanmean,AllData);
            % Figure
            hFig = figure('Color', 'w');
            hold on
            % Line at zero
            plot([0.75,length(bootMean)+0.75+0.25],[0 0],'k')
            % Iterate over ythe 4 cages
            for iStim = 1:length(bootMean)
                % Plot mean+-CI as bar plot
                rectangle('Position', [iStim, 0, 0.75, bootMean(iStim)], 'EdgeColor', 'none', 'FaceColor', [0 0 0])
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(2,iStim)],'k')
                plot([iStim, iStim+0.75],[bootCI(2,iStim), bootCI(2,iStim)],'k')
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(1,iStim)],'w')
                plot([iStim, iStim+0.75],[bootCI(1,iStim), bootCI(1,iStim)],'w')
                % Properties for swarm plot
                properties.Orientation =        'vertical';          %(Set how the box plot should be oriented 'vertical' or 'horizontal')
                properties.MarkerType =         'o';                %(Marker type)
                properties.MarkerFaceColor =    'k';                 %(Marker face color)
                properties.MarkerEdgeColor =    'none';              %(Marker edge color)
                properties.MarkerSize =         5;                   %(Marker size)
                % Only color is changeing
                properties.MarkerFaceColor = SET.Colors(iStim,:);
                % Plot swarm
                FourCageDecision_SubFcn.beeswarmplot_advanced(AllData(:,iStim), iStim+0.75/2, 0.75/2, properties)
            end%iStim
            % Cosmetics
            xlim([0.5,length(bootMean)+0.75+0.5])
            xticks([1+0.75/2 : length(bootMean)+0.75/2])
            set(gca, 'XTickLabel', {'ctr', 'food', 'colony', 'mixture'})
            ylabel('number of visits', 'interpreter', 'none')
            title(SET.Conditions.cages{iCond})
            export_fig(['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'_NumberOfVisits'], '-pdf', '-painters')
            % Clean
            close(hFig)
            clear AllData bootCI bootMean currData hFig iStim properties

            % ***** CREATE PizzaTime PLOTS *****
            % Collect data
            currData = PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaTime;
            AllData = [[currData.Ctr]', [currData.Leaves]', [currData.Lo]', [currData.LoLeaves]']*100;
            % Bootstrapping
            bootMean = mean(bootstrp(5000,@nanmean,AllData));
            bootCI = bootci(5000,@nanmean,AllData);
            % Figure
            hFig = figure('Color', 'w');
            hold on
            % Line at zero
            plot([0.75,length(bootMean)+0.75+0.25],[0 0],'k')
            % Iterate over ythe 4 cages
            for iStim = 1:length(bootMean)
                % Plot mean+-CI as bar plot
                rectangle('Position', [iStim, 0, 0.75, bootMean(iStim)], 'EdgeColor', 'none', 'FaceColor', [0 0 0])
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(2,iStim)],'k')
                plot([iStim, iStim+0.75],[bootCI(2,iStim), bootCI(2,iStim)],'k')
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(1,iStim)],'w')
                plot([iStim, iStim+0.75],[bootCI(1,iStim), bootCI(1,iStim)],'w')
                % Properties for swarm plot
                properties.Orientation =        'vertical';          %(Set how the box plot should be oriented 'vertical' or 'horizontal')
                properties.MarkerType =         'o';                %(Marker type)
                properties.MarkerFaceColor =    'k';                 %(Marker face color)
                properties.MarkerEdgeColor =    'none';              %(Marker edge color)
                properties.MarkerSize =         5;                   %(Marker size)
                % Only color is changeing
                properties.MarkerFaceColor = SET.Colors(iStim,:);
                % Plot swarm
                FourCageDecision_SubFcn.beeswarmplot_advanced(AllData(:,iStim), iStim+0.75/2, 0.75/2, properties)
            end%iStim
            % Cosmetics
            xlim([0.5,length(bootMean)+0.75+0.5])
            xticks([1+0.75/2 : length(bootMean)+0.75/2])
            set(gca, 'XTickLabel', {'ctr', 'food', 'colony', 'mixture'})
            ylabel('prop. time on slice (%)', 'interpreter', 'none')
            title(SET.Conditions.cages{iCond})
            export_fig(['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'_PizzaTime'], '-pdf', '-painters')
            % Clean
            close(hFig)
            clear AllData bootCI bootMean currData hFig iStim properties

            % ***** CREATE PizzaDensiy PLOTS *****
            % Collect data
            currData = PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity;
            AllData = [[currData.Ctr]', [currData.Leaves]', [currData.Lo]', [currData.LoLeaves]'];            
            % Get relative values
            AllData = AllData./sum(AllData,2);
            % Bootstrapping
            bootMean = bootstrp(5000, @nanmean, AllData);
            bootCI = bootci(5000, {@nanmean, AllData});
            % Pool this for dendorgram
            Data4Dendro.(SET.Phase{iPhase})(iCond,:) = mean(bootMean);
            Data4Dendro_std.(SET.Phase{iPhase})(iCond,:) = std(AllData);
            bootMean = mean(bootMean);
            % Figure
            hFig = figure('Color', 'w');
            hold on
            % Line at zero
            plot([0.75,length(bootMean)+0.75+0.25],[0 0],'k')
            % Iterate over ythe 4 cages
            for iStim = 1:length(bootMean)
                % Plot mean+-CI as bar plot
                rectangle('Position', [iStim, 0, 0.75, bootMean(iStim)], 'EdgeColor', 'none', 'FaceColor', [0 0 0])
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(2,iStim)],'k')
                plot([iStim, iStim+0.75],[bootCI(2,iStim), bootCI(2,iStim)],'k')
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(1,iStim)],'w')
                plot([iStim, iStim+0.75],[bootCI(1,iStim), bootCI(1,iStim)],'w')
                % Properties for swarm plot
                properties.Orientation =        'vertical';          %(Set how the box plot should be oriented 'vertical' or 'horizontal')
                properties.MarkerType =         'o';                %(Marker type)
                properties.MarkerFaceColor =    'k';                 %(Marker face color)
                properties.MarkerEdgeColor =    'none';              %(Marker edge color)
                properties.MarkerSize =         5;                   %(Marker size)
                % Only color is changeing
                properties.MarkerFaceColor = SET.Colors(iStim,:);
                % Plot swarm
                FourCageDecision_SubFcn.beeswarmplot_advanced(AllData(:,iStim), iStim+0.75/2, 0.75/2, properties)
            end%iStim
            % Cosmetics
            xlim([0.5,length(bootMean)+0.75+0.5])
            ylim([0 1])
            xticks([1+0.75/2 : length(bootMean)+0.75/2])
            set(gca, 'XTickLabel', {'ctr', 'food', 'colony', 'mixture'})
            ylabel('avg. relative animal density on slice', 'interpreter', 'none')
            title(SET.Conditions.cages{iCond})
            export_fig(['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'_PizzaDensity'], '-pdf', '-painters')
            % Clean
            close(hFig)
            clear AllData bootCI bootMean currData hFig iStim properties

            % ***** CREATE PreferenceDensiy PLOTS *****
            % Collect data
            currData = PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).PizzaDensity;
            AllData = [[currData.Ctr]'+[currData.Leaves]', [currData.Lo]'+[currData.LoLeaves]'];
            % Get preference index (density at locust-containing quarters minus density
            % at quarters that do not contian locusts)
            AllData = (AllData(:,2)-AllData(:,1))./(AllData(:,2)+AllData(:,1));
            % Bootstrapping
            bootMean = bootstrp(5000, @nanmean, AllData);
            bootMean = mean(bootMean);
            bootCI = bootci(5000, {@nanmean, AllData});
            % Stats
            Stats.PreferenceDensiy.N_Boot = 5e7;
            Stats.PreferenceDensiy.seed = 1234;
            Stats.PreferenceDensiy.nComparisons = 1;
            % Test whether the distribution is different from zero
            [...
                Stats.PreferenceDensiy.(SET.Conditions.cages{iCond}).p,...
                Stats.PreferenceDensiy.(SET.Conditions.cages{iCond}).s,...
                ~,...
                Stats.PreferenceDensiy.(SET.Conditions.cages{iCond}).c] = FourCageDecision_SubFcn.BootstrapHypothesisTesting('one-sample', AllData, 0, Stats.PreferenceDensiy.N_Boot, Stats.PreferenceDensiy.seed, Stats.PreferenceDensiy.nComparisons);
            save('Stats.mat', 'Stats')
            % Figure
            hFig = figure('Color', 'w');
            hold on
            % Line at zero
            plot([0.75,length(bootMean)+0.75+0.25],[0 0],'k')
            % Iterate over ythe 4 cages
            for iStim = 1:length(bootMean)
                % Plot mean+-CI as bar plot
                if bootMean(iStim)>0
                rectangle('Position', [iStim, 0, 0.75, bootMean(iStim)], 'EdgeColor', 'none', 'FaceColor', [0 0 0])
                else
                    rectangle('Position', [iStim, bootMean(iStim), 0.75, abs(bootMean(iStim))], 'EdgeColor', 'none', 'FaceColor', [0 0 0])
                end
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(2,iStim)],'k')
                plot([iStim, iStim+0.75],[bootCI(2,iStim), bootCI(2,iStim)],'k')
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(1,iStim)],'w')
                plot([iStim, iStim+0.75],[bootCI(1,iStim), bootCI(1,iStim)],'w')
                % Properties for swarm plot
                properties.Orientation =        'vertical';          %(Set how the box plot should be oriented 'vertical' or 'horizontal')
                properties.MarkerType =         'o';                %(Marker type)
                properties.MarkerFaceColor =    'k';                 %(Marker face color)
                properties.MarkerEdgeColor =    'none';              %(Marker edge color)
                properties.MarkerSize =         5;                   %(Marker size)
                % Only color is changeing
                properties.MarkerFaceColor = SET.phaseColors.(SET.Phase{iPhase});
                % Plot swarm
                FourCageDecision_SubFcn.beeswarmplot_advanced(AllData(:,iStim), iStim+0.75/2, 0.75/2, properties)
            end%iStim
            % Cosmetics
            xlim([0.5,length(bootMean)+0.75+0.5])
            ylim([-1 1])
            xticks([1+0.75/2 : length(bootMean)+0.75/2])
            set(gca, 'XTickLabel', ['p=',num2str(Stats.PreferenceDensiy.(SET.Conditions.cages{iCond}).p)])
            ylabel('PI', 'interpreter', 'none')
            title(SET.Conditions.cages{iCond})
            export_fig(['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'_PreferenceDensiy'], '-pdf', '-painters')
            % Clean
            close(hFig)
            clear AllData bootCI bootMean currData hFig iStim properties

            % ***** CREATE CageDensiy PLOTS *****
            % Collect data
            currData = PooledDATA.(SET.Phase{iPhase}).(SET.Conditions.cages{iCond}).CageDensity;
            AllData = [[currData.Ctr]', [currData.Leaves]', [currData.Lo]', [currData.LoLeaves]'];
            % Bootstrapping
            bootMean = mean(bootstrp(5000, @nanmean, AllData));
            bootCI = bootci(5000, {@nanmean, AllData});
            % Figure
            hFig = figure('Color', 'w');
            hold on
            % Line at zero
            plot([0.75,length(bootMean)+0.75+0.25],[0 0],'k')
            % Iterate over ythe 4 cages
            for iStim = 1:length(bootMean)
                % Plot mean+-CI as bar plot
                rectangle('Position', [iStim, 0, 0.75, bootMean(iStim)], 'EdgeColor', 'none', 'FaceColor', [0 0 0])
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(2,iStim)],'k')
                plot([iStim, iStim+0.75],[bootCI(2,iStim), bootCI(2,iStim)],'k')
                plot([iStim+0.75/2, iStim+0.75/2],[bootMean(iStim), bootCI(1,iStim)],'w')
                plot([iStim, iStim+0.75],[bootCI(1,iStim), bootCI(1,iStim)],'w')
                % Properties for swarm plot
                properties.Orientation =        'vertical';          %(Set how the box plot should be oriented 'vertical' or 'horizontal')
                properties.MarkerType =         'o';                 %(Marker type)
                properties.MarkerFaceColor =    'k';                 %(Marker face color)
                properties.MarkerEdgeColor =    'none';              %(Marker edge color)
                properties.MarkerSize =         5;                   %(Marker size)
                % Only color is changeing
                properties.MarkerFaceColor = SET.Colors(iStim,:);
                % Plot swarm
                FourCageDecision_SubFcn.beeswarmplot_advanced(AllData(:,iStim), iStim+0.75/2, 0.75/2, properties)
            end%iStim
            % Cosmetics
            xlim([0.5,length(bootMean)+0.75+0.5])
            ylim([0 0.7])
            xticks([1+0.75/2 : length(bootMean)+0.75/2])
            set(gca, 'XTickLabel', {'ctr', 'food', 'colony', 'mixture'})
            ylabel('avg. animal density at cage', 'interpreter', 'none')
            title(SET.Conditions.cages{iCond})
            export_fig(['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iCond},'_CageDensity'], '-pdf', '-painters')
            % Clean
            close(hFig)
            clear AllData bootCI bootMean currData hFig iStim properties

        end%if there is content
    end%iCond

    % ***** Dendrogram for similarity between conditions *****
    hFig = figure('Color', 'w');
    idx = find(any(~isnan(Data4Dendro.(SET.Phase{iPhase})),2));
    Y = pdist(Data4Dendro.(SET.Phase{iPhase})(idx,:), 'euclidean');
    Z = linkage(Y,'single');
    leafOrder = optimalleaforder(Z,Y);
    dendrogram(Z, 'Reorder', leafOrder)
    % set(gca, 'YTick', 0:7)
    set(gca, 'XTickLabels', SET.Conditions.cages(idx(str2num(get(gca, 'XTickLabel')))), 'view', [90 90])
    ylabel(['euclidean distance | c=', num2str(round(cophenet(Z,Y),2))])
    export_fig(['Output\', SET.Phase{iPhase},'_dendro'], '-pdf', '-painters')
    close(hFig)
    clear Y Z leafOrder
    % Save as table
    dendroTable = array2table(Data4Dendro.(SET.Phase{iPhase}), 'RowNames', SET.Conditions.cages, 'VariableNames', SET.Cages);
    writetable(dendroTable, ['Output/Data4Dendro_', SET.Phase{iPhase}, '.csv'])
    dendroTable = array2table(Data4Dendro_std.(SET.Phase{iPhase}), 'RowNames', SET.Conditions.cages, 'VariableNames', SET.Cages);
    writetable(dendroTable, ['Output/Data4Dendro_std_', SET.Phase{iPhase}, '.csv'])

end%iPhase

