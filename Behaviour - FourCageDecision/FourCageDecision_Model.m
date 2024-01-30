%% FourCageDecision_Model
% Applies a Bayesian decision-making model to the 4cage behavioral data in
% order to asses  the quality/reliability of different information classes.
% For these experiments, a single animal (gregarious or solitarious) was
% places in a circular arena (d=90cm) with four stimulus cages, placed on
% opposing sides. The four cages were (i) an empty control Ctr, (ii) a cage
% with fresh blackberry leaves Lvs, (iii) a cage with eitgh gregarious 5th
% instar hoppers Lct, and (iv) a cage with fresh blackberry leaves and
% eight gregarious 5th instar hoppers LvsLct. These cages were either
% transparent with holes (access to both visual and olfactory information),
% opaque (blocked visual inforamtion), or sealed (blocked olfactory
% inforamtion).
% The decision-making rule is an adaption of a Bayesian estimation for the
% integration of different information classes (Gunzel et al. 2023
% iScience, Arganda et al. 2012 PNAS). Here, we predict the proportion of
% animal density within one quarter of the arena based on socio-visual,
% socio-olfactory, nutri-visual, and nutri-olfactory information. For this,
% we first calculate to probability of a quarter beeing a 'good' choice
% based on each information class, e.g. socio-visual information:
%
%    P(A is good | socio-visual) = 1/(1+SV^(-w))
%        where:
%        SV is the quality of socio-visual information [1e-1 : 1e1]
%        w is the availability of this information class [0 or 1]
%
% Now, the predicted proportion of animal density in this quarter is the
% product of all i inforamtion classes for this quarter divided by the sum
% of products for the other 3 (j) quarters:
%
%     P_A = Product[i](P(A is good|C_i)) / Sum[j](Product[i](P(Q_j is good|C_i)))
%
% To fit the model, we find optimal values for the quality of inforamtion
% parameters (e.g., SV in the example above).
%
% Version:
% 20-April-2023 (R2023a)

% Prepare
clc; clear all; close all
mkdir('Output')
warning('off')

% Start processing pool
pPool = gcp;
pPool.IdleTimeout = 720;

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig) 
addpath(genpath('...\GitHub\export_fig'))

%% Settings

% Path to the data
SET.BasePath = '...';

% Phases that should be analyzed
SET.Phase = {'gregarious', 'solitarious'};

% The cage conditions that should be taken into account
SET.Conditions.cages = {...
    'VpOpFp',... visual cues (see-through),   feeding
    'VmOpFp',... no visual cues(opaque),      feeding
    'VpOmFp',... no olfactory cues (glass),   feeding
    };

% Hard-coded, which rows in the data correspond to the above-mentioned cage
% conditions
SET.Conditions.cages_rows = 1:2:6;

% Provide labels for the four different cages. Note, that the order has to
% stay the same
SET.Cages = {...
    'Ctr',...        empty control
    'Lvs',...        blackberry leaves
    'Lct'...         group of locusts
    'LvsLct',...     group of locusts + blackberry leaves
    };

% Settings for parameter fitting
% --- Number of global iterations
Settings.nRep_global = 5000;
% --- Iterations for bayesopt
Settings.nRep_fit = 100;
% --- Error measure maxRMSE, RMSE, maxRSS, RSS, cosine, correlation
SET.dist = 'maxRSS';

%% Model Settings

% Optimizable parameters
SET.bayesfit.socio_visual =    optimizableVariable('socio_visual',    [1e-1 1e1], 'Transform', 'log'); % socio-visual 
SET.bayesfit.nutri_visual =    optimizableVariable('nutri_visual',    [1e-1 1e1], 'Transform', 'log'); % nutri-visual 
SET.bayesfit.socio_olfactory = optimizableVariable('socio_olfactory', [1e-1 1e1], 'Transform', 'log'); % socio-olfactory
SET.bayesfit.nutri_olfactory = optimizableVariable('nutri_olfactory', [1e-1 1e1], 'Transform', 'log'); % nutri-olfactory

% Availability matrices
% socio-visual 
w.socio_visual = [...    | Ctr | Lvs | Lct | LvsLct
    0 0 1 1;...   VpOpFp | 0   | 0   | 1   | 1
    0 0 0 0;...   VmOpFp | 0   | 0   | 0   | 0
    0 0 1 1;...   VpOmFp | 0   | 0   | 1   | 1
    ];
% nutri-visual food
w.nutri_visual = [...    | Ctr | Lvs | Lct | LvsLc
    0 1 0 1;...   VpOpFp | 0   | 1   | 0   | 1
    0 0 0 0;...   VmOpFp | 0   | 0   | 0   | 0
    0 1 0 1;...   VpOmFp | 0   | 1   | 0   | 1
    ];
% socio-olfactory
w.socio_olfactory = [... | Ctr | Lvs | Lct | LvsLc
    0 0 1 1;...   VpOpFp | 0   | 0   | 1   | 1
    0 0 1 1;...   VmOpFp | 0   | 0   | 1   | 1
    0 0 0 0;...   VpOmFp | 0   | 0   | 0   | 0
    ];
% nutri-olfactory
w.nutri_olfactory = [... | Ctr | Lvs | Lct | LvsLc
    0 1 0 1;...   VpOpFp | 0   | 1   | 0   | 1
    0 1 0 1;...   VmOpFp | 0   | 1   | 0   | 1
    0 0 0 0;...   VpOmFp | 0   | 0   | 0   | 0
    ];

%% Fit decision model

% Iterate over each phase and fit parameters independently. This allows a
% coparison of which information class is important to which phase
for iPhase = 1:length(SET.Phase)

    % Get data
    % --- avg
    ground_truth_avg = table2array(readtable(['Output/Data4Dendro_', SET.Phase{iPhase}, '.csv']));
    ground_truth_avg = ground_truth_avg(SET.Conditions.cages_rows,:);
    ground_truth_avg = ground_truth_avg./sum(ground_truth_avg,2);
    % --- std
    ground_truth_std = table2array(readtable(['Output/Data4Dendro_std_', SET.Phase{iPhase}, '.csv']));
    ground_truth_std = ground_truth_std(SET.Conditions.cages_rows,:);

    % Draw from distributions with given avg and std to repeat the fitting
    % many times in order to estimate the variability of the fitted
    % parameters
    ground_truth = zeros([size(ground_truth_avg), Settings.nRep_global]);
    for iRow = 1:size(ground_truth_avg, 1)
        for iCol = 1:size(ground_truth_avg, 2)
            ground_truth(iRow, iCol, :) = random('Normal', ground_truth_avg(iRow, iCol), ground_truth_std(iRow, iCol), [Settings.nRep_global, 1]);
        end%iCol
    end%iRow
    ground_truth(ground_truth<0) = 0;
    ground_truth(ground_truth>1) = 1;
    % Normalize
    for iRep = 1:Settings.nRep_global
        ground_truth(:,:,iRep) = ground_truth(:,:,iRep)./sum(ground_truth(:,:,iRep),2);
    end%iRep

    % Run fitting many times
    % --- Preallocation
    Model.(SET.Phase{iPhase}) = zeros(Settings.nRep_global,4);
    model_fit.(SET.Phase{iPhase}) = zeros(size(ground_truth));
    RMSE.(SET.Phase{iPhase}) = zeros(Settings.nRep_global, 1);
    hWait = waitbar(0, 'Please wait ...');
    for iRep = 1:Settings.nRep_global
        waitbar(iRep/Settings.nRep_global, hWait, ['Please wait ...', num2str(iRep),' /', num2str(Settings.nRep_global)])
        % Optimization function
        % --- Pack variables
        vars = [...
            SET.bayesfit.socio_visual,...
            SET.bayesfit.nutri_visual,...
            SET.bayesfit.socio_olfactory,...
            SET.bayesfit.nutri_olfactory];
        % --- Prepare evaluation function
        fun = @(vars) FourCageDecision_SubFcn.evalfun(ground_truth_avg, w, SET.dist,...
            vars.socio_visual,    vars.nutri_visual,...
            vars.socio_olfactory, vars.nutri_olfactory);

        % Fit parameters using the Bayesian optimization algorithm
        result = bayesopt(fun, vars,...
            'Verbose',0,...
            'AcquisitionFunctionName','expected-improvement-plus',...
            'IsObjectiveDeterministic', true,...
            'UseParallel', true,...
            'ExplorationRatio',0.95,...
            'PlotFcn', [],...
            'MinWorkerUtilization', pPool.NumWorkers-1,...
            'MaxObjectiveEvaluations', Settings.nRep_fit);
        result = bestPoint(result);
        Model.(SET.Phase{iPhase})(iRep,:) = table2array(result);

        % socio_visual    nutri_visual    socio_olfactory    nutri_olfactory

        % Get model prediction and error
        [model_fit.(SET.Phase{iPhase})(:,:,iRep), RMSE.(SET.Phase{iPhase})(iRep)] = FourCageDecision_SubFcn.exefun(...
            ground_truth(:,:,iRep),...
            w,...
            SET.dist,...
            Model.(SET.Phase{iPhase})(iRep,1),...
            Model.(SET.Phase{iPhase})(iRep,2),...
            Model.(SET.Phase{iPhase})(iRep,3),...
            Model.(SET.Phase{iPhase})(iRep,4));
    end%iRep

    % Close waitbar
    close(hWait)

    % Save everything
    save('Output\FourCageDecision_Model_Result.mat', 'SET', 'w', 'Model', 'model_fit', 'RMSE')

end%iPhase

%% Visualize results

% ***** Parameters *****
% Get some nice colors
cols = [...
    227, 74, 51;...
    44, 162, 95]/255;
% Iterate over both phases and create a figure for each
for iPhase = 1:length(SET.Phase)

    hFig = figure('color', 'w'); hold on
    box on
    % Log-transform data to have similar scales for attraction and aversion
    vals = log10(Model.(SET.Phase{iPhase}));
    vals_avg = mean(vals);
    vals_ci = bootci(5000, @mean, vals);
    % Iterate over the different information classes and plot the result as
    % a bar
    for iV = 1:length(vals_avg)
        if vals_avg(iV)<0
            rectangle('position',[iV-0.25, vals_avg(iV), 0.5, abs(vals_avg(iV))], 'EdgeColor','none', 'FaceColor',cols(iPhase,:))
        else
            rectangle('position',[iV-0.25, 0, 0.5, vals_avg(iV)], 'EdgeColor','none', 'FaceColor',cols(iPhase,:))
        end
        plot([1 1]*iV, vals_ci(:,iV), 'k')
        plot([-0.25, 0.25]+iV, [1 1]*vals_ci(1,iV), 'k')
        plot([-0.25, 0.25]+iV, [1 1]*vals_ci(2,iV), 'k')
    end%iV
    % Cosmetics
    plot([0.5 4.5], [0 0], 'k')
    plot([2.5 2.5],[-1 1],'k:')
    xticks(1:4)
    xticklabels({'socio-vis','nutri-vis','socio-olf','nutri-olf'})
    ylabel('reliability score [a.u.]')
    xlabel('information class')
    % Save
    export_fig(['Output\FourCageDecision_Model_Confidence_',SET.Phase{iPhase}], '-pdf', '-painters')
    close(hFig)

end%iPhase


% ***** Results *****
% Get some nice colors
SET.Colors = [...
    200,200,200;...Ctr
    102,166,030;...Lvs
    217,095,002;...Lct
    231,041,138;...LvsLct
    ]/255;

% Iterate over both phases and create a figure for each
for iPhase = 1:length(SET.Phase)

    hFig = figure('color', 'w', 'Name', SET.Phase{iPhase}, 'Units','normalized','Position',[0.25,0.25,0.5,0.25]);
    
    % Get the model's results
    model_vals = model_fit.(SET.Phase{iPhase});    
    
    % Iterate over each condition
    subplot_pos = linspace(1,20,length(SET.Conditions.cages));
    for iC = 1:length(SET.Conditions.cages)
        
        % Get the raw data
        raw_data_file = ['Output\', SET.Phase{iPhase},'_',SET.Conditions.cages{iC},'PizzaDensity.csv'];
        raw_data = table2array(readtable(raw_data_file));
        % Get relative values
        raw_data = raw_data./sum(raw_data,2);
        
        % Plot each condition in its own subplot
        subplot(1, length(SET.Conditions.cages), iC); hold on

        % Iterate over the different information classes
        for iV = 1:size(model_vals,2)

            % Get the average and credible interval
            dataMean = mean(raw_data(:,iV));
            dataCI = bootci(5000, {@nanmean, raw_data(:,iV)});
            bootMean = mean(model_vals(iC,iV,:));
            CredInt = [quantile(model_vals(iC,iV,:), 2.5/100), quantile(model_vals(iC,iV,:), 97.5/100)];
            
            % Plot the data average +- CI
            rectangle('Position', [iV-0.35, 0, 0.35, dataMean], 'EdgeColor', 'none', 'FaceColor', [0.251 0.251 0.251])
            plot([iV-0.35 iV],[dataMean dataMean],'color',[0.5 0.5 0.5])
            plot([-0.175 -0.175]+iV, dataCI, 'Color', [0.5 0.5 0.5])
            % Plot model mean+-CI as bar plot
            rectangle('Position', [iV, 0, 0.35, bootMean], 'EdgeColor', 'none', 'FaceColor', [0.251 0.251 0.251])
            plot([iV iV+0.35],[bootMean bootMean],'color',[0.5 0.5 0.5])
            plot([0.175 0.175]+iV, CredInt, 'color', [0.5 0.5 0.5])            
            % Properties for swarm plot
            properties.Orientation =        'vertical';          %(Set how the box plot should be oriented 'vertical' or 'horizontal')
            properties.MarkerType =         'o';                 %(Marker type)
            properties.MarkerFaceColor =    'k';                 %(Marker face color)
            properties.MarkerEdgeColor =    'none';              %(Marker edge color)
            properties.MarkerSize =         3;                   %(Marker size)
            % Only color is changeing
            properties.MarkerFaceColor = SET.Colors(iV,:);
            % Plot swarm
            FourCageDecision_SubFcn.beeswarmplot_advanced(raw_data(:,iV), iV-0.175, 0.35, properties)
        end%iV
        set(gca,'Units', 'centimeters', 'Position', [subplot_pos(iC) 1 2.5 2.2])

        % Cosmetics
        ylim([0 1])
        xlim([0.5 size(model_vals,2)+0.5])
        title({SET.Conditions.cages{iC};...
            ['RMSE: ', num2str(round(sqrt(mean((mean(model_vals(iC,:,:),3)-mean(raw_data)).^2)),2))]})
        ylabel('avg. relative animal density')
        xticks(1:size(model_vals,2))
        xticklabels(SET.Cages)
        set(gca,'TickLabelInterpreter', 'none')

    end%iC
    
    % Save
    export_fig(['Output\FourCageDecision_Model_Result_',SET.Phase{iPhase}], '-pdf', '-painters')
    close(hFig)

end%iPhase


