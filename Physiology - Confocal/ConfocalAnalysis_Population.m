%% Confocal_ConfocalAnalysis_Population
% This script pools the data across all animals for both social phenotypes 
% and imaging depths. Divided into certain stimulus combinations, it plots 
% the respective activity time courses with corresponding swarm/violin 
% plots, the Bliss interaction score analysis results, and - for the cell 
% bodies only - bar plots for the proportion of active regions.
%
% Version:
% 05-Jan-2024 (R2023a)

% Prepare
clc; clear all; close all
warning('off')

% Add toolboxes
% --- A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig)
addpath(genpath('...\GitHub\export_fig'))
mkdir('ConfocalAnalysis_Population')

%% Settings

% Set paths
SET.main_path = '...';

% Give list of phases
SET.phases = {'gregarious'; 'solitarious'};

% Give layer names
SET.layers = {'layer01_top', 'layer02_middle'};

% Give list of stimulus combinations that should be analyzed together. For each
% triplet, we take the AL regions that were respondnig at least one of the
% the respective triplet's stimuli. Note that this results in different
% sample sizes per triplet as we do not have all combinations in all
% animals.
SET.StimCombiList = {...
    {'COL', 'ZHAE', 'ZHAECOL'};...
    {'COL', 'Z3HL', 'Z3HLCOL'};...
    {'COL', 'OCT', 'OCTCOL'};...
    {'ZHAE', 'Z3HL', 'Z3HLZHAE'};...
    {'ZHAE', 'OCT', 'OCTZHAE'};...
    {'COL', 'ZHAE', 'ZHAECOL', 'Z3HL', 'Z3HLCOL', 'OCT', 'OCTCOL', 'Z3HLZHAE', 'OCTZHAE'}};

% Add traces of the solvent
SET.RefStim = 'MOL';

% Set framerate
SET.fps = 1;

% Set the number of frames the final trace should have
SET.totalFrameNumber = 17;
% Set the analysis window. Note this is hard-coded previous knowledge for
% data that has been cut.
SET.AnalysisWindow = [...
    10 11;... layer 01
    8 9;...   layer 02
    ];

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
SET.TC_lim{1} = [-1 6; -1 8; -1 8; -1 8; -1 8; -1 8; -1 8];
SET.TC_lim{2} = [-1 6; -0.5 1.5; -1.5 1; -1 2; -1.5 1.5; -1.5 1.5];

SET.Swarm_lim{1} = [-1 10; -1 10; -1 10; -1 10; -1 10; -1 10; -1 6];
SET.Swarm_lim{2} = [-2 12; -2 2; -2 2; -2 2; -2 2; -2 2];

SET.Bliss_lim{1} = [-10 10; -8 8; -8 8; -8 8; -8 8; -8 8];
SET.Bliss_lim{2} = [-7 7; -1 1; -1 1; -1 1; -1 1; -1 1];
SET.Bliss_pdf_lim{1} = [0, 0.1; 0, 0.15; 0, 0.15; 0, 0.1; 0, 0.15; 0, 0.12];
SET.Bliss_pdf_lim{2} = [0, 0.25; 0, 0.25; 0, 0.25; 0, 0.25; 0, 0.25; 0, 0.25];

%% Get data
% Iterate over all animals and pool data from all granules

% Iterate over both phases
for iPhase = 1:length(SET.phases)
    % Iterate over all layers
    for iLayer = 1:length(SET.layers)
        % Iterate over all stimulus triplets
        for iTrip = 1:length(SET.StimCombiList)

            % Skip all but the first for glomeruli
            if (iLayer) == 2 && (iTrip > 1)
                break
            end

            % Get overview of animal folders
            path2animal = [SET.main_path, SET.phases{iPhase},'\',SET.layers{iLayer},'\'];
            curr.dir.all = dir(path2animal);

            % Prepare to pool a lot of data
            poolPropRegions = [];
            poolPropMixtureRegions = [];
            poolVennProp = [];
            poolVennProp_info = [];
            poolTC = [];
            poolAvg = [];
            poolBinaryBliss = [];
            poolAvgPosBliss = [];
            poolAvgNegBliss = [];
            poolPdfBliss = [];
            poolInfo = [];

            % Iterate over all animals
            for iAni = 1:size(curr.dir.all,1)
                % Check whether this is what we are looking for, i.e.
                % whether the string "Animal" is present
                if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir

                    % Get the path to the current animal
                    curr.dir.animal = [path2animal,curr.dir.all(iAni).name,'\'];

                    % Get meta data and segmentation
                    load([curr.dir.animal,'03_Data_Processed\meta_info.mat']);
                    Segmentation = load([curr.dir.animal,'03_Data_Processed\img_segmentation.mat']);
                    Segmentation = Segmentation.Segmentation;

                    % Get time vector to cut data
                    time_vec = meta.baseline(1):(meta.baseline(1)+SET.totalFrameNumber-1);

                    % Check whether the current triplet was tested
                    all_stim_available = 1;
                    for iStim = 1:length(SET.StimCombiList{iTrip})
                        if isempty(find(strcmp( meta.unique_stim, SET.StimCombiList{iTrip}{iStim})))
                            all_stim_available = 0;
                            disp([SET.StimCombiList{iTrip}{iStim}, ' missing in ', curr.dir.all(iAni).name])
                            break
                        end%if
                    end%iStim
                    % Continue if everything is available
                    if all_stim_available

                        % Get the list of stimuli
                        StimList = SET.StimCombiList{iTrip};
                        % Check whether the solvent was tested as well
                        if any(strcmp( meta.unique_stim, SET.RefStim))
                            StimList{end+1} = SET.RefStim;
                        end%if solvent

                        % Get list of valid pockets. Take all pockets that
                        % were active at least once.
                        mask_list = zeros(length(SET.StimCombiList{iTrip}), 1);
                        for iStim = 1:length(SET.StimCombiList{iTrip})
                            mask_list(iStim) = find(strcmp( meta.unique_stim, SET.StimCombiList{iTrip}{iStim}));
                            % Load data
                            currData = load([curr.dir.animal,'03_Data_Processed\',StimList{iStim},'.mat']);
                            currData = currData.ImageStream.(StimList{iStim});
                            currData = currData(:,:,time_vec);
                            [x, y, t] = size(currData);
                            currData = reshape(currData, [x*y, t]);
                            % Get the average mean projection during the
                            % active window for each pocket
                            pocket_list = unique(Segmentation.pockets_labeled(:));
                            avg = nan(x,y);
                            for iP = 1:length(pocket_list)
                                ind = find(Segmentation.pockets_labeled == pocket_list(iP));
                                avg(ind) = mean(mean(currData(ind,SET.AnalysisWindow(iLayer,:))));
                            end%iP
                            % Get the binary image of active pockets
                            avg(~Segmentation.manualROI) = NaN;
                            Segmentation.pocket_active_img(:,:,mask_list(iStim)) = imbinarize(avg);
                            clear currData x y t ind avg pocket_list
                        end%iStim
                        mask = sum(Segmentation.pocket_active_img(:,:,mask_list),3);
                        binary_mask = mask>0;
                        pocket_list_valid = unique(Segmentation.pockets_labeled.*double(binary_mask));
                        pocket_list_valid(pocket_list_valid==0) = [];

                        % Get the proportion of regions responding to a
                        % given number of stiuli
                        curr_propResp = nan(length(pocket_list_valid), 1);
                        for iP = 1:length(pocket_list_valid)
                            ind = find(Segmentation.pockets_labeled == pocket_list_valid(iP), 1);
                            curr_propResp(iP) = mask(ind);
                        end%iP
                        curr_propResp = hist(curr_propResp, 1:length(SET.StimCombiList{iTrip}));
                        poolPropRegions = [poolPropRegions; curr_propResp(:)'/length(pocket_list_valid)];

                        % Get proportion of regions that responded to the
                        % mixture, but not to the individual components.
                        % Note, this assumes that the third entry of a
                        % triplet is the mixture
                        if length(SET.StimCombiList{iTrip}) == 3
                            % mask = (...
                            %     Segmentation.pocket_active_img(:,:,mask_list(1))==0 &...
                            %     Segmentation.pocket_active_img(:,:,mask_list(2))==0 &...
                            %     Segmentation.pocket_active_img(:,:,mask_list(3))==1);
                            % p = length(unique(Segmentation.pockets_labeled(mask)))/curr_propResp(1);
                            % poolPropMixtureRegions = [poolPropMixtureRegions; p];

                            mask1 = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==1);
                            mask2 = (Segmentation.pocket_active_img(:,:,mask_list(3))==1);
                            p = length(unique(Segmentation.pockets_labeled(mask1))) / length(unique(Segmentation.pockets_labeled(mask2)));
                            poolPropMixtureRegions = [poolPropMixtureRegions; p];

                        else
                            p = [];
                            % --- zhaecol
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(4))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(5))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(6))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(7))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(8))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(9))==0);
                            p = [p, length(unique(Segmentation.pockets_labeled(mask)))/curr_propResp(1)];
                            % --- z3hlcol
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(4))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(5))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(6))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(7))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(8))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(9))==0);
                            p = [p, length(unique(Segmentation.pockets_labeled(mask)))/curr_propResp(1)];
                            % --- octcol
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(4))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(5))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(6))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(7))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(8))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(9))==0);
                            p = [p, length(unique(Segmentation.pockets_labeled(mask)))/curr_propResp(1)];
                            % --- zahez3hl
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(4))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(5))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(6))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(7))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(8))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(9))==0);
                            p = [p, length(unique(Segmentation.pockets_labeled(mask)))/curr_propResp(1)];
                            % --- zaheoct
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(4))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(5))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(6))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(7))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(8))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(9))==1);
                            p = [p, length(unique(Segmentation.pockets_labeled(mask)))/curr_propResp(1)];
                            poolPropMixtureRegions = [poolPropMixtureRegions; p];
                        end%if

                        % Get the logical relationships between the three
                        % stimuli (Venn diagram)
                        if iTrip==1
                            % Lct only
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==0);
                            p1 = length(unique(Segmentation.pockets_labeled(mask)));
                            % Laa only
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==0);
                            p2 = length(unique(Segmentation.pockets_labeled(mask)));
                            % LaaLct only
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==1);
                            p3 = length(unique(Segmentation.pockets_labeled(mask)));
                            % Laa and Lct
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==0);
                            p4 = length(unique(Segmentation.pockets_labeled(mask)));
                            % Laa and LaaLct
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==1);
                            p5 = length(unique(Segmentation.pockets_labeled(mask)));
                            % Lct and LaaLct
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==0 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==1);
                            p6 = length(unique(Segmentation.pockets_labeled(mask)));
                            % Lct and Laa and LaaLct
                            mask = (...
                                Segmentation.pocket_active_img(:,:,mask_list(1))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(2))==1 &...
                                Segmentation.pocket_active_img(:,:,mask_list(3))==1);
                            p7 = length(unique(Segmentation.pockets_labeled(mask)));
                            % Get proportion
                            p = [p1, p2, p3, p4, p5, p6, p7];
                            p = p/sum(p);
                            poolVennProp = [poolVennProp; p];
                            poolVennProp_info = {'Lct only','Laa only','LaaLct only','Laa and Lct','Laa and LaaLct','Lct and LaaLct','Lct and Laa and LaaLct'};
                        end%iTrip

                        % Iterate over all valid pockets of all stimuli
                        curr_poolTC = nan(length(pocket_list_valid)*length(StimList), length(time_vec));
                        curr_poolAvg = nan(length(pocket_list_valid)*length(StimList), 1);
                        curr_poolInfo = cell(length(pocket_list_valid)*length(StimList), 4);
                        cnt=1;
                        for iStim = 1:length(StimList)
                            % Load data
                            currData = load([curr.dir.animal,'03_Data_Processed\',StimList{iStim},'.mat']);
                            currData = currData.ImageStream.(StimList{iStim});
                            [x,y,t] = size(currData);
                            currData = reshape(currData, [x*y, t]);
                            % Cut data
                            currData = currData(:,time_vec);
                            % Iterate over all pockets
                            for iP = 1:length(pocket_list_valid)
                                % Get index positions of the current pocker
                                idx = find(Segmentation.pockets_labeled == pocket_list_valid(iP));
                                % Get timecourse, activity, and info
                                curr_poolTC(cnt,:) = nanmean(currData(idx,:),1);
                                curr_poolAvg(cnt,1) = nanmean(curr_poolTC(cnt, SET.AnalysisWindow(iLayer,1):SET.AnalysisWindow(iLayer,end) ));
                                curr_poolInfo{cnt,1} = iAni;
                                curr_poolInfo{cnt,2} = pocket_list_valid(iP);
                                curr_poolInfo{cnt,3} = iStim;
                                curr_poolInfo{cnt,4} = StimList{iStim};
                                % Counter
                                cnt=cnt+1;
                            end%ip
                        end%iStim

                        % Get the Bliss score to predict the combined effect of
                        % Laa and Lct and if they are acting independently of
                        % each other.
                        % --- Get the correct indices
                        idx_1 = find(strcmp(curr_poolInfo(:,4), SET.StimCombiList{iTrip}{1}));
                        idx_2 = find(strcmp(curr_poolInfo(:,4), SET.StimCombiList{iTrip}{2}));
                        idx_3 = find(strcmp(curr_poolInfo(:,4), SET.StimCombiList{iTrip}{3}));
                        m1 = curr_poolTC(idx_1,:)/100;
                        m2 = curr_poolTC(idx_2,:)/100;
                        m3 = curr_poolTC(idx_3,:)/100;
                        % --- Calculate the Bliss expected effect
                        f_expected = (m1/2 + m2/2) - ((m1/2).*(m2/2));
                        % --- Get the observed values
                        f_observed = m3;
                        % --- Compare the Bliss expected effect with the observed
                        bliss_score = 100*(f_observed - f_expected);
                        % --- Subtract baseline
                        bliss_score = bliss_score - mean(bliss_score(:,meta.baseline-meta.baseline(1)+1),2);
                        % Get Bliss scores in active region
                        active_bliss = mean(bliss_score(:, SET.AnalysisWindow(iLayer,1):SET.AnalysisWindow(iLayer,end) ),2);
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
                    end%if triplet
                end%if animal
            end%iAni

            % Get grand means for TCs and AVGs
            ani_info = cell2mat(poolInfo(:,1));
            animal_list = unique(ani_info);
            grandMean_animal = cell(length(animal_list), 1);
            for iStim = 1:length(StimList)
                grandMean_TC.(StimList{iStim}) = nan(length(animal_list), SET.totalFrameNumber);
                grandMean_AVG.(StimList{iStim}) = nan(length(animal_list), 1);
                for iAni = 1:length(animal_list)
                    idx = find(ani_info==animal_list(iAni) & strcmp(poolInfo(:,4), StimList{iStim}));
                    grandMean_TC.(StimList{iStim})(iAni,:) = nanmean(poolTC(idx,:),1);
                    grandMean_AVG.(StimList{iStim})(iAni,1) = nanmean(poolAvg(idx,:),1);
                    grandMean_animal{iAni,1} = curr.dir.all(animal_list(iAni)).name;
                end%iAni
            end%iStim

            % -------------------------------------------------------------
            % Save everything
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolPropRegions = poolPropRegions;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolPropMixtureRegions = poolPropMixtureRegions;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolVennProp = poolVennProp;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolVennProp_info = poolVennProp_info;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolTC = poolTC;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolAvg = poolAvg;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolBinaryBliss = poolBinaryBliss;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolAvgPosBliss = poolAvgPosBliss;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolAvgNegBliss = poolAvgNegBliss;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolPdfBliss = poolPdfBliss;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).poolInfo = poolInfo;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).grandMean_TC = grandMean_TC;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).grandMean_AVG = grandMean_AVG;
            PooledData.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).animal = grandMean_animal;
            save('ConfocalAnalysis_Population\PooledData.mat', 'PooledData', 'SET')

            % -------------------------------------------------------------
            % Plot time courses
            hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 1 1]);
            xvec = 0:(1/SET.fps):((SET.totalFrameNumber/SET.fps)-(1/SET.fps));
            hold on
            % Indicate location of active window
            rectangle('position', [xvec(SET.AnalysisWindow(iLayer,1)), 0, xvec(SET.AnalysisWindow(iLayer,end))-xvec(SET.AnalysisWindow(iLayer,1)), 1],...
                'FaceColor',[0.5 0.5 0.5],...
                'EdgeColor','none')
            % Iterate over all stimuli and plot mean +/- CI
            for iStim = 1:length(StimList)
                % Get the right color index
                col_idx = find(strcmp(SET.ColorNames, StimList{iStim}));
                % Get grand mean and CI
                avg = nanmean(grandMean_TC.(StimList{iStim}),1);
                if size(grandMean_TC.(StimList{iStim}),1) == 1
                    ci = nan(2, SET.totalFrameNumber);
                elseif size(grandMean_TC.(StimList{iStim}),1) == 2
                    ci = grandMean_TC.(StimList{iStim});
                else
                    ci = bootci(5000, {@nanmean, grandMean_TC.(StimList{iStim})});
                end
                % Plot
                plot(xvec, avg, 'color', SET.Colors(col_idx,:))
                plot(xvec, ci, 'color', SET.Colors(col_idx,:))
            end%iStim
            % Cosmetics
            axis square
            ylabel('\DeltaF/F (%)')
            xlabel('time (s)')
            ylim(SET.TC_lim{iLayer}(iTrip,:))
            xlim([xvec(1), xvec(end)])
            title({[SET.phases{iPhase},' | ', SET.layers{iLayer}];...
                strjoin(SET.StimCombiList{iTrip}, ' | ')}, 'Interpreter','none')
            fig_name = ['ConfocalAnalysis_Population\', SET.phases{iPhase}, '_', SET.layers{iLayer}, '_', 'TCs_StimCombiList', num2str(iTrip)];
            export_fig(fig_name, '-pdf', '-painters')
            close(hFig)

            % ---------------------------------------------------------------------
            % Plot violin/swarm plots
            hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 1 1]);
            % Iterate over all stimuli
            for iStim = 1:length(StimList)
                subplot(2,length(StimList), iStim); hold on
                % Get the right color index
                col_idx = find(strcmp(SET.ColorNames, StimList{iStim}));
                % Violin
                vioData = grandMean_AVG.(StimList{iStim});
                properties.MinVal = nanmin(vioData);
                properties.MaxVal = nanmax(vioData);
                properties.AvgType = 'mean';
                properties.MeanCol = 'k';
                properties.MeanSymbol = 'line';
                properties.MeanWidth = 2;
                properties.SeparateOutliers = 0;
                Confocal_SubFcn.violinplot_advanced(vioData, 0, 0.5, properties)
                plot([0,0], [nanmin(vioData), nanmax(vioData)], 'k')
                clear vioData properties
                % Swarm
                beeData = grandMean_AVG.(StimList{iStim});
                properties.MarkerType = 'o';
                properties.MarkerFaceColor = SET.Colors(col_idx,:);
                properties.MarkerSize = 5;
                Confocal_SubFcn.beeswarmplot_advanced(beeData, 0.25, 0.5, properties)
                % Add mean and CI
                avg = nanmean(beeData);
                plot([0.125 0.375], [avg avg], 'k')
                if length(beeData)>2
                    ci = bootci(5000, {@nanmean, beeData});
                    plot([0.25 0.25], ci, 'k')
                end
                clear beeData properties
                % Cosmetics
                ylim(SET.Swarm_lim{iLayer}(iTrip,:))
                xlim([-0.5 0.75])
                xticks([])
                ylabel('\DeltaF/F (%)')
                xlabel(StimList{iStim}, 'Interpreter','none')
                axis square
                % Connecting lines
                if iStim ~= length(StimList)
                    subplot(2, length(StimList), iStim+length(StimList)); hold on
                    plot([1 2], [grandMean_AVG.(StimList{iStim})(:), grandMean_AVG.(StimList{iStim+1})(:)], 'k')
                    % Cosmetics
                    ylim(SET.Swarm_lim{iLayer}(iTrip,:))
                    xlim([0.5 2.5])
                    xticks([])
                    ylabel('\DeltaF/F (%)')
                    xlabel([StimList{iStim},'-',StimList{iStim+1}], 'Interpreter','none')
                    title(signrank(grandMean_AVG.(StimList{iStim}), grandMean_AVG.(StimList{iStim+1})))
                    axis square
                end%if
            end%iStim
            fig_name = ['ConfocalAnalysis_Population\', SET.phases{iPhase}, '_', SET.layers{iLayer}, '_', 'VioSwarm_StimCombiList', num2str(iTrip)];
            export_fig(fig_name, '-pdf', '-painters')
            close(hFig)

            % ---------------------------------------------------------------------
            % Plot Bliss score
            hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 0.5 1]);

            subplot(2,2,1); hold on
            % Violin
            properties.MinVal = nanmin(poolBinaryBliss);
            properties.MaxVal = nanmax(poolBinaryBliss);
            properties.AvgType = 'mean';
            properties.MeanCol = 'k';
            properties.MeanSymbol = 'line';
            properties.MeanWidth = 2;
            properties.SeparateOutliers = 0;
            Confocal_SubFcn.violinplot_advanced(poolBinaryBliss, 0, 0.5, properties)
            plot([0,0], [min(poolBinaryBliss), max(poolBinaryBliss)], 'k')
            clear properties
            % Swarm
            properties.MarkerType = 'o';
            properties.MarkerFaceColor = SET.PhaseColors.(SET.phases{iPhase});
            properties.MarkerSize = 5;
            Confocal_SubFcn.beeswarmplot_advanced(poolBinaryBliss, 0.25, 0.5, properties)
            % Add mean and CI
            avg = nanmean(poolBinaryBliss);
            plot([0.125 0.375], [avg avg], 'k')
            if length(poolBinaryBliss)>2
                ci = bootci(5000, {@nanmean, poolBinaryBliss});
                plot([0.25 0.25], ci, 'k')
            end
            clear properties
            % --- Cosmetics
            xlim([-0.5 0.75])
            ylim([0 1])
            ylabel('prop. syn. gran.')
            axis square

            subplot(2,2,2); hold on
            % Violin
            properties.MinVal = nanmin(poolAvgPosBliss);
            properties.MaxVal = nanmax(poolAvgPosBliss);
            properties.AvgType = 'mean';
            properties.MeanCol = 'k';
            properties.MeanSymbol = 'line';
            properties.MeanWidth = 2;
            properties.SeparateOutliers = 0;
            Confocal_SubFcn.violinplot_advanced(poolAvgPosBliss, 0, 0.5, properties)
            plot([0,0], [min(poolAvgPosBliss), max(poolAvgPosBliss)], 'k')
            clear properties
            % Swarm
            properties.MarkerType = 'o';
            properties.MarkerFaceColor = SET.PhaseColors.(SET.phases{iPhase});
            properties.MarkerSize = 5;
            Confocal_SubFcn.beeswarmplot_advanced(poolAvgPosBliss, 0.25, 0.5, properties)
            % Add mean and CI
            avg = nanmean(poolAvgPosBliss);
            plot([0.125 0.375], [avg avg], 'k')
            if length(poolAvgPosBliss)>2
                ci = bootci(5000, {@nanmean, poolAvgPosBliss});
                plot([0.25 0.25], ci, 'k')
            end
            clear properties
            axis square

            subplot(2,2,2); hold on
            % Violin
            properties.MinVal = min(poolAvgNegBliss);
            properties.MaxVal = max(poolAvgNegBliss);
            properties.AvgType = 'mean';
            properties.MeanCol = 'k';
            properties.MeanSymbol = 'line';
            properties.MeanWidth = 2;
            properties.SeparateOutliers = 0;
            Confocal_SubFcn.violinplot_advanced(poolAvgNegBliss, 0, 0.5, properties)
            plot([0,0], [min(poolAvgNegBliss), max(poolAvgNegBliss)], 'k')
            clear properties
            % Swarm
            properties.MarkerType = 'o';
            properties.MarkerFaceColor = SET.PhaseColors.(SET.phases{iPhase});
            properties.MarkerSize = 5;
            Confocal_SubFcn.beeswarmplot_advanced(poolAvgNegBliss, 0.25, 0.5, properties)
            % Add mean and CI
            avg = nanmean(poolAvgNegBliss);
            plot([0.125 0.375], [avg avg], 'k')
            if length(poolAvgNegBliss)>2
                ci = bootci(5000, {@nanmean, poolAvgNegBliss});
                plot([0.25 0.25], ci, 'k')
            end
            clear properties
            % --- Cosmetics
            xlim([-0.5 0.75])
            ylim(SET.Bliss_lim{iLayer}(iTrip,:))
            ylabel('avg. Bliss of syn gran')
            axis square

            % Average of all animals during active regions
            subplot(2,2,3); hold on
            % --- Get data
            StandBliss = (poolAvgPosBliss + poolAvgNegBliss) ./ (poolAvgPosBliss + abs(poolAvgNegBliss));
            % Violin
            properties.MinVal = nanmin(StandBliss);
            properties.MaxVal = nanmax(StandBliss);
            properties.AvgType = 'mean';
            properties.MeanCol = 'k';
            properties.MeanSymbol = 'line';
            properties.MeanWidth = 2;
            properties.SeparateOutliers = 0;
            Confocal_SubFcn.violinplot_advanced(StandBliss, 0, 0.5, properties)
            plot([0,0], [min(StandBliss), max(StandBliss)], 'k')
            clear properties
            % Swarm
            properties.MarkerType = 'o';
            properties.MarkerFaceColor = SET.PhaseColors.(SET.phases{iPhase});
            properties.MarkerSize = 5;
            Confocal_SubFcn.beeswarmplot_advanced(StandBliss, 0.25, 0.5, properties)
            % Add mean and CI
            avg = nanmean(StandBliss);
            plot([0.125 0.375], [avg avg], 'k')
            if length(StandBliss)>2
                ci = bootci(5000, {@nanmean, StandBliss});
                plot([0.25 0.25], ci, 'k')
            end
            clear properties
            % --- Cosmetics
            xlim([-0.5 0.75])
            ylim([-1 1])
            ylabel('standardized Bliss score')
            axis square

            % Average PDF
            subplot(2,2,4); hold on
            % --- Get data
            xvec = linspace(-100,100,1000);
            % --- Plot
            avg = nanmean(poolPdfBliss,1);
            plot(xvec, avg, 'color', SET.PhaseColors.(SET.phases{iPhase}))
            if size(poolPdfBliss,1) == 1
                ci = nan(2, length(avg));
            elseif size(poolPdfBliss,1) == 2
                ci = poolPdfBliss;
            else
                ci = bootci(5000, {@nanmean, poolPdfBliss});
            end
            plot(xvec, ci', 'color', SET.PhaseColors.(SET.phases{iPhase}))
            % --- Cosmetics
            xlim([-15 15])
            ylim(SET.Bliss_pdf_lim{iLayer}(iTrip,:))
            xlabel('Bliss score')
            ylabel('PDF')
            axis square
            fig_name = ['ConfocalAnalysis_Population\', SET.phases{iPhase}, '_', SET.layers{iLayer}, '_', 'BlissScore_StimCombiList', num2str(iTrip)];
            export_fig(fig_name, '-pdf', '-painters')
            close(hFig)

            % ---------------------------------------------------------------------
            % Imaging statistics
            StimList = StimList(~strcmp(StimList, 'MOL'));
            Stats.(SET.phases{iPhase}).seed=1234;
            Stats.(SET.phases{iPhase}).n = 1e7;
            Stats.(SET.phases{iPhase}).nComp = length(combnk(1:length(StimList), 2));
            for iStim = 1:length(StimList)
                for jStim = 1:length(StimList)
                    if iStim>jStim
                        [...
                            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).([StimList{iStim},'_vs_',StimList{jStim}]).p,...
                            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).([StimList{iStim},'_vs_',StimList{jStim}]).s,...
                            ~,...
                            Stats.(SET.phases{iPhase}).(SET.layers{iLayer}).(['StimCombiList', num2str(iTrip)]).([StimList{iStim},'_vs_',StimList{jStim}]).c] = ...
                            Confocal_SubFcn.BootstrapHypothesisTesting('two-sample-pairs',...
                            grandMean_AVG.(StimList{iStim}),...
                            grandMean_AVG.(StimList{jStim}),...
                            Stats.(SET.phases{iPhase}).n,...
                            Stats.(SET.phases{iPhase}).seed,...
                            Stats.(SET.phases{iPhase}).nComp);
                    end%if
                end%jStim
            end%iStim
            % Save
            save(['ConfocalAnalysis_Population', '\', 'statistics.mat'], 'Stats')

            clearvars -except iPhase iLayer iTrip PooledData Stats SET
        end%iTrip
    end%iLayer
end%iPhase


%% Proportion of regions responing to a given number of stimuli
hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 0.5 0.5]);
hold on
% --- greg data
avg_greg = nanmean(PooledData.gregarious.layer01_top.StimCombiList6.poolPropRegions);
ci_greg = bootci(5000, {@nanmean, PooledData.gregarious.layer01_top.StimCombiList6.poolPropRegions});
% --- soli data
avg_soli = nanmean(PooledData.solitarious.layer01_top.StimCombiList6.poolPropRegions);
ci_soli = bootci(5000, {@nanmean, PooledData.solitarious.layer01_top.StimCombiList6.poolPropRegions});
% Prep for stats
Stats.prop_n_regions.seed = 1234;
Stats.prop_n_regions.n = 1e7;
Stats.prop_n_regions.nComp = 1;
% Plot
for iStim  = 1:length(avg_greg)
    % --- greg
    rectangle('Position', [iStim-0.25, 0, 0.25, avg_greg(iStim)], 'FaceColor', [191 0 0]/255, 'EdgeColor', 'none')
    % Swarm
    beeData = PooledData.gregarious.layer01_top.StimCombiList6.poolPropRegions(:,iStim);
    properties.MarkerType = 'o';
    properties.MarkerFaceColor = [0.21 0.21 0.21];
    properties.MarkerSize = 5;
    Confocal_SubFcn.beeswarmplot_advanced(beeData, iStim-0.125, 0.125, properties)
    plot([1 1]*(iStim-0.25/2), ci_greg(:,iStim), 'k')
    % --- soli
    rectangle('Position', [iStim, 0, 0.25, avg_soli(iStim)], 'FaceColor', [0 0 191]/255, 'EdgeColor', 'none')
    % Swarm
    beeData = PooledData.solitarious.layer01_top.StimCombiList6.poolPropRegions(:,iStim);
    properties.MarkerType = 'o';
    properties.MarkerFaceColor = [0.21 0.21 0.21];
    properties.MarkerSize = 5;
    Confocal_SubFcn.beeswarmplot_advanced(beeData, iStim+0.125, 0.125, properties)
    plot([1 1]*(iStim+0.25/2), ci_soli(:,iStim), 'k')
    % Stats
    [...
        Stats.prop_n_regions.p(iStim),...
        Stats.prop_n_regions.s(iStim),...
        ~,...
        Stats.prop_n_regions.c(iStim)] = ...
        Confocal_SubFcn.BootstrapHypothesisTesting('two-sample',...
        PooledData.gregarious.layer01_top.StimCombiList6.poolPropRegions(:,iStim),...
        PooledData.solitarious.layer01_top.StimCombiList6.poolPropRegions(:,iStim),...
        Stats.prop_n_regions.n,...
        Stats.prop_n_regions.seed,...
        Stats.prop_n_regions.nComp);
end%iStim
ylim([0 0.3])
xlim([0.5 length(avg_greg)+0.5])
xticks(1:length(avg_greg))
export_fig('ConfocalAnalysis_Population\prop_n_regions', '-pdf', '-painters')
close(hFig)
% Save
save(['ConfocalAnalysis_Population', '\', 'statistics.mat'], 'Stats')


%% Proportion of mixture-specific regions

% Get data
% --- greg
greg = {};
avg_greg = [];
ci_greg = [];
greg{1} = PooledData.gregarious.layer01_top.StimCombiList1.poolPropMixtureRegions;
greg{2} = PooledData.gregarious.layer01_top.StimCombiList3.poolPropMixtureRegions;
greg{3} = PooledData.gregarious.layer01_top.StimCombiList4.poolPropMixtureRegions;
greg{4} = PooledData.gregarious.layer01_top.StimCombiList5.poolPropMixtureRegions;
% --- soli
soli = {};
avg_soli = [];
ci_soli = [];
soli{1} = PooledData.solitarious.layer01_top.StimCombiList1.poolPropMixtureRegions;
soli{2} = PooledData.solitarious.layer01_top.StimCombiList3.poolPropMixtureRegions;
soli{3} = PooledData.solitarious.layer01_top.StimCombiList4.poolPropMixtureRegions;
soli{4} = PooledData.solitarious.layer01_top.StimCombiList5.poolPropMixtureRegions;
% Get average and CI
for iStim = 1:length(greg)
    % --- greg
    avg_greg = [avg_greg, nanmean(greg{iStim})];
    ci_greg = [ci_greg, bootci(5000, {@nanmean, greg{iStim}})];
    % --- soli
    avg_soli = [avg_soli, nanmean(soli{iStim})];
    ci_soli = [ci_soli, bootci(5000, {@nanmean, soli{iStim}})];
end%iCombi

% Prepare for stats
Stats.prop_mix_regions.seed=1234;
Stats.prop_mix_regions.n = 1e7;
Stats.prop_mix_regions.nComp = 1;
% Plot
hFig = figure('Color','w', 'Units','normalized', 'Position',[0 0 0.5 0.5]);
hold on
for iStim  = 1:length(avg_greg)
    % --- greg
    rectangle('Position', [iStim-0.25, 0, 0.25, avg_greg(iStim)], 'FaceColor', [191 0 0]/255, 'EdgeColor', 'none')    
    beeData = greg{iStim};
    properties.MarkerType = 'o';
    properties.MarkerFaceColor = [0.21 0.21 0.21];
    properties.MarkerSize = 5;
    Confocal_SubFcn.beeswarmplot_advanced(beeData, iStim-0.125, 0.125, properties)
    plot([1 1]*(iStim-0.25/2), ci_greg(:,iStim), 'k')
    % --- soli
    rectangle('Position', [iStim, 0, 0.25, avg_soli(iStim)], 'FaceColor', [0 0 191]/255, 'EdgeColor', 'none')
    beeData = soli{iStim};
    properties.MarkerType = 'o';
    properties.MarkerFaceColor = [0.21 0.21 0.21];
    properties.MarkerSize = 5;
    Confocal_SubFcn.beeswarmplot_advanced(beeData, iStim+0.125, 0.125, properties)
    plot([1 1]*(iStim+0.25/2), ci_soli(:,iStim), 'k')
    % Stats
    [...
        Stats.prop_mix_regions.p(iStim),...
        Stats.prop_mix_regions.s(iStim),...
        ~,...
        Stats.prop_mix_regions.c(iStim)] = ...
        Confocal_SubFcn.BootstrapHypothesisTesting('two-sample',...
        greg{iStim},...
        soli{iStim},...
        Stats.prop_mix_regions.n,...
        Stats.prop_mix_regions.seed,...
        Stats.prop_mix_regions.nComp);
end%iStim
set(gca, 'XTick', 1:length(avg_greg), 'XTickLabel', {'Laa Lct', 'Oct Lct', 'Z3hl Laa', 'Oct Laa'})
ylim([0 0.5])
xlim([0.5 length(avg_greg)+0.5])
export_fig('ConfocalAnalysis_Population\prop_mix_regions', '-pdf', '-painters')
close(hFig)
% Save
save(['ConfocalAnalysis_Population', '\', 'statistics.mat'], 'Stats')


%% Compare gregarious and solitarious responses to Lct, Laa, and LaaLct

Stats.greg_vs_soli.seed=1234;
Stats.greg_vs_soli.n = 1e7;
Stats.greg_vs_soli.nComp = 1;
for iTrip = 1:5
    % Get list of stimuli
    stim_list = SET.StimCombiList{iTrip};
    % Iterate over all stimuli
    for iStim = 1:length(stim_list)
        [...
            Stats.greg_vs_soli.(['StimCombiList',num2str(iTrip)]).(stim_list{iStim}).p(1),...
            Stats.greg_vs_soli.(['StimCombiList',num2str(iTrip)]).(stim_list{iStim}).s(1),...
            ~,...
            Stats.greg_vs_soli.(['StimCombiList',num2str(iTrip)]).(stim_list{iStim}).c(1)] = ...
            Confocal_SubFcn.BootstrapHypothesisTesting('two-sample',...
            PooledData.gregarious.(SET.layers{1}).(['StimCombiList',num2str(iTrip)]).grandMean_AVG.(stim_list{iStim}),...
            PooledData.solitarious.(SET.layers{1}).(['StimCombiList',num2str(iTrip)]).grandMean_AVG.(stim_list{iStim}),...
            Stats.greg_vs_soli.n,...
            Stats.greg_vs_soli.seed,...
            Stats.greg_vs_soli.nComp);
        if iTrip == 1
            [...
                Stats.greg_vs_soli.(['StimCombiList',num2str(iTrip)]).(stim_list{iStim}).p(2),...
                Stats.greg_vs_soli.(['StimCombiList',num2str(iTrip)]).(stim_list{iStim}).s(2),...
                ~,...
                Stats.greg_vs_soli.(['StimCombiList',num2str(iTrip)]).(stim_list{iStim}).c(2)] = ...
                Confocal_SubFcn.BootstrapHypothesisTesting('two-sample',...
                PooledData.gregarious.(SET.layers{2}).(['StimCombiList',num2str(iTrip)]).grandMean_AVG.(stim_list{iStim}),...
                PooledData.solitarious.(SET.layers{2}).(['StimCombiList',num2str(iTrip)]).grandMean_AVG.(stim_list{iStim}),...
                Stats.greg_vs_soli.n,...
                Stats.greg_vs_soli.seed,...
                Stats.greg_vs_soli.nComp);
        end%if
    end%iStim
end%iTrip
% Save
save(['ConfocalAnalysis_Population', '\', 'statistics.mat'], 'Stats')