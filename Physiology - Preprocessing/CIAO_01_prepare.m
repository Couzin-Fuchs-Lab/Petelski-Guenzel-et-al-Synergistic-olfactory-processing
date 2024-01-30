%% Calcium Imaging Analysis Operator (prepare)
% CIOA_01_prepare is the first of a series of scripts aiming for an
% easy-to-use calcium imaging analysis pipeline. Here, we collect and
% prepare raw data, apply spatial and temporal filters, and correct
% baseline shifts due to, for example, bleaching.
% CIOA only assumes a strict folder structure, and each analysis step is
% saved in a given folder. Thus, each module can be run independently, and
% results from different stages can be used without running a long
% pipeline.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% In the following, we lay out the folder structure needed for CIAO to run
% without errors:
%   > DATA                       : this is the folder you select at the beginning
%   |--> Animal01                : this folder can have any name that contains the string "Animal"
%   |  |--> 01_Data_raw          : CIAO assumes this folder contains the raw data, separated by trials
%   |  |  |--> Trial01           : CIAO assumes this folder contains the text file "protocol.txt" with meta information, together with individual tiff images
%   |  |  |--> Trial##           : same logic for all following trials
%   |  |--> 02_Data_Matlab       : folder will be created automatically and will contain preprocessed data
%   |  |  |--> Trial01           : each stimulus will be saved in the order it was given as a separate mat file
%   |  |  |--> Trial##           : same logic for all following trials
%   |  |--> 03_Data_Processed    : folder will be created automatically and will contain segmented data
%   |  |  |--> Trial01           : each stimulus will be saved in the order it was given as a separate mat file
%   |  |  |--> Trial##           : same logic for all following trials
%   |  |--> 04_Data_Summary      : folder will be created automatically and will contain segmented data
%   |  |  |--> Trial01           : each stimulus will be saved in the order it was given as a separate mat file
%   |  |  |--> Trial##           : same logic for all following trials
%   |--> Animal##                : same logic for all following animals
%   ...
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% CIAO requires the text file "protocol.txt", containing meta information,
% to load all data correctly. For this, three entries are crucial:
% "Stimuli", "S#_on", and "Duration". Under Stimuli, you list all stimulus
% names, separated by a comma and a blank ", ", for the respective trial.
% Under Duration, you provide the length (frames) of a stimulus file. Note,
% that currently CIOA can only process stimuli with equal lengths. If a
% stimulus is shorter or longer, it will be extended with NaNs or cut,
% respectively. Under S#_on, you provide the stimulus onset for each
% stimulus presentation (e.g., S1_on for the first, and S2_on for the
% second). For multiple presentations within the same trial, CIOA will
% use each periode before an onset indidividually. If you do not want this,
% consider splitting it into several trials (see folder structure above).
% An exemplary protocol.txt file could look like the following:
%
%   Date            : yymmdd
%   Animal          : 001
%   Trial           : 01
%   Duration        : 25
%   FPS             : 5
%   Stimuli         : stim01, stim02, stim03
%   S1_on           : 10
%   Note            : example data set
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Last, in oder to load the correct data file, the file name must contain
% the stiulus name separated by underscores. An examplary file name could
% look like the following:
%
%   yymmdd_Animal001_Trial01_stim01_ProtocolName.*
%   yymmdd_Animal001_Trial01_stim02_ProtocolName.*
%
% Note for ratiometric pst data (e.g. R340/380) the numerator stimulus name
% should end with a A, and the denominator stimulus name should end with
% a B, but do not add A and B in the protocol.txt. An examplary file name
% could look like the following:
%
%   yymmdd_Animal001_Trial01_stim01A_ProtocolName.pst
%   yymmdd_Animal001_Trial01_stim01B_ProtocolName.pst
%
% Different imaging protocols produce raw data in different data formats.
% Make sure the correct data format is selected in the section "Settings"
% below. This section is the only part that needs to be customized for your
% dataset.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Version: 23-Nov-23 (R2023a)

%% Clean and prepare the environment
clc; clear all; close all
disp('----- Ciao, welcome back! -----')
% Add paths:
% --- Current folder and subfolders
addpath(genpath(pwd))
% --- Non-rigid motion correction of calcium imaging data
%     https://github.com/flatironinstitute/NoRMCorre
addpath(genpath('...\GitHub\NoRMCorre'))
% --- Toolbox for exporting publication quality figures
%     https://github.com/altmany/export_fig
addpath(genpath('...\GitHub\export_fig'))
% --- Toolbox for nice colors
%     https://github.com/DrosteEffect/BrewerMap
addpath(genpath('...\GitHub\BrewerMap'))
% Start parallel processing
gcp;

%% Settings

% Set the file format the data is saved in. Currently the following formats
% are supported:
% - tiff-stack ('tiff_stack')
% - individual tiff files ('single_tiff')
% - lsm ('lsm')
% - pst ('pst').
SET.fileformat = 'lsm';

% Set target size of images to decrease file size.
SET.TargetSize = [256 256];

% Set the standard deviation for 2D Gaussian filtering ('gauss'), box
% size for 2D median filtering ('med'), or both [box, std] for 
% both ('med-gauss'), or vice versa for 'gauss-med'
% Note that spatial filtering is applied after resizing the image.
SET.SmoothData_spat_type = 'med-gauss';
SET.SmoothData_spat = [3 1];

% Set whether to detrend data. For a course trend correction, we subtract a
% trendline based on a linear fit to each stimulus' raw pixel intensity
% data.
SET.Detrend = true;

% Set whether to correct for motion
SET.CorrectMotion = true;

% Set the maximum number of iterations for image registration during motion
% correction. This will be ignored if SET.CorrectMotion is set to false.
SET.maxRegIter = inf;

% Set for how may iterations NoRMCorre should run
SET.NoRMCorreIter = 5;

% Set depth for 3D box filtering with no spatial filter, i.e. temporal
% fitlering.
SET.SmoothData_temp = 3;

% Set baseline length (frames). The window will end at the stimulus onset
SET.Baseline = 5;

% Set whether to run the script for all animals or just for new ones
SET.overwrite = false;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% No need to touch the code below

%% Process data

% Get the path to raw data
SET.path2animal = uigetdir(pwd,'Select folder listing all animals');
% Get an overview of animal folders
curr.dir.all = dir(SET.path2animal);
% Prepare to show raw traces
rawTraces = figure('Units','normalized','Position', [0.25 0.25 0.5 0.5]);
% Iterate over all animals
for iAni = 1:size(curr.dir.all,1)
    % Check whether this is what we are looking for, i.e. whether the
    % string "Animal" is present
    if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir
        % Get information about trials
        curr.path.animal = [SET.path2animal,'\',curr.dir.all(iAni).name,'\01_Data_raw'];
        curr.dir.animal = dir(curr.path.animal);
        SET.N = 0;
        for iItem = 1:size(curr.dir.animal,1)
            if ~isempty(strfind(curr.dir.animal(iItem).name, 'Trial')) && curr.dir.animal(iItem).isdir
                SET.N = SET.N+1;
            end%if
        end%iItem

        % ----- Load and align frames -----
        % Check whether to run preprocessing or not
        if SET.overwrite || ~isfolder([SET.path2animal,'\',curr.dir.all(iAni).name,'\','02_Data_Matlab'])

            % Load stacks
            % -------------------------------------------------------------
            switch SET.fileformat
                case 'tiff_stack'
                    [currStack, currStackInfo, SET, curr] = SubFcn.getStack_tiffstack(SET, curr, iAni);
                case 'single_tiff'
                    [currStack, currStackInfo, SET, curr] = SubFcn.getStack_singletiff(SET, curr, iAni);
                case 'pst'
                    [currStack, currStackInfo, SET, curr] = SubFcn.getStack_pst(SET, curr, iAni);
                case 'lsm'
                    [currStack, currStackInfo, SET, curr] = SubFcn.getStack_lsm(SET, curr, iAni);
            end% swtich data formats
            % Check for consistency among trials
            if length(unique(curr.duration))>1
                error('Different trials have different durations. We cannot handle this.')
            else
                curr.duration = unique(curr.duration);
            end% if missmatching durations


            % Filter images
            % -------------------------------------------------------------
            currStack = reshape(currStack, [size(currStack,1), size(currStack,2), 1, size(currStack,3)]);
            currStack = num2cell(currStack, [1 2]);
            switch SET.SmoothData_spat_type
                case 'gauss'
                    parfor iT = 1:size(currStack,4)
                        snap = currStack{1,1,1,iT};
                        currStack{1,1,1,iT} = imgaussfilt(snap, SET.SmoothData_spat);
                    end%iT
                case 'med'
                    parfor iT = 1:size(currStack,4)
                        snap = currStack{1,1,1,iT};
                        currStack{1,1,1,iT} = medfilt2(snap, [SET.SmoothData_spat, SET.SmoothData_spat]);
                    end%iT
                case 'med-gauss'
                    parfor iT = 1:size(currStack,4)
                        snap = currStack{1,1,1,iT};
                        snap = medfilt2(snap, [SET.SmoothData_spat(1), SET.SmoothData_spat(1)]);
                        currStack{1,1,1,iT} = imgaussfilt(snap, SET.SmoothData_spat(2));
                    end%iT
                case 'gauss-med'
                    parfor iT = 1:size(currStack,4)
                        snap = currStack{1,1,1,iT};
                        snap = imgaussfilt(snap, SET.SmoothData_spat(2));
                        currStack{1,1,1,iT} = medfilt2(snap, [SET.SmoothData_spat(1), SET.SmoothData_spat(1)]);                        
                    end%iT
            end%switch
            currStack = cell2mat(currStack);
            currStack = reshape(currStack, [size(currStack,1), size(currStack,2), size(currStack,4)]);


            % Detrend each stimulus
            % -------------------------------------------------------------
            if SET.Detrend
                avg_value = median(currStack(:));
                for iStim = 1:curr.duration:size(currStack, 3)                    
                    % Get the current time course
                    Y = squeeze(nanmean(nanmean(currStack(:, :, iStim : iStim+curr.duration-1))));
                    % Fit a line to the data
                    % --- Prepare data
                    [xData, yData] = prepareCurveData(...
                        [1:length(Y)]', Y(:));
                    % ---- Set up fittype and options.
                    ft = fittype( 'poly1' );
                    % ---- Fit model to data.
                    [fitresult, gof] = fit( xData, yData, ft );
                    % Get the (de)trendline
                    detrend = fitresult.p1*xData + fitresult.p2;
                    % Subtract it
                    currStack(:, :, iStim : iStim+curr.duration-1) = currStack(:, :, iStim : iStim+curr.duration-1) - reshape(detrend, [1, 1, curr.duration]);
                    % Subtract glow, i.e. min projection, and add offset
                    Y = currStack(:, :, iStim : iStim+curr.duration-1);
                    currStack(:, :, iStim : iStim+curr.duration-1) = Y - min(Y,[],3) + avg_value;                     
                end%iStim
            end%if detrend


            % Motion correction
            % -------------------------------------------------------------
            % Prepare to keep track of the mask regions
            Mask_region = nan(2,2,SET.N);
            Mask = ones([SET.TargetSize, SET.N]);
            if SET.CorrectMotion
                for iTrial = 1:SET.N
                    idx = find(currStackInfo(:,1) == iTrial);
                    [currStack(:,:,idx), Mask(:,:,iTrial)] = SubFcn.correctMotion(currStack(:,:,idx), currStackInfo(idx,:), SET);
                end
            end% if correct motion


            % Ok, we have corrected for movement. Now it is time to
            % separate data by trial and stimulus
            % Iterate over trials
            for iTrial = 1:SET.N

                % Get info on the valid region
                [row, col] = find(Mask(:,:,iTrial)==1);
                Mask_region(:,:,iTrial) = [min(row), min(col); max(row), max(col)];

                % Create folder
                mkdir([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\Trial',sprintf('%02d', iTrial)])
                % Assign a color to each stimulus
                cols = parula(length(curr.stimuli_seq{iTrial}));

                % Preallocation
                BF_std = zeros(SET.TargetSize);
                BF_avg = zeros(SET.TargetSize);
                cnt_stim = 0;
                % Iterate over stimuli
                for iStim = 1:length(curr.stimuli_seq{iTrial})
                    % Get current data and kick out from stack
                    currData = currStack(:,:,1:curr.duration);
                    if size(currStack, 3) > curr.duration
                        currStack = currStack(:,:,curr.duration+1:end);
                    else
                        clear currStack
                    end

                    % Temporal filtering
                    currData(Mask_region(1,1,iTrial):Mask_region(2,1,iTrial), Mask_region(1,2,iTrial):Mask_region(2,2,iTrial), :) =...
                        imboxfilt3(currData(Mask_region(1,1,iTrial):Mask_region(2,1,iTrial), Mask_region(1,2,iTrial):Mask_region(2,2,iTrial), :), [1 1 SET.SmoothData_temp]);

                    % Get average activity before stimulus onset as the baseline
                    if length(curr.stimulus_onset) == 1
                        BL = curr.stimulus_onset(1)-SET.Baseline:curr.stimulus_onset(1);
                        BL_avg = nanmean(currData(:,:,BL),3);
                        SET.BaselineWindow = BL;
                        % Iterate over each frame, get dF/F values
                        for jFrame = 1:size(currData,3)
                            % Get dF/F
                            currData(:, :, jFrame) = 100 * ((currData(:, :, jFrame) - BL_avg) ./ BL_avg);
                        end%iFrame
                    else
                        SET.BaselineWindow = [];
                        BL_avg = nan([SET.TargetSize, length(curr.stimulus_onset)]);
                        BL_id = ones(1,size(currData,3));
                        % Get BL windows
                        for iBL = 1:length(curr.stimulus_onset)
                            BL = curr.stimulus_onset(iBL)-SET.Baseline:curr.stimulus_onset(iBL);
                            BL_avg(:,:,iBL) = nanmean(currData(:,:,BL),3);
                            SET.BaselineWindow = [SET.BaselineWindow; BL];
                            % Update which BL to use
                            if iBL > 1
                                BL_id(BL(1):end) = BL_id(BL(1):end)+1;
                            end%if next BL
                        end%iBL
                        % Get dF/F
                        for jFrame = 1:size(currData,3)
                            % Get dF/F
                            currData(:, :, jFrame) = 100 * ((currData(:, :, jFrame) - BL_avg(:,:,BL_id(jFrame))) ./ BL_avg(:,:,BL_id(jFrame)));
                        end%iFrame
                    end%if multiple presentations
                    currData(currData==inf) = 0;

                    % Get average response
                    avg_response = squeeze(nanmedian(nanmedian(currData)));

                    % Generate BF as average nanstd or avg  over all
                    % snapshots
                    BF_std = BF_std+squeeze(nanstd(currData,[],3));
                    BF_avg = BF_avg+squeeze(nanmedian(currData,3));

                    % Show raw traces
                    figure(rawTraces);
                    subplot(1,2,1); hold on
                    plot(avg_response, 'Color',cols(iStim,:), 'LineWidth', 2)
                    subplot(1,2,2); hold on
                    if mean(avg_response) > 0
                        rectangle('Position', [iStim-0.25, 0, 0.5 mean(avg_response)], 'FaceColor',cols(iStim,:), 'EdgeColor','k')
                    elseif mean(avg_response) < 0
                        rectangle('Position', [iStim-0.25, nanmean(avg_response), 0.5 abs(mean(avg_response))], 'FaceColor',cols(iStim,:), 'EdgeColor','k')
                    end
                    drawnow

                    % Mean projection
                    avgBF = nanmedian(currData,3);
                    avgBG(~Mask(:,:,iTrial)) = NaN;
                    % Normalized std projection
                    stdBF = nanstd(currData,[],3);                   
                    stdBF = log(sqrt(SubFcn.normalizeImage(stdBF))); 
                    stdBF(~Mask(:,:,iTrial)) = NaN;
                    % --- Now depict everything
                    BFs = figure('color', 'w');
                    set(gcf,'units', 'normalized', 'position', [0 0 1 1])
                    ax = subplot(1,2,1);
                    imagesc(avgBF);
                    colormap(ax, flipud(brewermap(1000,'RdBu')))
                    clim([-nanmax(abs(avgBF(:))), nanmax(abs(avgBF(:)))])
                    axis equal tight off
                    title('avg')
                    clear avgBF
                    ax = subplot(1,2,2);
                    imagesc(stdBF);
                    colormap(ax, 'gray')
                    clim([nanmin(stdBF(:)), nanmax(stdBF(:))])
                    axis equal tight off
                    title('std')
                    clear stdBF
                    export_fig([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\Trial',sprintf('%02d', iTrial),'\',curr.stimuli_seq{iTrial}{iStim},'_BFs'], '-pdf', '-painters')
                    close(BFs)

                    % Put everything together
                    ImageStream.(curr.stimuli_seq{iTrial}{iStim}) = currData; clear currData
                    avgTC.(curr.stimuli_seq{iTrial}{iStim}) = avg_response; clear avg_response

                    % Save everything in the order of stimulus
                    % presentation. Note, this takes ages.
                    save([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\Trial',sprintf('%02d', iTrial),'\',curr.stimuli_seq{iTrial}{iStim},'.mat'], 'ImageStream', 'avgTC', '-v7.3')

                    % Prepare for next round
                    clear ImageStream avgTC BL_avg
                    avg_response = zeros(1,curr.duration);
                end%iStim

                % Export raw traces as pdf
                figure(rawTraces);
                % Avg time course across the whole image
                subplot(1,2,1); hold on
                set(gcf, 'color', 'white')
                box on
                legend(curr.stimuli_seq{iTrial}, 'Location', 'southoutside', 'interpreter', 'none','NumColumns',2)
                ylabel('\DeltaF/F (%)', 'interpreter', 'tex')
                xlabel('time (frames)')
                title([curr.dir.all(iAni).name,'\Trial',sprintf('%02d', iTrial)], 'interpreter', 'none')
                % Avg of the time course
                subplot(1,2,2); hold on
                ylabel('avg. \DeltaF/F (%)', 'interpreter', 'tex')
                title([curr.dir.all(iAni).name,'\Trial',sprintf('%02d', iTrial)], 'interpreter', 'none')
                set(gca, 'XTick', 1:length(curr.stimuli_seq{iTrial}), 'XTickLabels', curr.stimuli_seq{iTrial}, 'TickLabelInterpreter', 'none')
                plot([0 length(curr.stimuli_seq{iTrial})+1],[0 0],'k')
                export_fig([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\Trial',sprintf('%02d', iTrial),'\rawTraces'], '-pdf', '-painters')
                % Depict and save BF images
                BF_std = BF_std/length(curr.stimuli_seq{iTrial});
                BF_avg = BF_avg/length(curr.stimuli_seq{iTrial});
                save([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\Trial',sprintf('%02d', iTrial),'\BF.mat'], 'BF_std','BF_avg','-v7.3')
                % --- nanstd ---
                clf
                % Show normalized BF image (std)
                % --- Normalize
                BF_std = log(sqrt(SubFcn.normalizeImage(BF_std))); 
                % --- show image
                imagesc(BF_std)
                colormap gray
                clim([nanmin(BF_std(:)), nanmax(BF_std(:))])
                axis equal tight off
                title([curr.dir.all(iAni).name,'\Trial',sprintf('%02d', iTrial)], 'interpreter', 'none')
                % Save as pdf
                export_fig([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\Trial',sprintf('%02d', iTrial),'\BF_std'], '-pdf', '-painters')
                % --- AVG ---
                clf
                % Show image                
                imagesc(BF_avg)
                colormap(flipud(brewermap(1001,'RdBu')))
                clim([-nanmax(abs(BF_avg(:))), nanmax(abs(BF_avg(:)))])
                axis equal tight off
                title([curr.dir.all(iAni).name,'\Trial',sprintf('%02d', iTrial)], 'interpreter', 'none')
                % Save as pdf
                export_fig([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\Trial',sprintf('%02d', iTrial),'\BF_avg'], '-pdf', '-painters')
                clear BF_*
                clf
            end%iTrial

            % Export meta information
            meta.animal=curr.dir.all(iAni).name;
            meta.duration=curr.duration;
            meta.stimulus_onset = curr.stimulus_onset;
            meta.baseline = SET.BaselineWindow;
            meta.n_trials=SET.N;
            meta.resolution=SET.TargetSize;
            meta.Mask = Mask;
            meta.Mask_region = Mask_region;
            meta.protocol = curr.protocol;
            meta.SmoothData_spat_type=SET.SmoothData_spat_type;
            meta.SmoothData_spat=SET.SmoothData_spat;
            meta.SmoothData_temp=SET.SmoothData_temp;
            meta.detrend = SET.Detrend;
            meta.stimuli=curr.stimuli_seq;
            meta.valid_trials=1:SET.N;

            save([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\meta_info.mat'], 'meta','-v7.3')

        end%if overwrite
    end%if animal folder
end%iAni
%% Clean and prepare the environment
clc; clear all; close all
disp('----- Ciao, see you next time! -----')






