%% Calcium Imaging Analysis Operator (prepare)
% CIOA_01_prepare is the first of a series of scripts aiming for an
% easy-to-use calcium imaging analysis pipeline. Here, we collect and
% prepare raw, apply spatial and temporal filters, and correct baseline
% shifts due to, for example, bleaching.
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
%   |  |  |--> Trial##           : same logic for all following trials
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
% Version: 03-Mar-23 (R2023a)

%% Clean and prepare environment
clc; close all
disp('----- Ciao, welcome back! -----')
% Add paths:
% --- Current folder and subfolders
addpath(genpath(pwd))

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig)
addpath(genpath('...\GitHub\export_fig'))

% Matlab routines for online non-rigid motion correction of
% calcium imaging data
% (https://github.com/flatironinstitute/NoRMCorre)
addpath(genpath('...\GitHub\NoRMCorre'))

%% Settings
% Set basename for individual frames
SET.basename = 'ss_single_';
SET.filetype = '.tiff';

% Set nanstd for 2D Gaussian filtering ('gauss'), or box size for
% 2D median filtering (med).
% Note that spatial filtering is applied before resizing the image.
SET.SmoothData_spat_type = 'med';
SET.SmoothData_spat = 11;

% Set target size of images to decreasefilesize.
SET.TargetSize = [256 256];

% Set how much of the edge should be cut off (pixels)
SET.edge = 5;

% Set whether to correct for motion
SET.CorrectMotion = true;

% Set the maximum number of iterations for image registration
SET.maxRegIter = 10;

% Set depth for 3D box filtering with no spatial filter
SET.SmoothData_temp = 5;

% Set baseline length (frames). The window will end at stimulus onset
SET.Baseline = 20;

% Set whether to detrend data. For a course trend correction, we subtract a
% trendline based on an linear fit to each stimulus' raw pixel intensity
% data.
SET.Detrend = true;

% Set whether to run script for all animals, or just for new ones
SET.overwrite = true;

%% Iterate over all animals and their corresponding trials
% Get path to raw data
SET.path2animal = uigetdir(pwd,'Select folder listing all animals');

% Get overview of animal folders
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
            % Preallocation
            currStack = [];
            poolEdge = zeros([SET.TargetSize, SET.N]);
            currStackInfo = [];
            curr.stimuli_seq = cell(SET.N,1);
            % Iterate over trials
            for iTrial = 1:SET.N
                % Keep path
                curr.path.trial = [SET.path2animal,'\',curr.dir.all(iAni).name,'\01_Data_raw','\Trial',sprintf('%02d', iTrial)];

                % Get meta inforamtion
                opts = delimitedTextImportOptions("NumVariables", 2);
                opts.Delimiter = ":";
                opts.VariableNames = ["Parameter", "Value"];
                opts.VariableTypes = ["string", "string"];
                curr.protocol{iTrial} = readtable([curr.path.trial,'\protocol.txt'], opts);

                % Get onset of stimuli
                % -- Check whether the same stimulus has been given
                %    multiple times throughout the trial
                nPresentation = sum(contains(curr.protocol{iTrial}.Parameter, '_on'));
                curr.stimulus_onset = nan(nPresentation,1);
                for iPresentation = 1:nPresentation
                    curr.stimulus_onset(iPresentation) = str2double(curr.protocol{iTrial}.Value(strcmp(curr.protocol{iTrial}.Parameter,['S',num2str(iPresentation),'_on'])));
                end%iPresentation
                clear nPresentation

                % Itentify sequence of stimuli
                curr.stimuli_seq{iTrial} = split(curr.protocol{iTrial}.Value(find(strcmp(curr.protocol{iTrial}.Parameter, 'Stimuli'))), ', ');

                % Get duration of each stimulus [frames]
                curr.duration(iTrial) = str2double(curr.protocol{iTrial}.Value(find(strcmp(curr.protocol{iTrial}.Parameter, 'Duration'))));

                % Identify prefix
                SET.prefix = [];
                prefix = 1;
                while isempty(SET.prefix)
                    try
                        snap = imread([curr.path.trial,'\',SET.basename,sprintf(['%0',num2str(prefix),'d'],1), SET.filetype]);
                        SET.prefix = num2str(prefix);
                        SET.originalSize = size(snap);
                    catch
                        prefix = prefix+1;
                    end%try
                    if prefix>100
                        error('Unable to load frames.')
                    end%if
                end%while
                clear snap prefix

                % Prealloccation
                Stack = nan([SET.TargetSize, (length(curr.stimuli_seq{iTrial})*curr.duration(iTrial))]);
                % Iterate over all frames across all trials and create a
                % stack
                for iFrame = 1:(length(curr.stimuli_seq{iTrial})*curr.duration)
                    % Get current snapshot
                    snap = imread([curr.path.trial,'\',SET.basename,sprintf(['%0',SET.prefix,'d'],iFrame),SET.filetype]);
                    snap = im2double(snap);
                    % Spatial filtering
                    switch SET.SmoothData_spat_type
                        case 'gauss'
                            snap = imgaussfilt(snap, SET.SmoothData_spat);
                        case 'med'
                            snap = medfilt2(snap, [SET.SmoothData_spat SET.SmoothData_spat]);
                    end
                    % Resize snap
                    snap = imresize(snap, SET.TargetSize, 'box');
                    % Stack
                    Stack(:,:,iFrame) = snap;
                    currStackInfo = [currStackInfo; iTrial];
                end%iFrame
                % Stack trials togehter
                currStack = cat(3, currStack, Stack);
                clear Stack
            end%iTrial

            % Check for consistency among trials
            if length(unique(curr.duration))>1
                error('Different trials have different durations. We cannot handle this.')
            else
                curr.duration = unique(curr.duration);
            end

            % Motion correction
            % -------------------------------------------------------------
            % Prepare to keep track of the mask regions
            Mask_region = nan(2,2,SET.N);
            Mask = ones([SET.TargetSize, SET.N]);
            if SET.CorrectMotion
                % First, align frames within stimuli using a non-rigid motion
                % correction
                for iStim = 1:curr.duration:size(currStack, 3)
                    Y = currStack(:, :, iStim : iStim+curr.duration-1);
                    options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[size(Y,1)/8,size(Y,2)/8],'mot_uf',4,'bin_width',20,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',curr.duration/2,'print_msg',false,'use_parallel','true','iter',2);
                    Y = normcorre_batch(Y, options_nonrigid);
                    currStack(:, :, iStim : iStim+curr.duration-1) = Y;
                    clear Y
                end% iStim

                % Second, shift all frames that belong together in order to
                % match different stimuli from the same trial
                % The registration algorithm works in the Fourier domain (look
                % up the convolution theorem if you want to know why), so we
                % start by computing the 2D Fourier transform of the raw images
                % and registration template
                cnt_while = 1;
                while true
                    % Get std-projection per stimulus for aligning individual
                    % stimuli
                    % --- Preallocation
                    stimStack_std = [];
                    trial_list = zeros(length(vertcat(curr.stimuli_seq{:})),1);
                    cnt = 1;
                    cnt_trial = 1;
                    for iStim = 1:curr.duration:size(currStack, 3)
                        stimStack_std = cat(3, stimStack_std, nanstd(currStack(:, :, iStim : iStim+curr.duration-1), [], 3));
                        trial_list(cnt) = currStackInfo(iStim,1);
                        cnt = cnt+1;
                    end
                    % Determine the needed shift to match the first stimulus as
                    % template-stimulus
                    % --- Preallocation of x/y shift
                    xyshifts = zeros(2, size(stimStack_std, 3));
                    % Iterate over stimuli
                    for iStim = 1:size(stimStack_std, 3)
                        % Get fft2 from template
                        idx1 = find(trial_list==trial_list(iStim));
                        [~, idx2] = min(abs(squeeze(mean(mean(stimStack_std(:,:,idx1)))) - median(squeeze(mean(mean(stimStack_std(:,:,idx1)))))));
                        fft_template_avg = fft2(stimStack_std(:,:,idx1(idx2)));
                        % Get fft2 from current stimulus
                        fft_frame = fft2(stimStack_std(:,:,iStim));
                        % Determine needed translation
                        xyshifts(:, iStim) = SubFcn.dftregister(fft_template_avg, fft_frame, []);
                    end%iFrame
                    % Correct frame shifts by translation
                    cnt=1;
                    for iStim = 1:curr.duration:size(currStack, 3)
                        [currStack(:,:,iStim : iStim+curr.duration-1), edge] = ...
                            SubFcn.shiftframe(...
                            currStack(:,:,iStim : iStim+curr.duration-1), ...
                            xyshifts(1,cnt),...
                            xyshifts(2,cnt));
                        poolEdge(:,:,trial_list(cnt)) = poolEdge(:,:,trial_list(cnt))+edge;
                        cnt = cnt+1;
                    end%iStim
                    % Stop when all done
                    if sum(abs(xyshifts(:))) == 0 || cnt_while>SET.maxRegIter
                        break
                    end% if done
                    cnt_while = cnt_while+1;
                end%while
                % Set the edge to nan
                poolEdge = poolEdge>0;
                Mask(find(poolEdge)) = NaN;
            end%if correct motion


            % Ok, we have corrected for movement. Now it is time to separate
            % data by trial and stimulus
            % Iterate over trials
            for iTrial = 1:SET.N

                % Get info on the valid region
                [row, col] = find(Mask(:,:,iTrial)==1);
                Mask_region(:,:,iTrial) = [min(row)+SET.edge, min(col)+SET.edge; max(row)-SET.edge, max(col)-SET.edge];
                Mask(:,:,iTrial) = NaN;
                Mask(Mask_region(1,1,iTrial):Mask_region(2,1,iTrial), Mask_region(1,2,iTrial):Mask_region(2,2,iTrial), iTrial) = 1;

                % Create folder
                mkdir([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\Trial',sprintf('%02d', iTrial)])

                % Assign a color to each stimulus
                cols = parula(length(curr.stimuli_seq{iTrial}));

                % Preallocation
                BF_std = zeros(SET.TargetSize);
                BF_avg = zeros(SET.TargetSize);

                % Iterate over stimuli
                for iStim = 1:length(curr.stimuli_seq{iTrial})

                    % Get current data and kick out from stack
                    currData = currStack(:,:,1:curr.duration);
                    if size(currStack, 3) > curr.duration
                        currStack = currStack(:,:,curr.duration+1:end);
                    else
                        clear currStack
                    end

                    % Apply mask
                    currData = currData.*Mask(:,:,iTrial);

                    % If this is the first stimulus, get the median pixel
                    % intensity as baseline
                    if iStim == 1
                        baseline_pixel = nanmedian(currData(:));
                    end

                    % Subtract 'glow' which is the min-projection across
                    % all frames
                    currData = currData-nanmin(currData,[],3);

                    % Detrend signal. For this, fit a line to the pixel
                    % intensity data. This will reuslt in an offset.
                    % Correct this by using the value stored in the
                    % variable baseline_pixel.
                    if SET.Detrend
                        % Fit a line to the data
                        % --- Prepare data
                        [xData, yData] = prepareCurveData(...
                            [1:size(currData,3)]',...
                            squeeze(nanmedian(nanmedian(currData(Mask_region(1,1,iTrial):Mask_region(2,1,iTrial), Mask_region(1,2,iTrial):Mask_region(2,2,iTrial),:)))) );
                        % ---- Set up fittype and options.
                        ft = fittype( 'poly1' );
                        % ---- Fit model to data.
                        [fitresult, gof] = fit( xData, yData, ft );
                        % Get the (de)trendline
                        detrend = fitresult.p1*xData + fitresult.p2;
                        detrend = reshape(detrend, [1,1,size(currData,3)]);
                        % Correct trend
                        currData = currData - detrend;
                        % Correct resulting offset
                        currData = currData - nanmedian(currData(:)) + baseline_pixel;
                        % Clean up
                        clear detrend fitresult gof ft xData yData
                    end

                    % Temporal filtering
                    currData(Mask_region(1,1,iTrial):Mask_region(2,1,iTrial), Mask_region(1,2,iTrial):Mask_region(2,2,iTrial), :) =...
                        imboxfilt3(currData(Mask_region(1,1,iTrial):Mask_region(2,1,iTrial), Mask_region(1,2,iTrial):Mask_region(2,2,iTrial), :), [1 1 SET.SmoothData_temp]);

                    % Get average activty before stimulus onset as baseline as baseline
                    BL = curr.stimulus_onset(1)-SET.Baseline:curr.stimulus_onset(1);
                    BL_avg = mean(currData(:,:,BL),3);
                    SET.BaselineWindow = BL;

                    % Iterate over each frame, get dF/F values
                    for jFrame = 1:size(currData,3)
                        % Get dF/F
                        currData(:, :, jFrame) = 100 * ((currData(:, :, jFrame) - BL_avg) ./ BL_avg);
                    end%iFrame
                    currData(currData==inf) = 0;

                    % Get average response
                    avg_response = squeeze(nanmean(nanmean(currData)));

                    % Generate BF as average nanstd or avg  over all
                    % snapshots
                    BF_std = BF_std+squeeze(nanstd(currData,[],3));
                    BF_avg = BF_avg+squeeze(nanmean(currData,3));

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

                    % Show BFs
                    avgBF = nanmean(currData,3);
                    % --- Normalize
                    avgBF = avgBF - nanmin(avgBF(:));
                    avgBF = sqrt(avgBF);
                    avgBF = avgBF - nanmean(avgBF(:));
                    avgBF = avgBF / nanstd(avgBF(:));
                    % --- Now the std projection
                    stdBF = nanstd(currData,[],3);
                    % --- Normalize
                    stdBF = stdBF - nanmin(stdBF(:));
                    stdBF = sqrt(stdBF);
                    stdBF = stdBF - nanmean(stdBF(:));
                    stdBF = stdBF / nanstd(stdBF(:));
                    % --- Now depict everything
                    BFs = figure('color', 'w');
                    set(gcf,'units', 'normalized', 'position', [0 0 1 1])
                    subplot(1,2,1)
                    imagesc(avgBF); clear avgBF
                    colormap(SubFcn.ColMapInferno)
                    axis equal
                    axis off
                    title('avg')
                    subplot(1,2,2)
                    imagesc(stdBF); clear stdBF
                    colormap(SubFcn.ColMapInferno)
                    axis equal
                    axis off
                    title('std')
                    export_fig([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\Trial',sprintf('%02d', iTrial),'\',curr.stimuli_seq{iTrial}{iStim},'_BFs'], '-pdf', '-painters')
                    close(BFs)

                    % Put severything together
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
                % Avg time course across whole image
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
                BF_std = BF_std - nanmin(BF_std(:));
                BF_std = sqrt(BF_std);
                BF_std = BF_std - nanmean(BF_std(:));
                BF_std = BF_std / nanstd(BF_std(:));
                % --- show image
                imagesc(BF_std)
                colormap(SubFcn.ColMapInferno)
                axis equal
                axis off
                title([curr.dir.all(iAni).name,'\Trial',sprintf('%02d', iTrial)], 'interpreter', 'none')
                % Save as pdf
                export_fig([SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\Trial',sprintf('%02d', iTrial),'\BF_std'], '-pdf', '-painters')
                % --- AVG ---
                clf
                % Show image
                % --- Normalize
                BF_avg = BF_avg - nanmin(BF_avg(:));
                BF_avg = sqrt(BF_avg);
                BF_avg = BF_avg - nanmean(BF_avg(:));
                BF_avg = BF_avg / nanstd(BF_avg(:));
                % --- show image
                imagesc(BF_avg)
                colormap(SubFcn.ColMapInferno)
                axis equal
                axis off
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
            meta.edge = SET.edge;
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

%% Clean and prepare environment
clc; clear all; close all
disp('----- Ciao, see you next time! -----')
