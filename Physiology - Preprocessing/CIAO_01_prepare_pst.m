%% Calcium Imaging Analysis Operator (prepare)
% CIOA_01_prepare is the first of a series of scripts aiming for an
% easy-to-use calcium imaging analysis pipeline. Here, we collect and
% prepare raw, apply spatial and temporal filters, and correct baseline
% shifts due to, for example, bleaching.
% CIOA assumes a strict folder structure, and each analysis step is saved in
% a given folder. Thus, each module can be run independently, and results
% from different stages can be used without running a long pipeline.
%
% In the following, we will lay out the folder structure needed for CIAO to
% run without errors:
%   > DATA                       : this is the folder you select at the beginning
%   |--> Animal01                : this folder can have any name that contains the string "Animal01"
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
%
% Version: 03-Mar-23 (R2022a)
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
%% Settings
% Set nanstd for 2D Gaussian filtering ('gauss'), or box size for
% 2D median filtering (med).
% Note that spatial filtering is applied before resizing the image.
SET.SmoothData_spat_type = 'med';
SET.SmoothData_spat = 3;
% Set whether to correct for motion
SET.CorrectMotion = true;
% Set the maximum number of iterations for image registration
SET.maxRegIter = 10;
% Set how much of the edge should be cut off (pixels)
SET.edge = 2;
% Set depth for 3D box filtering with no spatial filter
SET.SmoothData_temp = 3;
% Set baseline length (frames). The window will end at the stimulus onset
SET.Baseline = 9;
% Set whether to detrend data
SET.Detrend = true;
% Set whether to run the script for all animals or just for new ones
SET.overwrite = true;
%% Iterate over all animals and their corresponding trials
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
            currStack = [];
            currStackInfo = [];
            for iTrial = 1:SET.N
                % Keep path
                curr.path.trial = [SET.path2animal,'\',curr.dir.all(iAni).name,'\01_Data_raw','\Trial',sprintf('%02d', iTrial)];
                % Get meta inforamtion
                curr.protocol{iTrial} = readtable([curr.path.trial,'\protocol.csv']);
                % Identify sequence of stimuli
                curr.stimuli_seq{iTrial} = curr.protocol{iTrial}.Odour;
                % Get duration of each stimulus [frames]
                curr.duration = curr.protocol{iTrial}.num_frames(1);
                % Get the dimensions of the frames
                SET.TargetSize = [curr.protocol{iTrial}.img_width(1), curr.protocol{iTrial}.img_height(1)];
                % Preallocation
                poolEdge = zeros([SET.TargetSize, SET.N]);
                % Iterate over all stacks across all trials and create a
                % large stack
                for iStim = 1:length(curr.stimuli_seq{iTrial})
                    % Get stacks of current stimulus
                    fr340 = SubFcn.pstread(SET.TargetSize(1), SET.TargetSize(2), curr.duration, [curr.path.trial, '\', curr.protocol{iTrial}.DBB1{iStim}, '.pst']);
                    fr380 = SubFcn.pstread(SET.TargetSize(1), SET.TargetSize(2), curr.duration, [curr.path.trial, '\', curr.protocol{iTrial}.dbb2{iStim}, '.pst']);
                    % Iterate over frames, smooth and form the ratio
                    for iFrame = 1:curr.duration
                        % SMooth before forming the ratio
                        switch SET.SmoothData_spat_type
                            case 'gauss'
                                snap_340 = imgaussfilt(fr340(:,:,iFrame), SET.SmoothData_spat);
                                snap_380 = imgaussfilt(fr380(:,:,iFrame), SET.SmoothData_spat);
                            case 'med'
                                snap_340 = imgaussfilt(fr340(:,:,iFrame), [SET.SmoothData_spat, SET.SmoothData_spat]);
                                snap_380 = imgaussfilt(fr380(:,:,iFrame), [SET.SmoothData_spat, SET.SmoothData_spat]);
                        end
                        % Ratio
                        snap = snap_340 ./ snap_380;
                        snap = snap*100;
                        % Stack
                        currStack = cat(3, currStack, snap);
                        currStackInfo = [currStackInfo; iTrial, iStim, iFrame];
                        cnt = cnt+1;
                    end%iFrame
                end%iStim
            end%iTrial
            clear fr340 fr380 snap_340 snap_380
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
                    % Correct frameshifts by translation
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
            end
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
                    % Get avg image during baseline
                    curr.stimulus_onset = curr.protocol{iTrial}.StimON(iStim);
                    SET.BaselineWindow = curr.stimulus_onset-SET.Baseline:curr.stimulus_onset;
                    BL_avg = nanmean(currData(:,:,SET.BaselineWindow),3);
                    % Subtract baseline
                    currData = currData-BL_avg;
                    % Get average response
                    avg_response = squeeze(nanmean(nanmean(currData)));
                    % Generate BF as average nanstd or avg  over all
                    % snapshots
                    BF_std = BF_std + squeeze(nanstd(currData,[],3));
                    BF_avg = BF_avg + squeeze(nanmean(currData,3));
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
                    % Put snaps together with baseline already subtracted
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
                % Show image
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