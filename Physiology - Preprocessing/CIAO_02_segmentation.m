%% Calcium Imaging Analysis Operator (segmentation)
% CIOA_02_segmentation is the second of a series of scripts aiming for an 
% easy-to-use calcium imaging analysis pipeline. Here, we pool trials, 
% subtract a reference stimulus, segment the images into granules, and 
% provide a first assessment of whether a granule was active or not.
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
% Version: 13-Dec-22 (R2022a)

%% Clean and prepare environment
clc; clear all; close all
disp('----- Ciao, welcome back! -----')
% Add paths:
% --- Current folder and subfolders
addpath(genpath(pwd))

%% Settings

% Set reference stimulus, for example the solvent of an odorant. Leave
% empty (SET.ref_stim = []) to skip.
SET.ref_stim = [];

% Set how regional maxima should be determined. Either just based on
% the std-projection ('raw'), a version fitlered with the movmax
% function ('filtered') or a combination of both ('both');
SET.regmax_method = 'filtered';

% Set how many iterations the refinement process during segmentation
% should run
SET.SegmentationRep = 100;

% Set the minimum number of pixels a area should have
SET.minPixCount = 5;

% Set active window used to get the average activity
SET.activeRegion = [55:84];

% Set zScore for a first assessment of which regions responded to the
% stimulus. Note, that this still can be fine-tuned for each animal with
% the script CIAO_02b_customZscore
SET.active_Zscore = 3.5;

% Set whether to run script for all animals, or just for new ones
SET.overwrite = false;

%% Iterate over all animals and their corresponding trials
% Get path to raw data
SET.path2animal = uigetdir(pwd,'Select folder listing all animals');

% Get overview of animal folders
curr.dir.all = dir(SET.path2animal);

% Iterate over all animals
hWait = waitbar(0, {'Please wait...'});
for iAni = 1:size(curr.dir.all,1)
    % Check whether this is what we are looking for, i.e. whether the
    % string "Animal" is present
    if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir

        % Check whether to run preprocessing or not
        if SET.overwrite || ~isfolder([SET.path2animal,'\',curr.dir.all(iAni).name,'\03_Data_Processed'])

            % Update waitbar
            waitbar(iAni/size(curr.dir.all,1), hWait, {['Please wait...   ', num2str(iAni),'/',num2str(size(curr.dir.all,1))]});

            % Get information about trials
            curr.path.data = [SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab'];

            % Get meta info
            load([curr.path.data, '\meta_info.mat']);

            % Prepare to save results
            curr.path.save = [SET.path2animal,'\',curr.dir.all(iAni).name,'\03_Data_Processed'];
            mkdir(curr.path.save)

            % Get a template from the first trial to align data from
            % multiple trials (if there were more than one)
            stimStack = nan([meta.resolution, length(vertcat(meta.stimuli{:}))]);
            stimStack_info = zeros([length(meta.valid_trials)*length(meta.stimuli), 1]);
            meta.unique_stim = unique(vertcat(meta.stimuli{:}));
            cnt = 1;
            for iTrial = 1:length(meta.valid_trials)
                for iStim = 1:length(meta.stimuli{iTrial})
                    if isfile([curr.path.data, '\Trial', sprintf('%02d', meta.valid_trials(iTrial)), '\', meta.stimuli{iTrial}{iStim}, '.mat'])
                        temp = load([curr.path.data, '\Trial', sprintf('%02d', meta.valid_trials(iTrial)), '\', meta.stimuli{iTrial}{iStim}, '.mat']);
                        stimStack(:,:,cnt) = nanstd(temp.ImageStream.(meta.stimuli{iTrial}{iStim}), [], 3);
                        stimStack_info(cnt) = meta.valid_trials(iTrial);                        
                    end
                    cnt = cnt+1;
                end%iStim
            end%iTrial
            % Clear NaN
            stimStack(isnan(stimStack)) = 0;
            % Get template
            idx_template = find(stimStack_info == meta.valid_trials(1));
            fft_template = fft2(nanmax(stimStack(:,:,idx_template),[],3));
            clear temp idx

            % Get overall average and prepare subtraction of reference
            curr.refResponse = zeros([meta.resolution, meta.duration]);
            if ~isempty(SET.ref_stim)
                for iTrial = 1:length(meta.valid_trials)
                    % Load reference
                    ref = load([curr.path.data, '\Trial', sprintf('%02d', meta.valid_trials(iTrial)), '\', SET.ref_stim, '.mat']);
                    ref = ref.ImageStream.(SET.ref_stim);
                    % Register with template
                    idx = find(stimStack_info == meta.valid_trials(iTrial));
                    fft_frame = fft2(nanmax(stimStack(:,:,idx),[],3));
                    xyshifts = SubFcn.dftregister(fft_template, fft_frame, []);
                    ref = SubFcn.shiftframe(ref, xyshifts(1), xyshifts(2));
                    % Pool
                    curr.refResponse(:,:,1:size(ref,3)) = curr.refResponse(:,:,1:size(ref,3)) + ref;
                end%iTrial
            end%if subtract ref
            % Average
            curr.refResponse = curr.refResponse/length(meta.valid_trials);
            clear fft_frame xyshifts ref iTrial

            % Now, load each stimulus, form the average over all trials and
            % subtract the response to the reference.
            % Preallocation
            poolEdge = zeros(meta.resolution);
            % Iterate over each stimulus            
            for iStim = 1:length(meta.unique_stim)
                % Preallocation
                clear ImageStream
                ImageStream.(meta.unique_stim{iStim}) = zeros(size(curr.refResponse));
                cnt_stim = 0;
                % Iterate over all valid trials
                for iTrial = 1:length(meta.valid_trials)
                    % Load data
                    path2file = [curr.path.data, '\Trial', sprintf('%02d', meta.valid_trials(iTrial)), '\', meta.unique_stim{iStim}, '.mat'];
                    if isfile(path2file)
                        temp_img = load([curr.path.data, '\Trial', sprintf('%02d', meta.valid_trials(iTrial)), '\', meta.unique_stim{iStim}, '.mat']);
                        temp_img = temp_img.ImageStream.(meta.unique_stim{iStim});
                        % Register with template
                        idx = find(stimStack_info == meta.valid_trials(iTrial));
                        fft_frame = fft2(nanmax(stimStack(:,:,idx), [], 3));
                        xyshifts = SubFcn.dftregister(fft_template, fft_frame, []);
                        [temp_img, edge] = SubFcn.shiftframe(temp_img, xyshifts(1), xyshifts(2));
                        poolEdge = poolEdge+edge;
                        % Pool data
                        ImageStream.(meta.unique_stim{iStim}) = ImageStream.(meta.unique_stim{iStim}) + temp_img;
                        % Update counter
                        cnt_stim = cnt_stim+1;
                    end
                    clear temp_img fft_frame xyshifts
                end%iTrial
                % Form average and subtract the reference stimulus
                ImageStream.(meta.unique_stim{iStim}) = ImageStream.(meta.unique_stim{iStim})/cnt_stim;
                ImageStream.(meta.unique_stim{iStim}) = ImageStream.(meta.unique_stim{iStim}) - curr.refResponse;
                % Save all
                save([curr.path.save,'\',meta.unique_stim{iStim},'.mat'], 'ImageStream', '-v7.3')
                clear ImageStream responseActivity baselineActivity
            end%iStim


            % Prepare to set the edge to nan
            poolEdge = poolEdge>0;
            Mask = ones(meta.resolution);
            Mask(find(poolEdge)) = NaN;
            meta.Mask = mean(cat(3,Mask, meta.Mask),3);
            meta.Mask = double(meta.Mask==1);
            [row, col] = find(meta.Mask==1);
            if isempty(row) || isempty(col)
                meta.Mask_region = zeros(2);
            else
                % Crop mask region by one to be on the save side
                meta.Mask_region = [min(row)+1, min(col)+1; max(row)-1, max(col)-1];
                meta.Mask = zeros(size(meta.Mask));
                meta.Mask(meta.Mask_region(1,1):meta.Mask_region(2,1), meta.Mask_region(1,2):meta.Mask_region(2,2)) = 1;
            end

            %% 
            if meta.baseline(end) == 20
                SET.activeRegion = 21:27;
            elseif meta.baseline(end) == 10
                SET.activeRegion = 11:17;
            end
            %%

            % Update and save meta info
            meta.activeRegion = SET.activeRegion;
            meta.active_Zscore = SET.active_Zscore;
            meta.ref_stim = SET.ref_stim;
            meta.segmentation.regmax_method = SET.regmax_method;
            meta.segmentation.SegmentationRep = SET.SegmentationRep;
            clear row col poolEdge
            % Save new meta info
            save([curr.path.save, '\meta_info.mat'], 'meta');

            % Apply mask. For this, load each stimulus again. At the same
            % time, use this to get a stack pooling all stimuli.
            PoolStream = [];
            PoolStream_info = [];
            for iStim = 1:length(meta.unique_stim)
                % Get data
                load([curr.path.save,'\',meta.unique_stim{iStim},'.mat'])
                % Multiply with mask
                ImageStream.(meta.unique_stim{iStim}) = ImageStream.(meta.unique_stim{iStim}).*meta.Mask;
                % Pool
                PoolStream = cat(3, PoolStream, ImageStream.(meta.unique_stim{iStim}));
                PoolStream_info = [PoolStream_info; ones(size(ImageStream.(meta.unique_stim{iStim}), 3), 1) * iStim];
            end%iStim

            % Segment AL by defining activity 'pockets'. For this, identify
            % local maxima in the std-projection and
            SET.Zscore_props.active_Zscore = SET.active_Zscore;
            SET.Zscore_props.Baseline = meta.baseline;
            SET.Zscore_props.Window = SET.activeRegion;

            [Segmentation.pockets_bw, Segmentation.pockets, Segmentation.pockets_labeled, Segmentation.pocket_Zscore, Segmentation.pocket_Zscore_img, Segmentation.pocket_Corr, Segmentation.pocket_Corr_img] =...
                SubFcn.stackSegmentation(PoolStream, PoolStream_info, meta.Mask, SET.regmax_method, SET.SegmentationRep, SET.minPixCount, SET.Zscore_props);
            
            % Save this
            save([curr.path.save, '\img_segmentation.mat'], 'Segmentation', 'SET');

        end%if overwrite
    end%if animal folder
end%iAni

%% Clean and prepare environment
close(hWait)
clc; clear all; close all
disp('----- Ciao, see you next time! -----')