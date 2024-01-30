
classdef SubFcn

    properties
    end

    methods(Static)

        function xyShifts = dftregister(fft_template, fft_frames, maxShift)
            % Register frames using Kuglin and Hines phase correlation method
            %
            % function xyShifts = dftregister(fft_template, fft_frames, maxShift)
            %
            % Purpose
            % Calculate the x/y translation shift between one or more frames and a template frame.
            %
            % Inputs
            % fft_template - The template (fixed) 2-D image to which we will register fft_frames.
            % fft_frames - One more more (moving) frames that will be registered to the template
            % maxShift - Shifts larger than "maxShift" pixels are not allowed. If empty, no
            %            maximum shift constraint is applied.
            %
            % Outputs
            % xyShifts - An array of shifts. First row are x shifts and second row y shifts.
            %
            % Adapted from Maxime Rio, 2017

            % weighting coefficient to balance between phase correlation (alpha = 1)
            % and normal correlation (alpha = 0, no normalization)
            alpha = 0.5;

            % compute phase correlation from the normalized cross-spectrum
            cs = fft_template .* conj(fft_frames);
            cc = ifft2(cs ./ abs(cs).^alpha, 'symmetric');

            % constrain maximum shifts found
            if ~isempty(maxShift)
                cc(:, maxShift+2:end-maxShift, :) = NaN;
                cc(maxShift+2:end-maxShift, :) = NaN;
            end

            % split input dimensions
            cc_dims = num2cell(size(cc));
            [nx, ny] = deal(cc_dims{1:2});

            if numel(cc_dims) > 2
                other_dims = cc_dims(3:end);
            else
                other_dims = {1};
            end

            % deduce (x,y)-shift from the maximum in phase correlation
            [~, idx] = max(reshape(cc, [], other_dims{:}), [], 1);
            [x, y] = ind2sub([nx, ny], idx);

            % compensate for circular shifts
            x_mask = x > fix(nx / 2);
            x(x_mask) = x(x_mask) - nx;

            y_mask = y > fix(ny / 2);
            y(y_mask) = y(y_mask) - ny;

            % compensate for 1-based indexing
            x = x - 1;
            y = y - 1;

            % concatenate (x,y)-shifts
            xyShifts = cat(1, x, y);
        end%FCN:dftregister

        function [trFrame, edge] = shiftframe(frame, sx, sy)
            % Fast frame translation
            %
            % function trFrame = shiftframe(frame, sx, sy)
            %
            % Purpose
            % Translate a frame by filling an output frame subregion with the input
            % frame, which is faster than imstranslate or imdilate-based translation,
            % for integer shifts
            %
            % Inputs
            % frame - the frame to shift
            % sx - x shift
            % sy - y shift
            %
            % Outputs
            % trFrame - the translated frame.
            %
            %
            % Adapted from Maxime Rio, 2017
            % preallocated result
            trFrame = nan(size(frame), 'like', frame);
            % fill part of output frame with input frame
            [nx, ny, ~] = size(frame);
            trFrame(max(1, 1+sx):min(nx, nx+sx), max(1, 1+sy):min(ny, ny+sy), :) = ...
                frame(max(1, 1-sx):min(nx-sx, nx), max(1, 1-sy):min(ny-sy, ny), :);
            % Keep track of the edge
            edge = isnan(nanmean(trFrame, 3));
            trFrame(isnan(trFrame)) = 0;
        end%FCN:shiftframes

        function [Data, varargout] = lsmread(fileName, varargin)
            %LSMREAD returns information and image data of a ZEISS lsm file.
            % [Data] = lsmread(filename) will return image data in its original
            % dimensions, organized as [T(ime), C(hannel), Z, X, Y]
            % [Data] = lsmread(filename,'Range',[T0 T1 C0 C1 Z0 Z1]) will return image data with specified dimension.
            % [Data LSMinfo] = lsmread(filename) will return image data as well as LSMinfo
            % [LSMinfo] = lsmread(filename,'InfoOnly') will only return LSM infos.
            %---------------------------------------------------------------------------------
            % Notes for the Zeiss LSM format.
            % The first 8 bytes are TIFF file header.
            % Each 12-byte IFD entry has the following format: Bytes 0-1: The Tag that
            % identifies the field. Byte 2-3: The field type. Bytes 4-7: Count.
            % Bytes 8-11: Value or offset.Byte count of the indicated types.
            % Types 1 = BYTE, 2 =ASCII (8-bit),3= SHORT (16-bit), 4= LONG(32-bit).
            % If the count is 1, then byes 8-11 stores the value, otherwise it's offset
            % pointing to the position in file where value is stored.
            % Each optical section at every time point (no matter how many channels)
            % will have two image file directories(IFDs). The first one contains
            % information for the real image data. The second IFD contains
            % information for thumbnail images. The LSM file is structured as such :
            % header -> IFDs -> (real image/thumbnail)s.
            % The very first IFD has an entry numbered 34412 and points to LSM-specific data.
            %----------------------------------------------------------------------------------
            % Version 1.7
            % Copyright Chao-Yuan Yeh, 2016.
            % The script is inspired by LSM File Toolbox by Peter Li and tiffrd by Francois Nedelec.

            fID = fopen(fileName);
            byteOrder = fread(fID,2,'*char')';
            if (strcmp(byteOrder,'II'))
                byteOrder = 'ieee-le';
            elseif (strcmp(byteOrder,'MM'))
                byteOrder = 'ieee-be';
            else
                error('This is not a correct TIFF file');
            end

            tiffID = fread(fID, 1, 'uint16', byteOrder);
            if (tiffID ~= 42)
                error('This is not a correct TIFF file');
            end

            fseek(fID, 4, 'bof');
            ifdPos = fread(fID, 1, 'uint32', byteOrder);
            fseek(fID, ifdPos, 'bof');

            IFDIdx = 0;
            while ifdPos ~= 0
                IFDIdx = IFDIdx+1;
                fseek(fID,ifdPos, 'bof');
                % The first two bytes of each IFD specify number of IFD entries.
                numEntries = fread(fID,1,'uint16',byteOrder);
                entryPos = ifdPos+2;
                % Each IFD entry is 12-byte long.
                fseek(fID, ifdPos+12*numEntries+2, 'bof');
                % The last four bytes of an IFD specifies offset to next IFD.
                % If this is zero, it means there's no other IFDs.
                ifdPos = fread(fID, 1, 'uint32', byteOrder);
                % IFD is structured like this: bytes 1-2 : tag, bytes 3-4: type,
                % bytes 5-8: count, bytes 9-12: value/offset
                for ii = 1:numEntries
                    fseek(fID, entryPos+12*(ii-1), 'bof');
                    IFD{IFDIdx}(1,ii) = fread(fID, 1, 'uint16', byteOrder);
                    IFD{IFDIdx}(2,ii) = fread(fID, 1, 'uint16', byteOrder);
                    IFD{IFDIdx}(3,ii) = fread(fID, 1, 'uint32', byteOrder);
                    IFD{IFDIdx}(4,ii) = fread(fID, 1, 'uint32', byteOrder);
                end
            end

            %Reading LSMinfo
            if IFD{1}(3, IFD{1}(1, :) == 258) == 1
                LSMinfo.bitDepth = IFD{1}(4, IFD{1}(1, :) == 258);
            else
                fseek(fID, IFD{1}(4,IFD{1}(1, :) == 258),'bof');
                LSMinfo.bitDepth = fread(fID, 1, 'uint16', byteOrder);
            end

            offsetLSMinfo = IFD{1}(4,IFD{1}(1, :) == 34412)+8;
            fseek(fID, offsetLSMinfo, 'bof');
            LSMinfo.dimX = fread(fID, 1, 'uint32', byteOrder);
            LSMinfo.dimY = fread(fID, 1, 'uint32', byteOrder);
            LSMinfo.dimZ = fread(fID, 1, 'uint32', byteOrder);
            LSMinfo.dimC = fread(fID, 1, 'uint32', byteOrder);
            LSMinfo.dimT = fread(fID, 1, 'uint32', byteOrder);
            fseek(fID, 12, 'cof');
            LSMinfo.voxSizeX = fread(fID, 1, 'float64', byteOrder);
            LSMinfo.voxSizeY = fread(fID, 1, 'float64', byteOrder);
            LSMinfo.voxSizeZ = fread(fID, 1, 'float64', byteOrder);
            fseek(fID, 26, 'cof');
            LSMinfo.specScan = fread(fID, 1, 'uint16', byteOrder);

            if any(strcmpi(varargin, 'InfoOnly'))
                Data = LSMinfo;
                fclose(fID);
            else
                % Creating offset list for every IFD
                offset = zeros(1, IFDIdx/2);
                for ii= 1 : IFDIdx/2
                    if IFD{2*ii - 1}(3,7) == 1 && IFD{2 * ii - 1}(4, 7) < 4294967296
                        offset(ii) = IFD{2 * ii - 1}(4, 7);
                    else
                        fseek(fID, IFD{2 * ii - 1}(4, 7), 'bof');
                        offset(ii) = fread(fID, 1, 'uint32', byteOrder);
                    end
                end

                if any(strcmpi(varargin, 'Range'))
                    Range = varargin{find(strcmpi(varargin, 'Range')) + 1};
                    if any([Range(1) < 1, Range(2) > LSMinfo.dimT, Range(3) < 1, ...
                            Range(4) > LSMinfo.dimC, Range(5) < 1, ...
                            Range(6)>LSMinfo.dimZ])
                        error('Input range exceeds data range');
                    else
                        dimT0 = Range(1);
                        dimT1 = Range(2);
                        dimC0 = Range(3);
                        dimC1 = Range(4);
                        dimZ0 = Range(5);
                        dimZ1 = Range(6);
                    end
                else
                    dimT0 = 1; dimT1 = LSMinfo.dimT;
                    dimZ0 = 1; dimZ1 = LSMinfo.dimZ;
                    dimC0 = 1; dimC1 = LSMinfo.dimC;
                end

                bitDepth = strcat('uint',num2str(LSMinfo.bitDepth));
                Data = zeros(dimT1-dimT0+1, dimC1-dimC0+1, dimZ1-dimZ0+1,...
                    LSMinfo.dimY, LSMinfo.dimX, bitDepth);

                for indT = dimT0 : dimT1
                    for indZ = dimZ0 : dimZ1
                        fseek(fID, offset((indT-1) * LSMinfo.dimZ + indZ),'bof');
                        switch bitDepth
                            case 'uint16'
                                fseek(fID, (dimC0-1)*LSMinfo.dimX*LSMinfo.dimY*2, 'cof');
                                for indC = dimC0 : dimC1
                                    Data(indT-dimT0+1, indC-dimC0+1, indZ-dimZ0+1, :, :) ...
                                        = reshape(uint16(...
                                        fread(fID, LSMinfo.dimX * LSMinfo.dimY, ...
                                        bitDepth, byteOrder)), LSMinfo.dimX, LSMinfo.dimY)';
                                end
                            case 'uint8'
                                fseek(fID, (dimC0-1)*LSMinfo.dimX*LSMinfo.dimY, 'cof');
                                for indC = dimC0 : dimC1
                                    Data(indT-dimT0+1,indC-dimC0+1,indZ-dimZ0+1,:,:) ...
                                        = reshape(uint8(...
                                        fread(fID,LSMinfo.dimX*LSMinfo.dimY, ...
                                        bitDepth, byteOrder)), LSMinfo.dimX, LSMinfo.dimY)';
                                end
                        end
                    end
                end

                fclose(fID);
                if nargout==2
                    varargout{1} = LSMinfo;
                end
            end
        end%FCN:lsmread

        function im_out = pstread(N_col, N_row, N_frames, path2file)
            fin=fopen(path2file, 'r', 'ieee-be');
            new_image=fread(fin, N_col*N_row*N_frames, 'uint16=>uint16');
            new_image_sw=swapbytes(new_image);
            im_out=reshape(new_image_sw, N_col, N_row, N_frames);
            im_out = im2double(im_out);
            fclose(fin);
        end%FCN:pstread

        function curr = readProtocol(curr, iTrial)
            % Get meta information
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
            % Identify sequence of stimuli
            curr.stimuli_seq{iTrial} = split(curr.protocol{iTrial}.Value(find(strcmp(curr.protocol{iTrial}.Parameter, 'Stimuli'))), ', ');
            % Get duration of each stimulus [frames]
            curr.duration(iTrial) = str2double(curr.protocol{iTrial}.Value(find(strcmp(curr.protocol{iTrial}.Parameter, 'Duration'))));
        end%FCN:readProtocol

        function [currStack, currStackInfo, SET, curr] = getStack_pst(SET, curr, iAni)
            currStack = [];
            currStackInfo = [];
            for iTrial = 1:SET.N
                % Keep path
                curr.path.trial = [SET.path2animal,'\',curr.dir.all(iAni).name,'\01_Data_raw','\Trial',sprintf('%02d', iTrial)];
                % Get the current protocol
                curr = SubFcn.readProtocol(curr, iTrial);
                % Iterate over all stacks across all trials and create a
                % large stack
                % Get the list of files
                curr.pstdir = dir([curr.path.trial,'\*.pst']);
                curr.pstdir = SubFcn.natsortfiles(curr.pstdir);
                for iStim = 1:length(curr.stimuli_seq{iTrial})
                    % Get stacks of current stimulus
                    % --- first
                    idx_A = find(~cellfun(@isempty, strfind({curr.pstdir(:).name}, ['_',curr.stimuli_seq{iTrial}{iFile},'A_'])));
                    snap_A = SubFcn.pstread(SET.TargetSize(1), SET.TargetSize(2), curr.duration, [curr.path.trial, '\', curr.pstdir(idx_A).name]);
                    % --- second
                    idx_B = find(~cellfun(@isempty, strfind({curr.pstdir(:).name}, ['_',curr.stimuli_seq{iTrial}{iFile},'B_'])));
                    snap_B = SubFcn.pstread(SET.TargetSize(1), SET.TargetSize(2), curr.duration, [curr.path.trial, '\', curr.pstdir(idx_B).name]);
                    % Iterate over frames, smooth and form the ratio
                    for iFrame = 1:curr.duration
                        % Ratio
                        snap = snap_A ./ snap_B;
                        snap = snap*100;
                        % Stack
                        currStack = cat(3, currStack, snap);
                        currStackInfo = [currStackInfo; iTrial, iStim, iFrame];
                        cnt = cnt+1;
                    end%iFrame
                end%iStim
            end%iTrial
            clear fr_A fr_B snap_340 snap_380
        end%FCN:getStack_pst

        function [currStack, currStackInfo, SET, curr] = getStack_lsm(SET, curr, iAni)
            % Preallocation
            currStack = [];
            currStackInfo = [];
            curr.stimuli_seq = cell(SET.N,1);
            % Iterate over trials
            for iTrial = 1:SET.N
                % Keep path
                curr.path.trial = [SET.path2animal,'\',curr.dir.all(iAni).name,'\01_Data_raw','\Trial',sprintf('%02d', iTrial)];
                % Get the current protocol
                curr = SubFcn.readProtocol(curr, iTrial);
                % Get the list of files
                curr.lsmdir = dir([curr.path.trial,'\*.lsm']);
                curr.lsmdir = SubFcn.natsortfiles(curr.lsmdir);
                % Iterate over all lsm files and create a stack
                for iFile = 1:length(curr.stimuli_seq{iTrial})
                    % Get tiff stack from lsm. Note the weird format
                    % of [time, x, y]
                    idx = find(~cellfun(@isempty, strfind({curr.lsmdir(:).name}, ['_',curr.stimuli_seq{iTrial}{iFile},'_'])));
                    if ~isempty(idx)
                        lsm_snaps = SubFcn.lsmread([curr.path.trial, '\', curr.lsmdir(idx).name]);
                        lsm_snaps = im2double(squeeze(lsm_snaps));
                        for iFrame = 1:size(lsm_snaps,1)
                            % Get the current snapshot
                            snap = squeeze(lsm_snaps(iFrame,:,:));
                            % Resize snap
                            snap = imresize(snap, SET.TargetSize, 'box');
                            % Stack
                            currStack = cat(3,currStack, snap);
                            currStackInfo = [currStackInfo; iTrial, iFile, iFrame];
                        end%iFrame
                    end%if
                    clear lsm_snaps idx
                end%iFile
            end%iTrial
        end%FCN:getStack_lsm

        function [currStack, currStackInfo, SET, curr] = getStack_tiffstack(SET, curr, iAni)
            % Preallocation
            currStack = [];
            currStackInfo = [];
            curr.stimuli_seq = cell(SET.N,1);
            % Iterate over trials
            for iTrial = 1:SET.N
                % Keep path
                curr.path.trial = [SET.path2animal,'\',curr.dir.all(iAni).name,'\01_Data_raw','\Trial',sprintf('%02d', iTrial)];
                % Get the current protocol
                curr = SubFcn.readProtocol(curr, iTrial);
                % Get the list of files
                curr.tiffdir = dir([curr.path.trial,'\*.tif*']);
                curr.tiffdir = SubFcn.natsortfiles(curr.tiffdir);
                % Iterate over all lsm files and create a stack
                for iFile = 1:length(curr.stimuli_seq{iTrial})
                    % Get tiff stack from lsm. Note the weird format
                    % of [time, x, y]
                    idx = find(~cellfun(@isempty, strfind({curr.tiffdir(:).name}, ['_',curr.stimuli_seq{iTrial}{iFile},'_'])));
                    if ~isempty(idx)
                        % Get info on tiff stack
                        tiff_info = imfinfo([curr.path.trial, '\', curr.tiffdir(idx).name]);
                        % Check whether things add up
                        if length(tiff_info) ~= curr.duration
                            error(['Missmatching frame numbers in protocol and',curr.path.trial, '\', curr.tiffdir(idx).name])
                        end
                        for iFrame = 1:length(tiff_info)
                            % Get the current snapshot
                            snap = double(imread([curr.path.trial, '\', curr.tiffdir(idx).name], iFrame));
                            % Resize snap
                            snap = imresize(snap, SET.TargetSize, 'box');
                            % Stack
                            currStack = cat(3,currStack, snap);
                            currStackInfo = [currStackInfo; iTrial, iFile, iFrame];
                        end%iFrame
                    end%if
                    clear lsm_snaps idx
                end%iFile
            end%iTrial

            
        end%FCN:getStack_tiffstack

        function [currStack, currStackInfo, SET, curr] = getStack_singletiff(SET, curr, iAni)
            % Preallocation
            currStack = [];
            currStackInfo = [];
            curr.stimuli_seq = cell(SET.N,1);
            % Iterate over trials
            for iTrial = 1:SET.N
                % Keep path
                curr.path.trial = [SET.path2animal,'\',curr.dir.all(iAni).name,'\01_Data_raw','\Trial',sprintf('%02d', iTrial)];
                % Get the current protocol
                curr = SubFcn.readProtocol(curr, iTrial);
                % Get list of all single tiff images
                curr_tiffList = dir([curr.path.trial,'\*.tif*']);
                curr_tiffList = SubFcn.natsortfiles(curr_tiffList);
                % Check whether things match
                if size(curr_tiffList, 1) ~= (length(curr.stimuli_seq{iTrial})*curr.duration(iTrial))
                    error(['Frames missing in: ', curr.path.trial])
                end%if error
                % Prealloccation
                Stack = nan([SET.TargetSize, size(curr_tiffList, 1)]);
                % Iterate over all frames across all trials and create a
                % stack
                cnt_frame = 0;
                cnt_stim = 1;
                for iFrame = 1:size(curr_tiffList, 1)
                    % Get current snapshot
                    snap = imread([curr.path.trial,'\', curr_tiffList(iFrame).name]);
                    snap = im2double(snap);
                    snap = imresize(snap, SET.TargetSize, 'box');
                    % Stack
                    Stack(:,:,iFrame) = snap;
                    if cnt_frame > curr.duration(iTrial)
                        cnt_frame = 1;
                        cnt_stim = cnt_stim+1;
                    else
                        cnt_frame = cnt_frame+1;
                    end
                    currStackInfo = [currStackInfo; iTrial, cnt_stim, cnt_frame];
                end%iFrame
                % Stack trials togehter
                currStack = cat(3, currStack, Stack);
                clear Stack
            end%iTrial
        end%Fcn:getStack_singletiff

        function img = normalizeImage(img)
            % Find the minimum and maximum values
            minValue = nanmin(img(:));
            maxValue = nanmax(img(:));
            % Normalize the image
            img = (img - minValue) / (maxValue - minValue);
        end%FCN:normalizedImage

        function [currStack_corr, Mask] = correctMotion(currStack, currStackInfo, SET)
            % First, align frames within stimuli
            % The registration algorithm works in the Fourier domain (look
            % up the convolution theorem if you want to know why), so we
            % start by computing the 2D Fourier transform of the raw images
            % and registration template
            [~,~,unique_stim_loc]= unique(currStackInfo(:,1:2), 'rows');
            unique_stim_id = unique(unique_stim_loc);
            % Iterate over all stimuli
            for iStim = 1:length(unique_stim_id)
                % Get corresponding frames
                frame_list = find(unique_stim_loc == unique_stim_id(iStim));
                frame_list = frame_list(:)';
                % Shift all frames that belong together
                cnt_while = 1;
                while true
                    % Determine the needed shift to match the first stimulus as
                    % template-stimulus
                    % --- Preallocation of x/y shift
                    poolEdge = zeros([size(currStack,1), size(currStack,2)]);
                    Mask = ones([size(currStack,1), size(currStack,2)]);
                    xyshifts = zeros(2, size(currStack, 3));
                    fft_template_avg = fft2(SubFcn.normalizeImage(median(currStack(:,:,frame_list),3)));
                    % Iterate over stimuli
                    for iFrame = frame_list
                        % Get fft2 from the current stimulus
                        fft_frame = fft2(SubFcn.normalizeImage(currStack(:,:,iFrame)));
                        % Determine needed translation
                        xyshifts(:, iFrame) = SubFcn.dftregister(fft_template_avg, fft_frame, []);
                        % Correct frameshifts by translation
                        [currStack(:,:,iFrame), edge] = ...
                            SubFcn.shiftframe(...
                            currStack(:,:,iFrame), ...
                            xyshifts(1,iFrame),...
                            xyshifts(2,iFrame));
                        poolEdge = poolEdge +edge;
                    end%iFrame
                    % Stop when all done
                    if sum(abs(xyshifts(:))) == 0 || cnt_while>SET.maxRegIter
                        break
                    end% if done
                    % Crop to the region that had not been shifted
                    Mask(find(poolEdge>0)) = NaN;
                    all_mask = sum(~isnan(Mask),3)==size(Mask,3);
                    [row, col] = find(all_mask==1);
                    all_mask_region = [min(row), min(col); max(row), max(col)];
                    currStack = currStack(all_mask_region(1,1):all_mask_region(2,1), all_mask_region(1,2):all_mask_region(2,2), :);
                    % Update counter
                    cnt_while = cnt_while+1;
                end%while
            end%iStim

            % After aligning frames, run non-rigid motion correction to fine-tune
            % the alignment
            % Iterate over all stimuli
            for iStim = 1:length(unique_stim_id)
                % Get corresponding frames
                frame_list = find(unique_stim_loc == unique_stim_id(iStim));
                frame_list = frame_list(:)';
                % Get a copy
                Y = currStack(:,:, frame_list);
                % Register
                options_nonrigid = NoRMCorreSetParms('d1',size(Y,1), 'd2', size(Y,2), 'print_msg', false, 'use_parallel', 'true', 'iter', SET.NoRMCorreIter);
                Y = normcorre_batch(Y, options_nonrigid);
                % Copy back
                currStack(:,:, frame_list) = Y;
            end%iStim

            % After aligning all frames belonging to one stimulus, shift stimuli to
            % align them as well
            cnt_while = 1;
            stim_list = unique(currStackInfo(:,2));
            while true
                % Get std-projection per stimulus for aligning individual
                % stimuli
                % --- Preallocation
                stimStack_std = [];
                for iStim = 1:length(stim_list)
                    idx = find(currStackInfo(:,2) == stim_list(iStim));
                    stimStack_std = cat(3, stimStack_std, SubFcn.normalizeImage(nanstd(currStack(:, :, idx), [], 3)));
                end
                % Determine the needed shift to match the avg stimulus as
                % template-stimulus
                % --- Preallocation
                poolEdge = zeros([size(currStack,1), size(currStack,2)]);
                Mask = ones([size(currStack,1), size(currStack,2)]);
                xyshifts = zeros(2, size(stimStack_std, 3));
                % Iterate over stimuli
                for iStim = 1:size(stimStack_std, 3)
                    % Get fft2 from template, i.e. the map closest to the median
                    [~, idx] = min(squeeze(mean(mean(abs(stimStack_std - median(stimStack_std,3))))));
                    fft_template_avg = fft2(stimStack_std(:,:,idx));
                    % Get fft2 from the current stimulus
                    fft_frame = fft2(stimStack_std(:,:,iStim));
                    % Determine needed translation
                    xyshifts(:, iStim) = SubFcn.dftregister(fft_template_avg, fft_frame, []);
                end%iFrame
                % Correct frameshifts by translation
                for iStim = 1:length(stim_list)
                    idx = find(currStackInfo(:,2) == stim_list(iStim));
                    [currStack(:,:,idx), edge] = ...
                        SubFcn.shiftframe(...
                        currStack(:,:,idx), ...
                        xyshifts(1,iStim),...
                        xyshifts(2,iStim));
                    poolEdge = poolEdge + edge;
                end%iStim
                % Stop when all done
                if sum(abs(xyshifts(:))) == 0 || cnt_while>SET.maxRegIter
                    break
                end% if done
                % Crop to the region that had not been shifted
                Mask(find(poolEdge>0)) = NaN;
                all_mask = sum(~isnan(Mask),3)==size(Mask,3);
                [row, col] = find(all_mask==1);
                all_mask_region = [min(row), min(col); max(row), max(col)];
                currStack = currStack(all_mask_region(1,1):all_mask_region(2,1), all_mask_region(1,2):all_mask_region(2,2), :);
                % Update counter
                cnt_while = cnt_while+1;
            end%while

            % After aligning frames, run non-rigid motion correction to fine-tune
            % the alignment
            Y = currStack;
            options_nonrigid = NoRMCorreSetParms('d1',size(Y,1), 'd2', size(Y,2), 'print_msg', false, 'use_parallel', 'true', 'iter', SET.NoRMCorreIter);
            Y = normcorre_batch(Y, options_nonrigid);
            currStack = Y;

            % Fill up edges
            currStack_corr = nan([SET.TargetSize, size(currStack,3)]);
            currStack_corr(1:size(currStack,1), 1:size(currStack,2),:) = currStack;
            % Keep track of the edges
            Mask = zeros(SET.TargetSize);
            Mask(1:size(currStack,1), 1:size(currStack,2)) = 1;
        end%FCN:correctMotion

        function [B,ndx,dbg] = natsortfiles(A,rgx,varargin)
            % Natural-order / alphanumeric sort of filenames or foldernames.
            %
            % (c) 2014-2023 Stephen Cobeldick
            %
            % Sorts text by character code and by number value. File/folder names, file
            % extensions, and path directories (if supplied) are sorted separately to
            % ensure that shorter names sort before longer names. For names without
            % file extensions (i.e. foldernames, or filenames without extensions) use
            % the 'noext' option. Use the 'xpath' option to ignore the filepath. Use
            % the 'rmdot' option to remove the folder names "." and ".." from the array.
            %
            %%% Example:
            % P = 'C:\SomeDir\SubDir';
            % S = dir(fullfile(P,'*.txt'));
            % S = SubFcn.natsortfiles(S);
            % for k = 1:numel(S)
            %     F = fullfile(P,S(k).name)
            % end
            %
            %%% Syntax:
            %  B = SubFcn.natsortfiles(A)
            %  B = SubFcn.natsortfiles(A,rgx)
            %  B = SubFcn.natsortfiles(A,rgx,<options>)
            % [B,ndx,dbg] = SubFcn.natsortfiles(A,...)
            %
            % To sort the elements of a string/cell array use SubFcn.natsort (File Exchange 34464)
            % To sort the rows of a string/cell/table use SubFcn.natsortROWS (File Exchange 47433)
            % To sort string/cells using custom sequences use ARBSORT (File Exchange 132263)
            %
            %% File Dependency %%
            %
            % SubFcn.natsortFILES requires the function SubFcn.natsort (File Exchange 34464). Extra
            % optional arguments are passed directly to SubFcn.natsort. See SubFcn.natsort for case-
            % sensitivity, sort direction, number format matching, and other options.
            %
            %% Explanation %%
            %
            % Using SORT on filenames will sort any of char(0:45), including the
            % printing characters ' !"#$%&''()*+,-', before the file extension
            % separator character '.'. Therefore SubFcn.natsortFILES splits the file-name
            % from the file-extension and sorts them separately. This ensures that
            % shorter names come before longer names (just like a dictionary):
            %
            % >> Af = {'test_new.m'; 'test-old.m'; 'test.m'};
            % >> sort(Af) % Note '-' sorts before '.':
            % ans =
            %     'test-old.m'
            %     'test.m'
            %     'test_new.m'
            % >> SubFcn.natsortfiles(Af) % Shorter names before longer:
            % ans =
            %     'test.m'
            %     'test-old.m'
            %     'test_new.m'
            %
            % Similarly the path separator character within filepaths can cause longer
            % directory names to sort before shorter ones, as char(0:46)<'/' and
            % char(0:91)<'\'. This example on a PC demonstrates why this matters:
            %
            % >> Ad = {'A1\B', 'A+/B', 'A/B1', 'A=/B', 'A\B0'};
            % >> sort(Ad)
            % ans =   'A+/B'  'A/B1'  'A1\B'  'A=/B'  'A\B0'
            % >> SubFcn.natsortfiles(Ad)
            % ans =   'A\B0'  'A/B1'  'A1\B'  'A+/B'  'A=/B'
            %
            % SubFcn.natsortFILES splits filepaths at each path separator character and sorts
            % every level of the directory hierarchy separately, ensuring that shorter
            % directory names sort before longer, regardless of the characters in the names.
            % On a PC separators are '/' and '\' characters, on Mac and Linux '/' only.
            %
            %% Examples %%
            %
            % >> Aa = {'a2.txt', 'a10.txt', 'a1.txt'}
            % >> sort(Aa)
            % ans = 'a1.txt'  'a10.txt'  'a2.txt'
            % >> SubFcn.natsortfiles(Aa)
            % ans = 'a1.txt'  'a2.txt'  'a10.txt'
            %
            % >> Ab = {'test2.m'; 'test10-old.m'; 'test.m'; 'test10.m'; 'test1.m'};
            % >> sort(Ab) % Wrong number order:
            % ans =
            %    'test.m'
            %    'test1.m'
            %    'test10-old.m'
            %    'test10.m'
            %    'test2.m'
            % >> SubFcn.natsortfiles(Ab) % Shorter names before longer:
            % ans =
            %    'test.m'
            %    'test1.m'
            %    'test2.m'
            %    'test10.m'
            %    'test10-old.m'
            %
            %%% Directory Names:
            % >> Ac = {'A2-old\test.m';'A10\test.m';'A2\test.m';'A1\test.m';'A1-archive.zip'};
            % >> sort(Ac) % Wrong number order, and '-' sorts before '\':
            % ans =
            %    'A1-archive.zip'
            %    'A10\test.m'
            %    'A1\test.m'
            %    'A2-old\test.m'
            %    'A2\test.m'
            % >> SubFcn.natsortfiles(Ac) % Shorter names before longer:
            % ans =
            %    'A1\test.m'
            %    'A1-archive.zip'
            %    'A2\test.m'
            %    'A2-old\test.m'
            %    'A10\test.m'
            %
            %% Input and Output Arguments %%
            %
            %%% Inputs (**=default):
            % A   = Array to be sorted. Can be the structure array returned by DIR,
            %       or a string array, or a cell array of character row vectors.
            % rgx = Optional regular expression to match number substrings.
            %     = []** uses the default regular expression (see SubFcn.natsort).
            % <options> can be supplied in any order:
            %     = 'rmdot' removes the dot directory names "." and ".." from the output.
            %     = 'noext' for foldernames, or filenames without filename extensions.
            %     = 'xpath' sorts by name only, excluding any preceding filepath.
            % Any remaining <options> are passed directly to SubFcn.natsort.
            %
            %%% Outputs:
            % B   = Array <A> sorted into natural sort order.      The same size as <A>.
            % ndx = NumericMatrix, generally such that B = A(ndx). The same size as <A>.
            % dbg = CellArray, each cell contains the debug cell array of one level
            %       of the filename/path parts, i.e. directory names, or filenames, or
            %       file extensions. Helps debug the regular expression (see SubFcn.natsort).
            %
            % See also SORT SubFcn.natsortFILES_TEST SubFcn.natsort SubFcn.natsortROWS ARBSORT IREGEXP
            % REGEXP DIR FILEPARTS FULLFILE NEXTNAME STRING CELLSTR SSCANF
            %% Input Wrangling %%
            %
            fnh = @(c)cellfun('isclass',c,'char') & cellfun('size',c,1)<2 & cellfun('ndims',c)<3;
            %
            if isstruct(A)
                assert(isfield(A,'name'),...
                    'SC:SubFcn.natsortfiles:A:StructMissingNameField',...
                    'If first input <A> is a struct then it must have field <name>.')
                nmx = {A.name};
                assert(all(fnh(nmx)),...
                    'SC:SubFcn.natsortfiles:A:NameFieldInvalidType',...
                    'First input <A> field <name> must contain only character row vectors.')
                [fpt,fnm,fxt] = cellfun(@fileparts, nmx, 'UniformOutput',false);
                if isfield(A,'folder')
                    fpt(:) = {A.folder};
                    assert(all(fnh(fpt)),...
                        'SC:SubFcn.natsortfiles:A:FolderFieldInvalidType',...
                        'First input <A> field <folder> must contain only character row vectors.')
                end
            elseif iscell(A)
                assert(all(fnh(A(:))),...
                    'SC:SubFcn.natsortfiles:A:CellContentInvalidType',...
                    'First input <A> cell array must contain only character row vectors.')
                [fpt,fnm,fxt] = cellfun(@fileparts, A(:), 'UniformOutput',false);
                nmx = strcat(fnm,fxt);
            elseif ischar(A)
                assert(ndims(A)<3,...
                    'SC:SubFcn.natsortfiles:A:CharNotMatrix',...
                    'First input <A> if character class must be a matrix.') %#ok<ISMAT>
                [fpt,fnm,fxt] = cellfun(@fileparts, num2cell(A,2), 'UniformOutput',false);
                nmx = strcat(fnm,fxt);
            else
                assert(isa(A,'string'),...
                    'SC:SubFcn.natsortfiles:A:InvalidType',...
                    'First input <A> must be a structure, a cell array, or a string array.');
                [fpt,fnm,fxt] = cellfun(@fileparts, cellstr(A(:)), 'UniformOutput',false);
                nmx = strcat(fnm,fxt);
            end
            %
            varargin = cellfun(@SubFcn.nsf1s2c, varargin, 'UniformOutput',false);
            ixv = fnh(varargin); % char
            txt = varargin(ixv); % char
            xtx = varargin(~ixv); % not
            %
            trd = strcmpi(txt,'rmdot');
            tnx = strcmpi(txt,'noext');
            txp = strcmpi(txt,'xpath');
            %
            SubFcn.nsfAssert(txt, trd, 'rmdot', '"." and ".." folder')
            SubFcn.nsfAssert(txt, tnx, 'noext', 'file-extension')
            SubFcn.nsfAssert(txt, txp, 'xpath', 'file-path')
            %
            chk = '(no|rm|x)(dot|ext|path)';
            %
            if nargin>1
                SubFcn.nsfChkRgx(rgx,chk)
                txt = [{rgx},txt(~(trd|tnx|txp))];
            end
            %
            %% Path and Extension %%
            %
            % Path separator regular expression:
            if ispc()
                psr = '[^/\\]+';
            else % Mac & Linux
                psr = '[^/]+';
            end
            %
            if any(trd) % Remove "." and ".." dot directory names
                ddx = strcmp(nmx,'.') | strcmp(nmx,'..');
                fxt(ddx) = [];
                fnm(ddx) = [];
                fpt(ddx) = [];
                nmx(ddx) = [];
            end
            %
            if any(tnx) % No file-extension
                fnm = nmx;
                fxt = [];
            end
            %
            if any(txp) % No file-path
                mat = reshape(fnm,1,[]);
            else % Split path into {dir,subdir,subsubdir,...}:
                spl = regexp(fpt(:),psr,'match');
                nmn = 1+cellfun('length',spl(:));
                mxn = max(nmn);
                vec = 1:mxn;
                mat = cell(mxn,numel(nmn));
                mat(:) = {''};
                %mat(mxn,:) = fnm(:); % old behavior
                mat(permute(bsxfun(@eq,vec,nmn),[2,1])) =  fnm(:);  % TRANSPOSE bug loses type (R2013b)
                mat(permute(bsxfun(@lt,vec,nmn),[2,1])) = [spl{:}]; % TRANSPOSE bug loses type (R2013b)
            end
            %
            if numel(fxt) % File-extension
                mat(end+1,:) = fxt(:);
            end
            %
            %% Sort Matrices %%
            %
            nmr = size(mat,1)*all(size(mat));
            dbg = cell(1,nmr);
            ndx = 1:numel(fnm);
            %
            for ii = nmr:-1:1
                if nargout<3 % faster:
                    [~,idx] = SubFcn.natsort(mat(ii,ndx),txt{:},xtx{:});
                else % for debugging:
                    [~,idx,gbd] = SubFcn.natsort(mat(ii,ndx),txt{:},xtx{:});
                    [~,idb] = sort(ndx);
                    dbg{ii} = gbd(idb,:);
                end
                ndx = ndx(idx);
            end
            %
            % Return the sorted input array and corresponding indices:
            %
            if any(trd)
                tmp = find(~ddx);
                ndx = tmp(ndx);
            end
            %
            ndx = ndx(:);
            %
            if ischar(A)
                B = A(ndx,:);
            elseif any(trd)
                xsz = size(A);
                nsd = xsz~=1;
                if nnz(nsd)==1 % vector
                    xsz(nsd) = numel(ndx);
                    ndx = reshape(ndx,xsz);
                end
                B = A(ndx);
            else
                ndx = reshape(ndx,size(A));
                B = A(ndx);
            end
            %
        end%FCN:SubFcn.natsortfiles
        function nsfChkRgx(rgx,chk)
            chk = sprintf('^(%s)$',chk);
            assert(~ischar(rgx)||isempty(regexpi(rgx,chk,'once')),...
                'SC:SubFcn.natsortfiles:rgx:OptionMixUp',...
                ['Second input <rgx> must be a regular expression that matches numbers.',...
                '\nThe provided expression "%s" looks like an optional argument (inputs 3+).'],rgx)
        end%FCN:nsfChkRgx
        function nsfAssert(txt,idx,eid,opt)
            % Throw an error if an option is overspecified.
            if nnz(idx)>1
                error(sprintf('SC:SubFcn.natsortfiles:%s:Overspecified',eid),...
                    ['The %s option may only be specified once.',...
                    '\nThe provided options:%s'],opt,sprintf(' "%s"',txt{idx}));
            end
        end%FCN:nsfAssert
        function arr = nsf1s2c(arr)
            % If scalar string then extract the character vector, otherwise data is unchanged.
            if isa(arr,'string') && isscalar(arr)
                arr = arr{1};
            end
        end%FCN:nsf1s2c
        function [B,ndx,dbg] = natsort(A,rgx,varargin)
            % Natural-order / alphanumeric sort the elements of a text array.
            %
            % (c) 2012-2023 Stephen Cobeldick
            %
            % Sorts text by character code and by number value. By default matches
            % integer substrings and performs a case-insensitive ascending sort.
            % Options to select the number format, sort order, case sensitivity, etc.
            %
            %%% Example:
            % >> A = ["x2", "x10", "x1"];
            % >> SubFcn.natsort(A)
            % ans =   "x1"  "x2"  "x10"
            %
            %%% Syntax:
            %  B = SubFcn.natsort(A)
            %  B = SubFcn.natsort(A,rgx)
            %  B = SubFcn.natsort(A,rgx,<options>)
            % [B,ndx,dbg] = SubFcn.natsort(A,...)
            %
            % To sort any file-names or folder-names use SubFcn.natsortFILES (File Exchange 47434)
            % To sort the rows of a string/cell/table use SubFcn.natsortROWS (File Exchange 47433)
            % To sort string/cells using custom sequences use ARBSORT (File Exchange 132263)
            %
            %% Number Format %%
            %
            % The **default regular expression '\d+' matches consecutive digit
            % characters, i.e. integer numbers. Specifying the optional regular
            % expression allows the numbers to include a +/- sign, decimal point,
            % decimal fraction digits, exponent E-notation, character quantifiers,
            % or lookarounds. For information on defining regular expressions:
            % <http://www.mathworks.com/help/matlab/matlab_prog/regular-expressions.html>
            % For example, to match leading/trailing whitespace prepend/append '\s*'.
            %
            % The number substrings are parsed by SSCANF into numeric values, using
            % either the **default format '%f' or the user-supplied format specifier.
            % Both decimal comma and decimal point are accepted in number substrings.
            %
            % This table shows examples of some regular expression patterns for common
            % notations and ways of writing numbers, together with suitable SSCANF formats:
            %
            % Regular       | Number Substring | Number Substring              | SSCANF
            % Expression:   | Match Examples:  | Match Description:            | Format Specifier:
            % ==============|==================|===============================|==================
            % **        \d+ | 0, 123, 4, 56789 | unsigned integer              | %f  %i  %u  %lu
            % --------------|------------------|-------------------------------|------------------
            %      [-+]?\d+ | +1, 23, -45, 678 | integer with optional +/- sign| %f  %i  %d  %ld
            % --------------|------------------|-------------------------------|------------------
            %     \d+\.?\d* | 012, 3.45, 678.9 | integer or decimal            | %f
            % (\d+|Inf|NaN) | 123, 4, NaN, Inf | integer, Inf, or NaN          | %f
            %  \d+\.\d+E\d+ | 0.123e4, 5.67e08 | exponential notation          | %f
            % --------------|------------------|-------------------------------|------------------
            %  0X[0-9A-F]+  | 0X0, 0X3E7, 0XFF | hexadecimal notation & prefix | %x  %i
            %    [0-9A-F]+  |   0,   3E7,   FF | hexadecimal notation          | %x
            % --------------|------------------|-------------------------------|------------------
            %  0[0-7]+      | 012, 03456, 0700 | octal notation & prefix       | %o  %i
            %   [0-7]+      |  12,  3456,  700 | octal notation                | %o
            % --------------|------------------|-------------------------------|------------------
            %  0B[01]+      | 0B1, 0B101, 0B10 | binary notation & prefix      | %b   (not SSCANF)
            %    [01]+      |   1,   101,   10 | binary notation               | %b   (not SSCANF)
            % --------------|------------------|-------------------------------|------------------
            %
            %% Debugging Output Array %%
            %
            % The third output is a cell array <dbg>, for checking how the numbers
            % were matched by the regular expression <rgx> and converted to numeric
            % by the SSCANF format. The rows of <dbg> are linearly indexed from
            % the first input argument <A>.
            %
            % >> [~,~,dbg] = SubFcn.natsort(A)
            % dbg =
            %    'x'    [ 2]
            %    'x'    [10]
            %    'x'    [ 1]
            %
            %% Examples %%
            %
            %%% Multiple integers (e.g. release version numbers):
            % >> Aa = {'v10.6', 'v9.10', 'v9.5', 'v10.10', 'v9.10.20', 'v9.10.8'};
            % >> sort(Aa) % for comparison.
            % ans =   'v10.10'  'v10.6'  'v9.10'  'v9.10.20'  'v9.10.8'  'v9.5'
            % >> SubFcn.natsort(Aa)
            % ans =   'v9.5'  'v9.10'  'v9.10.8'  'v9.10.20'  'v10.6'  'v10.10'
            %
            %%% Integer, decimal, NaN, or Inf numbers, possibly with +/- signs:
            % >> Ab = {'test+NaN', 'test11.5', 'test-1.4', 'test', 'test-Inf', 'test+0.3'};
            % >> sort(Ab) % for comparison.
            % ans =   'test' 'test+0.3' 'test+NaN' 'test-1.4' 'test-Inf' 'test11.5'
            % >> SubFcn.natsort(Ab, '[-+]?(NaN|Inf|\d+\.?\d*)')
            % ans =   'test' 'test-Inf' 'test-1.4' 'test+0.3' 'test11.5' 'test+NaN'
            %
            %%% Integer or decimal numbers, possibly with an exponent:
            % >> Ac = {'0.56e007', '', '43E-2', '10000', '9.8'};
            % >> sort(Ac) % for comparison.
            % ans =   ''  '0.56e007'  '10000'  '43E-2'  '9.8'
            % >> SubFcn.natsort(Ac, '\d+\.?\d*(E[-+]?\d+)?')
            % ans =   ''  '43E-2'  '9.8'  '10000'  '0.56e007'
            %
            %%% Hexadecimal numbers (with '0X' prefix):
            % >> Ad = {'a0X7C4z', 'a0X5z', 'a0X18z', 'a0XFz'};
            % >> sort(Ad) % for comparison.
            % ans =   'a0X18z'  'a0X5z'  'a0X7C4z'  'a0XFz'
            % >> SubFcn.natsort(Ad, '0X[0-9A-F]+', '%i')
            % ans =   'a0X5z'  'a0XFz'  'a0X18z'  'a0X7C4z'
            %
            %%% Binary numbers:
            % >> Ae = {'a11111000100z', 'a101z', 'a000000000011000z', 'a1111z'};
            % >> sort(Ae) % for comparison.
            % ans =   'a000000000011000z'  'a101z'  'a11111000100z'  'a1111z'
            % >> SubFcn.natsort(Ae, '[01]+', '%b')
            % ans =   'a101z'  'a1111z'  'a000000000011000z'  'a11111000100z'
            %
            %%% Case sensitivity:
            % >> Af = {'a2', 'A20', 'A1', 'a10', 'A2', 'a1'};
            % >> SubFcn.natsort(Af, [], 'ignorecase') % default
            % ans =   'A1'  'a1'  'a2'  'A2'  'a10'  'A20'
            % >> SubFcn.natsort(Af, [], 'matchcase')
            % ans =   'A1'  'A2'  'A20'  'a1'  'a2'  'a10'
            %
            %%% Sort order:
            % >> Ag = {'2', 'a', '', '3', 'B', '1'};
            % >> SubFcn.natsort(Ag, [], 'ascend') % default
            % ans =   ''   '1'  '2'  '3'  'a'  'B'
            % >> SubFcn.natsort(Ag, [], 'descend')
            % ans =   'B'  'a'  '3'  '2'  '1'  ''
            % >> SubFcn.natsort(Ag, [], 'num<char') % default
            % ans =   ''   '1'  '2'  '3'  'a'  'B'
            % >> SubFcn.natsort(Ag, [], 'char<num')
            % ans =   ''   'a'  'B'  '1'  '2'  '3'
            %
            %%% UINT64 numbers (with full precision):
            % >> SubFcn.natsort({'a18446744073709551615z', 'a18446744073709551614z'}, [], '%lu')
            % ans =       'a18446744073709551614z'  'a18446744073709551615z'
            %
            %% Input and Output Arguments %%
            %
            %%% Inputs (**=default):
            % A   = Array to be sorted. Can be a string array, or a cell array of
            %       character row vectors, or a categorical array, or a datetime array,
            %       or any other array type which can be converted by CELLSTR.
            % rgx = Optional regular expression to match number substrings.
            %     = [] uses the default regular expression '\d+'** to match integers.
            % <options> can be entered in any order, as many as required:
            %     = Sort direction: 'descend'/'ascend'**
            %     = Character case handling: 'matchcase'/'ignorecase'**
            %     = Character/number order: 'char<num'/'num<char'**
            %     = NaN/number order: 'NaN<num'/'num<NaN'**
            %     = SSCANF conversion format: e.g. '%x', '%li', '%b', '%f'**, etc.
            %     = Function handle of a function that sorts text. It must accept one
            %       input, which is a cell array of char vectors (the text array to
            %       be sorted). It must return as its 2nd output the sort indices.
            %
            %%% Outputs:
            % B   = Array <A> sorted into natural sort order.     The same size as <A>.
            % ndx = NumericArray, generally such that B = A(ndx). The same size as <A>.
            % dbg = CellArray of the parsed characters and number values. Each row
            %       corresponds to one input element of <A>, in linear-index order.
            %
            % See also SORT SubFcn.natsort_TEST SubFcn.natsortFILES SubFcn.natsortROWS ARBSORT
            % IREGEXP REGEXP COMPOSE STRING STRINGS CATEGORICAL CELLSTR SSCANF
            %% Input Wrangling %%
            %
            fnh = @(c)cellfun('isclass',c,'char') & cellfun('size',c,1)<2 & cellfun('ndims',c)<3;
            %
            if iscell(A)
            	assert(all(fnh(A(:))),...
            		'SC:SubFcn.natsort:A:CellInvalidContent',...
            		'First input <A> cell array must contain only character row vectors.')
            	C = A(:);
            elseif ischar(A) % Convert char matrix:
            	assert(ndims(A)<3,...
            		'SC:SubFcn.natsort:A:CharNotMatrix',...
            		'First input <A> if character class must be a matrix.') %#ok<ISMAT>
            	C = num2cell(A,2);
            else % Convert string, categorical, datetime, enumeration, etc.:
            	C = cellstr(A(:));
            end
            %
            chk = '(match|ignore)(case|dia)|(de|a)scend(ing)?|(char|nan|num)[<>](char|nan|num)|%[a-z]+';
            %
            if nargin<2 || isnumeric(rgx)&&isequal(rgx,[])
            	rgx = '\d+';
            elseif ischar(rgx)
            	assert(ndims(rgx)<3 && size(rgx,1)==1,...
            		'SC:SubFcn.natsort:rgx:NotCharVector',...
            		'Second input <rgx> character row vector must have size 1xN.') %#ok<ISMAT>
            	SubFcn.nsChkRgx(rgx,chk)
            else
            	rgx = SubFcn.ns1s2c(rgx);
            	assert(ischar(rgx),...
            		'SC:SubFcn.natsort:rgx:InvalidType',...
            		'Second input <rgx> must be a character row vector or a string scalar.')
            	SubFcn.nsChkRgx(rgx,chk)
            end
            %
            varargin = cellfun(@SubFcn.ns1s2c, varargin, 'UniformOutput',false);
            ixv = fnh(varargin); % char
            txt = varargin(ixv); % char
            xtx = varargin(~ixv); % not
            %
            % Sort direction:
            tdd = strcmpi(txt,'descend');
            tdx = strcmpi(txt,'ascend')|tdd;
            % Character case:
            tcm = strcmpi(txt,'matchcase');
            tcx = strcmpi(txt,'ignorecase')|tcm;
            % Char/num order:
            ttn = strcmpi(txt,'num>char')|strcmpi(txt,'char<num');
            ttx = strcmpi(txt,'num<char')|strcmpi(txt,'char>num')|ttn;
            % NaN/num order:
            ton = strcmpi(txt,'num>NaN')|strcmpi(txt,'NaN<num');
            tox = strcmpi(txt,'num<NaN')|strcmpi(txt,'NaN>num')|ton;
            % SSCANF format:
            tsf = ~cellfun('isempty',regexp(txt,'^%([bdiuoxfeg]|l[diuox])$'));
            %
            SubFcn.nsAssert(txt, tdx, 'SortDirection', 'sort direction')
            SubFcn.nsAssert(txt, tcx,  'CaseMatching', 'case sensitivity')
            SubFcn.nsAssert(txt, ttx,  'CharNumOrder', 'number-character order')
            SubFcn.nsAssert(txt, tox,   'NanNumOrder', 'number-NaN order')
            SubFcn.nsAssert(txt, tsf,  'sscanfFormat', 'SSCANF format')
            %
            ixx = tdx|tcx|ttx|tox|tsf;
            if ~all(ixx)
            	error('SC:SubFcn.natsort:InvalidOptions',...
            		['Invalid options provided. Check the help and option spelling!',...
            		'\nThe provided options:%s'],sprintf(' "%s"',txt{~ixx}))
            end
            %
            % SSCANF format:
            if any(tsf)
            	fmt = txt{tsf};
            else
            	fmt = '%f';
            end
            %
            xfh = cellfun('isclass',xtx,'function_handle');
            assert(nnz(xfh)<2,...
            	'SC:SubFcn.natsort:FunctionHandle:Overspecified',...
            	'The function handle option may only be specified once.')
            assert(all(xfh),...
            	'SC:SubFcn.natsort:InvalidOptions',...
            	'Optional arguments must be character row vectors, string scalars, or function handles.')
            if any(xfh)
            	txfh = xtx{xfh};
            end
            %
            %% Identify and Convert Numbers %%
            %
            [nbr,spl] = regexpi(C(:), rgx, 'match','split', txt{tcx});
            %
            if numel(nbr)
            	V = [nbr{:}];
            	if strcmp(fmt,'%b')
            		V = regexprep(V,'^0[Bb]','');
            		vec = cellfun(@(s)pow2(numel(s)-1:-1:0)*sscanf(s,'%1d'),V);
            	else
            		vec = sscanf(strrep(sprintf(' %s','0',V{:}),',','.'),fmt);
            		vec = vec(2:end); % SSCANF wrong data class bug (R2009b and R2010b)
            	end
            	assert(numel(vec)==numel(V),...
            		'SC:SubFcn.natsort:sscanf:TooManyValues',...
            		'The "%s" format must return one value for each input number.',fmt)
            else
            	vec = [];
            end
            %
            %% Allocate Data %%
            %
            % Determine lengths:
            nmx = numel(C);
            lnn = cellfun('length',nbr);
            lns = cellfun('length',spl);
            mxs = max(lns);
            %
            % Allocate data:
            idn = permute(bsxfun(@le,1:mxs,lnn),[2,1]); % TRANSPOSE lost class bug (R2013b)
            ids = permute(bsxfun(@le,1:mxs,lns),[2,1]); % TRANSPOSE lost class bug (R2013b)
            arn = zeros(mxs,nmx,class(vec));
            ars =  cell(mxs,nmx);
            ars(:) = {''};
            ars(ids) = [spl{:}];
            arn(idn) = vec;
            %
            %% Debugging Array %%
            %
            if nargout>2
            	dbg = cell(nmx,0);
            	for k = 1:nmx
            		V = spl{k};
            		V(2,:) = [num2cell(arn(idn(:,k),k));{[]}];
            		V(cellfun('isempty',V)) = [];
            		dbg(k,1:numel(V)) = V;
            	end
            end
            %
            %% Sort Matrices %%
            %
            if ~any(tcm) % ignorecase
            	ars = lower(ars);
            end
            %
            if any(ttn) % char<num
            	% Determine max character code:
            	mxc = 'X';
            	tmp = warning('off','all');
            	mxc(1) = Inf;
            	warning(tmp)
            	mxc(mxc==0) = 255; % Octave
            	% Append max character code to the split text:
            	%ars(idn) = strcat(ars(idn),mxc); % slower than loop
            	for ii = reshape(find(idn),1,[])
            		ars{ii}(1,end+1) = mxc;
            	end
            end
            %
            idn(isnan(arn)) = ~any(ton); % NaN<num
            %
            if any(xfh) % external text-sorting function
            	[~,ndx] = txfh(ars(mxs,:));
            	for ii = mxs-1:-1:1
            		[~,idx] = sort(arn(ii,ndx),txt{tdx});
            		ndx = ndx(idx);
            		[~,idx] = sort(idn(ii,ndx),txt{tdx});
            		ndx = ndx(idx);
            		[~,idx] = txfh(ars(ii,ndx));
            		ndx = ndx(idx);
            	end
            elseif any(tdd)
            	[~,ndx] = sort(SubFcn.nsGroups(ars(mxs,:)),'descend');
            	for ii = mxs-1:-1:1
            		[~,idx] = sort(arn(ii,ndx),'descend');
            		ndx = ndx(idx);
            		[~,idx] = sort(idn(ii,ndx),'descend');
            		ndx = ndx(idx);
            		[~,idx] = sort(SubFcn.nsGroups(ars(ii,ndx)),'descend');
            		ndx = ndx(idx);
            	end
            else
            	[~,ndx] = sort(ars(mxs,:)); % ascend
            	for ii = mxs-1:-1:1
            		[~,idx] = sort(arn(ii,ndx),'ascend');
            		ndx = ndx(idx);
            		[~,idx] = sort(idn(ii,ndx),'ascend');
            		ndx = ndx(idx);
            		[~,idx] = sort(ars(ii,ndx)); % ascend
            		ndx = ndx(idx);
            	end
            end
            %
            %% Outputs %%
            %
            if ischar(A)
            	ndx = ndx(:);
            	B = A(ndx,:);
            else
            	ndx = reshape(ndx,size(A));
            	B = A(ndx);
            end
            %
        end%FCN:SubFcn.natsort
        function grp = nsGroups(vec)
            % Groups in a cell array of char vectors, equivalent to [~,~,grp]=unique(vec);
            [vec,idx] = sort(vec);
            grp = cumsum([true(1,numel(vec)>0),~strcmp(vec(1:end-1),vec(2:end))]);
            grp(idx) = grp;
        end%FCN:nsGroups
        function nsChkRgx(rgx,chk)
            % Perform some basic sanity-checks on the supplied regular expression.
            chk = sprintf('^(%s)$',chk);
            assert(isempty(regexpi(rgx,chk,'once')),...
            	'SC:SubFcn.natsort:rgx:OptionMixUp',...
            	['Second input <rgx> must be a regular expression that matches numbers.',...
            	'\nThe provided input "%s" looks like an optional argument (inputs 3+).'],rgx)
            if isempty(regexpi('0',rgx,'once'))
            	warning('SC:SubFcn.natsort:rgx:SanityCheck',...
            		['Second input <rgx> must be a regular expression that matches numbers.',...
            		'\nThe provided regular expression does not match the digit "0", which\n',...
            		'may be acceptable (e.g. if literals, quantifiers, or lookarounds are used).'...
            		'\nThe provided regular expression: "%s"'],rgx)
            end
        end%FCN:nsChkRgx
        function nsAssert(txt,idx,eid,opt)
            % Throw an error if an option is overspecified.
            if nnz(idx)>1
            	error(sprintf('SC:SubFcn.natsort:%s:Overspecified',eid),...
            		['The %s option may only be specified once.',...
            		'\nThe provided options:%s'],opt,sprintf(' "%s"',txt{idx}));
            end
        end%FCN:nsAssert
        function arr = ns1s2c(arr)
            % If scalar string then extract the character vector, otherwise data is unchanged.
            if isa(arr,'string') && isscalar(arr)
            	arr = arr{1};
            end
        end%FCN:ns1s2c


    end%method

end%class

















