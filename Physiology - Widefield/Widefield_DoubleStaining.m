%% Widefield_DoubleStaining
% This script loads the anatomical double staining data to visualize it.
%
% Version:
% 05-Jan-2024 (R2023a)

% Prepare
clc; clear all; close all
warning('off')

% Add toolboxes
% A MATLAB toolbox for exporting publication quality figures
% (https://github.com/altmany/export_fig)
addpath(genpath('...\export_fig'))

%% Settings

% Set base path to find the data
SET.Stack_path = {...
    '...\220909_soli_PN_ORN\';...
    '...\220915_soli_PN_ORN\'};

% Set names of folders that contain individual tiff Stacks for each channel
SET.Stack_chan = {'tiff_C1'; 'tiff_C2'};

% Set the voxel size
SET.VoxelSize = [...
    0.1186094, 0.1186094, 0.6507171;...
    0.1186094, 0.1186094, 0.6507171];

% Set nanstd for 2D Gaussian filtering ('gauss'), or box size for
% 2D median filtering (med).
% Note that spatial filtering is applied before resizing the image.
SET.SmoothData_spat_type = 'med';
SET.SmoothData_spat = 11;

% Set target size of images to decreasefilesize.
SET.OriginalSize = [4096, 4096];
SET.TargetSize = [512 512];

% Correct voxel size
SET.VoxelSize(:,1) = SET.VoxelSize(:,1)*(SET.OriginalSize(1)/SET.TargetSize(1));
SET.VoxelSize(:,2) = SET.VoxelSize(:,2)*(SET.OriginalSize(2)/SET.TargetSize(2));

% Number of sub-Stacks (mini max-projections)
SET.SubStackSteps = 11;

%% Pool images

% Iterate over all stacks
for iStack = 1:length(SET.Stack_path)

    % Iterate over both color channels
    for iChan = 1:length(SET.Stack_chan)

        % Get content of current directory
        curr_dir = dir([SET.Stack_path{iStack}, SET.Stack_chan{iChan}]);

        % Preallocation
        Stack.(SET.Stack_chan{iChan}) = nan([SET.TargetSize,length(curr_dir)]);
        cnt = 1;

        % Iterate over all files
        hWait = waitbar(0, 'Please wait');
        for iFile = 1:length(curr_dir)
            waitbar(iFile/length(curr_dir), hWait)
            % Check whether it is a tif file
            if contains(curr_dir(iFile).name, '.tif')
                % Get raw image
                img = double(imread([SET.Stack_path{iStack}, SET.Stack_chan{iChan},'\',curr_dir(iFile).name]));
                % Spatial filtering
                switch SET.SmoothData_spat_type
                    case 'gauss'
                        img = imgaussfilt(img, SET.SmoothData_spat);
                    case 'med'
                        img = medfilt2(img, [SET.SmoothData_spat SET.SmoothData_spat]);
                end%switch
                % Resize image
                img = imresize(img, SET.TargetSize, 'box');
                % Pool
                Stack.(SET.Stack_chan{iChan})(:,:,cnt) = img;
                cnt = cnt+1;
            end%if
        end%iFile
        % Trim
        Stack.(SET.Stack_chan{iChan}) = Stack.(SET.Stack_chan{iChan})(:,:,1:cnt-1);

        % Normalize each frame
        Stack_smooth.(SET.Stack_chan{iChan}) = Stack.(SET.Stack_chan{iChan});
        for iFrame = 1:size(Stack.(SET.Stack_chan{iChan}),3)
            temp = Stack_smooth.(SET.Stack_chan{iChan})(:,:,iFrame);
            temp = temp - nanmin(temp(:));
            temp = sqrt(temp);
            temp = temp - nanmean(temp(:));
            temp = temp / nanstd(temp(:));
            Stack_smooth.(SET.Stack_chan{iChan})(:,:,iFrame) = temp;
        end%iFrame

        % Show every n-th subStack
        SubStackRegions = 1:SET.SubStackSteps:size(Stack.(SET.Stack_chan{iChan}),3);
        SubStack.(SET.Stack_chan{iChan}) = nan([SET.TargetSize, SET.SubStackSteps]);
        for iS = 1:length(SubStackRegions)-1
            hFig = figure('units', 'normalized', 'Position', [0 0 1 1], 'Color', 'w');
            subplot(1,2,1)
            SubStack.(SET.Stack_chan{iChan})(:,:,iS) = max(Stack.(SET.Stack_chan{iChan})(:,:,SubStackRegions(iS):SubStackRegions(iS+1)-1),[],3);
            imagesc(SubStack.(SET.Stack_chan{iChan})(:,:,iS))
            axis equal tight off
            colormap gray
            img = SubStack.(SET.Stack_chan{iChan})(:,:,iS);
            clim([quantile(img(:), 0.01), quantile(img(:), 0.99)])
            title(num2str(SubStackRegions(iS):SubStackRegions(iS+1)-1))

            subplot(1,2,2)
            SubStack_smooth.(SET.Stack_chan{iChan})(:,:,iS) = max(Stack_smooth.(SET.Stack_chan{iChan})(:,:,SubStackRegions(iS):SubStackRegions(iS+1)-1),[],3);
            imagesc(SubStack_smooth.(SET.Stack_chan{iChan})(:,:,iS))
            axis equal tight off
            colormap gray
            img = SubStack_smooth.(SET.Stack_chan{iChan})(:,:,iS);
            clim([quantile(img(:), 0.01), quantile(img(:), 0.99)])
            title(num2str(SubStackRegions(iS):SubStackRegions(iS+1)-1))

            % Save
            export_fig([SET.Stack_path{iStack}, SET.Stack_chan{iChan},'\Stack_',sprintf('%03d', iS)], '-pdf', '-painters')

            close(hFig)
        end%iS

        % Show the max-projection in different dimensions
        hFig = figure('units', 'normalized', 'Position', [0 0 1 1], 'Color', 'w');
        % --- condense depth
        subplot(2,2,1)
        img_top = max(Stack_smooth.(SET.Stack_chan{iChan}), [], 3);
        imagesc(img_top)
        clim([quantile(img_top(:), 0.01), quantile(img_top(:), 0.99)])
        xticks([])
        yticks([])
        xlabel(['width (',num2str(round(SET.VoxelSize(iStack,1)*size(Stack_smooth.(SET.Stack_chan{iChan}),1),2)), ' um)'])
        ylabel(['height (',num2str(round(SET.VoxelSize(iStack,2)*size(Stack_smooth.(SET.Stack_chan{iChan}),2),2)), ' um)'])
        axis equal
        ax = gca;
        ax.DataAspectRatio = [SET.VoxelSize(1), SET.VoxelSize(2), 1]/max(SET.VoxelSize([1,2]));
        colormap gray

        % --- condense height
        subplot(2,2,2)
        img_right = squeeze(max(Stack_smooth.(SET.Stack_chan{iChan}), [], 1))';
        imagesc(img_right)
        clim([quantile(img_right(:), 0.01), quantile(img_right(:), 0.99)])
        xticks([])
        yticks([])
        ylabel(['depth (',num2str(round(SET.VoxelSize(iStack,3)*size(Stack_smooth.(SET.Stack_chan{iChan}),3),2)), ' um)'])
        xlabel(['width (',num2str(round(SET.VoxelSize(iStack,1)*size(Stack_smooth.(SET.Stack_chan{iChan}),1),2)), ' um)'])
        ax = gca;
        ax.DataAspectRatio = [SET.VoxelSize(iStack,3), SET.VoxelSize(iStack,1), 1];
        colormap gray
        clear ax

        % --- condense width
        subplot(2,2,3)
        img_left = squeeze(max(Stack_smooth.(SET.Stack_chan{iChan}), [], 2))';
        imagesc(img_left)
        clim([quantile(img_left(:), 0.01), quantile(img_left(:), 0.99)])
        xticks([])
        yticks([])
        ylabel(['depth (',num2str(round(SET.VoxelSize(iStack,3)*size(Stack_smooth.(SET.Stack_chan{iChan}),3),2)), ' um)'])
        xlabel(['height (',num2str(round(SET.VoxelSize(iStack,2)*size(Stack_smooth.(SET.Stack_chan{iChan}),2),2)), ' um)'])
        ax = gca;
        ax.DataAspectRatio = [SET.VoxelSize(iStack,3), SET.VoxelSize(iStack,2), 1];
        colormap gray
        clear ax

        % Display the image
        subplot(2,2,4); hold on
        % Define vertices for the surfaces of the cube
        % --- Top
        x_top = [0 1; 0 1];
        y_top = [1 1; 0 0];
        z_top = [1 1; 1 1];
        % --- Left
        x_left = [0 0; 0 0];
        y_left = [1 0; 1 0];
        z_left = [1 1; 0 0];
        % --- Right
        x_right = [0 1; 0 1];
        y_right = [0 0; 0 0];
        z_right = [1 1; 0 0];
        % Add the left image to the cube
        h_left = surf(x_left, y_left, z_left, img_left,'FaceColor','texturemap','EdgeColor','k');
        % Add the right image to the cube
        h_right = surf(x_right, y_right, z_right, img_right,'FaceColor','texturemap','EdgeColor','k');
        % Add the top image to the cube
        h_top = surf(x_top,y_top,z_top,img_top,'FaceColor','texturemap','EdgeColor','k');
        % Cosmetics
        view(3);
        ax = gca;
        ax.DataAspectRatio = [SET.VoxelSize(iStack,3), SET.VoxelSize(iStack,3), SET.VoxelSize(iStack,1)];
        view([-45 22.5])
        xlabel(['width (',num2str(round(SET.VoxelSize(iStack,1)*size(Stack_smooth.(SET.Stack_chan{iChan}),1),2)), ' um)'])
        ylabel(['height (',num2str(round(SET.VoxelSize(iStack,2)*size(Stack_smooth.(SET.Stack_chan{iChan}),2),2)), ' um)'])
        zlabel(['depth (',num2str(round(SET.VoxelSize(iStack,3)*size(Stack_smooth.(SET.Stack_chan{iChan}),3),2)), ' um)'])
        xticks([])
        yticks([])
        zticks([])
        colormap gray;

        % Save
        export_fig([SET.Stack_path{iStack}, SET.Stack_chan{iChan},'\Stack_projections'], '-pdf', '-painters')

        % Close waitbar
        close(hWait)

        % Save
        save([SET.Stack_path{iStack}, SET.Stack_chan{iChan},'\Stack.mat'], 'Stack', 'SubStack', 'SET')

    end%iChan

    % Combine both channels
    c1 = load([SET.Stack_path{iStack}, SET.Stack_chan{1},'\Stack.mat']);
    c2 = load([SET.Stack_path{iStack}, SET.Stack_chan{2},'\Stack.mat']);

    figure('units', 'normalized', 'Position', [0 0 1 1], 'Color', 'w')
    cols1 = [linspace(0,0,1000)', linspace(0,1,1000)', linspace(0,0,1000)'];
    cols2 = [linspace(0,1,1000)', linspace(0,0,1000)', linspace(0,1,1000)'];
    stack_id = [1 5 10 15 20];
    % -------------------------------------------------------------------------
    for iS = 1:length(stack_id)
        % --- PNs
        subplot(3,5,1+iS-1)
        img1 = c1.SubStack.tiff_C1(:,:,stack_id(iS));
        img1(img1>quantile(img1(:),0.99)) = quantile(img1(:),0.99);
        img1(img1<quantile(img1(:),0.1)) = quantile(img1(:),0.1);
        imagesc(img1)
        axis equal tight off
        ax = gca;
        colormap(ax, cols1)
        % --- ORNs
        subplot(3,5,6+iS-1)
        img2 = c2.SubStack.tiff_C2(:,:,stack_id(iS));  
        img2(img2>quantile(img2(:),0.99)) = quantile(img2(:),0.99);
        img2(img2<quantile(img2(:),0.1)) = quantile(img2(:),0.1);
        imagesc(img2)
        % clim([quantile(img2(:), 0.1), quantile(img2(:), 0.99)])
        axis equal tight off
        ax = gca;
        colormap(ax, cols2)
        % --- Overlay
        subplot(3,5,11+iS-1)
        imagesc(imfuse(img1, img2))
        axis equal tight off
    end

    export_fig([SET.Stack_path{iStack},'StackOverlay'], '-pdf', '-painters')

    %%
    figure('units', 'normalized', 'Position', [0 0 1 1], 'Color', 'w')
    imagesc(c1.SubStack.tiff_C1(:,:,5))
    axis equal tight off
    colormap parula
    hold on
    locMax.x = [];
    locMax.y = [];
    while true
        pos = ginput(1);
        if any(pos>size(c1.SubStack.tiff_C1,1))
            break
        elseif any(pos<1)
            break
        else
            locMax.x = [locMax.x; pos(1)];
            locMax.y = [locMax.y; pos(2)];
            plot(locMax.x(end), locMax.y(end), 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none')
        end%if
    end%while
    save([SET.Stack_path{iStack},'PN_loc.mat'], 'locMax')

    %% Clean
    clearvars -except iStack SET

end%iStack