%% Iterate over all animals and kick out invalid trials
% CIAO_01b_kickout prevents the pooling of different trials in later steps. 
% For this, the respective trial number is deleted from the 
% meta.valid_trials variable. 
% Provide a list with invalid trials in a textfile (invalidTrials.txt) and
% locate it in the main data folder. The file should contain two columns
% named "Animal"and "Trial". Animal IDs should be proveded as A###, trial
% numbers as a simple integer.
% An examplary table could look like the following:
%     Animal,Trial
%     A026,1
%     A038,1
%     A123,1
%
%
% Version: 03-Mar-23 (R2022a)
% Add paths:
% --- Current folder and subfolders
addpath(genpath(pwd))

%%
% Get path to raw data
SET.path2animal = uigetdir(pwd,'Select folder listing all animals');

% Set which delimiter to use to read table
SET.deli = '';

% Get overview of animal folders
curr.dir.all = dir(SET.path2animal);

% Load list of invalid trials
invalidTrials = readtable([SET.path2animal, '\invalidTrials.txt']);

% Iterate over all animals
for iAni = 1:size(curr.dir.all,1)

    % Check whether this is what we are looking for, i.e. whether the
    % string "Animal" is present
    if ~isempty(strfind(curr.dir.all(iAni).name, 'Animal')) && curr.dir.all(iAni).isdir

        % Get information about trials
        curr.path.animal = [SET.path2animal,'\',curr.dir.all(iAni).name,'\01_Data_raw'];
        curr.dir.animal = dir(curr.path.animal);
        idx_start = strfind(curr.dir.all(iAni).name, 'Animal')+length('Animal');
        idx_stop = strfind(curr.dir.all(iAni).name, '_');
        idx_stop = idx_stop(find(idx_stop > idx_start,1))-1;
        curr.animal = curr.dir.all(iAni).name(idx_start:idx_stop);

        % Load meta
        curr.path.data = [SET.path2animal,'\',curr.dir.all(iAni).name,'\02_Data_Matlab\'];
        load([curr.path.data,'meta_info.mat']);

        % Check whether there are entries for the current animal in the
        % list of invalid trials
        idx = find(strcmp(invalidTrials.Animal, ['A', curr.animal]));
        % The same animal could have been tested in different depth (e.g.
        % for confocal microscopy experiments). Check whether this is given
        if any(strcmp('Depth',invalidTrials.Properties.VariableNames))
            curr.depth  = curr.dir.all(iAni).name(strfind(curr.dir.all(iAni).name, '_d')+2 : strfind(curr.dir.all(iAni).name, '_d')+4);
            idx = intersect(idx, find(strcmp(invalidTrials.Depth, ['d', curr.depth])));
        end

        % Iterate over invalid trials and kick them out
        for iTrial = 1:length(idx)
            % Get trial
            T = invalidTrials.Trial(idx(iTrial));
            % kick out of the list of valid trials
            meta.valid_trials(meta.valid_trials == T) = [];
        end%iTrial
        
        % Save meta
        save([curr.path.data,'meta_info.mat'], 'meta');

    end%if animal
end%iAni