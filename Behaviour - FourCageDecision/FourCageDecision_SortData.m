%% SortData_Behaviour
%  This script sorts behavioural data according to the different modalities
%  investigated:
%       'VpOpFp',... visual cues (see-through),   feeding
%       'VpOpFm',... visual cues (see-through),   not feeding
%       'VmOpFp',... no visual cues(opaque),      feeding
%       'VmOpFm',... no visual cues(opaque),      not feeding
%       'VpOmFp',... no olfactory cues (glass),   feeding
%       'VpOmFm',... no olfactory cues (glass),   not feeding
%  Version:
%  10-Feb-2021 (R2020a) Yannick Günzel

clear all; close all; clc
warning('off')

%% Conditions + Abbreviations
% --- Phase
SET.Phase_Abbreviation = {...
    'Gr',... gregarious
    'So'... solitarious
    };
SET.Phase = {'gregarious', 'solitarious'};
% --- Shelter
SET.ShelterCondition_Abbreviation = {...
    'S',... see-through
    'N',... opaque
    'G',... glass sealed
    'GS'... glass sealed + see-through shelter
    };
% --- Feeding
SET.FeedingConditions_Abbreviation = {...
    'F',...  feeding
    'NF',... not feeding
    };

%% Target folders
% Set different combinations of conditions (excluding Gr and So)
SET.Conditions_cages = {...
    'VpOpFp',... visual cues (see-through),   feeding
    'VpOpFm',... visual cues (see-through),   not feeding
    'VmOpFp',... no visual cues(opaque),      feeding
    'VmOpFm',... no visual cues(opaque),      not feeding
    'VpOmFp',... no olfactory cues (glass),   feeding
    'VpOmFm',... no olfactory cues (glass),   not feeding
    };
SET.Conditions_cotton = {...
    'LoLeavesCo',... 
    'LoLeavesZ',... 
    };
% Set paths
SET.BasePath = '...\Data\Behaviour\';


% Create folders
for iPhase = 1:length(SET.Phase)
    for iCondition = 1:length(SET.Conditions_cages)
        % Experiemnts with 4 cages
        mkdir([SET.BasePath, SET.Phase{iPhase},'\cages\',SET.Conditions_cages{iCondition}])        
    end%iCondition
    for iCondition = 1:length(SET.Conditions_cotton)
        % Experiemnts with 4 cages
        mkdir([SET.BasePath, SET.Phase{iPhase},'\cotton\',SET.Conditions_cotton{iCondition}])        
    end%iCondition
end%iPhase

%% Sort data

% -------------------------------------------------------------------------
%                      Iterate over content of the folder "unsorted_cages"
% -------------------------------------------------------------------------
% Waitbar
hWait = waitbar(0,'Please wait...');
% Get content
Content = dir([SET.BasePath,'unsorted_cages\']);
% Iterate over files
for iFile = 1:size(Content,1)
    % Waitbar
    hWait = waitbar(iFile/size(Content,1),hWait);
    % Check whether it is a file
    if ~Content(iFile).isdir
        % Get current file
        currFile = Content(iFile).name;
        % Find "_" in string to navigate
        idx = strfind(currFile, '_');
        % Get feeding condition
        cond.Feeding = currFile(idx(3)+1:idx(4)-1);
        % Get cage condition
        cond.Cage = currFile(idx(9)+1:idx(10)-1);
        % Get phase
        cond.Phase = currFile(idx(10)+1:idx(11)-1);
        
        % Differentiate between conditions and copy files accordingly
        if strcmp(cond.Phase, 'Gr')
            if strcmp(cond.Cage, 'S') && strcmp(cond.Feeding, 'F')%                                     'VpOpFp'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'gregarious\cages\','VpOpFp'])
            elseif strcmp(cond.Cage, 'S') && strcmp(cond.Feeding, 'NF')%                                'VpOpFm'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'gregarious\cages\','VpOpFm'])
            elseif strcmp(cond.Cage, 'N') && strcmp(cond.Feeding, 'F')%                                 'VmOpFp'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'gregarious\cages\','VmOpFp'])
            elseif strcmp(cond.Cage, 'N') && strcmp(cond.Feeding, 'NF')%                                'VmOpFm'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'gregarious\cages\','VmOpFm'])
            elseif (strcmp(cond.Cage, 'G') || strcmp(cond.Cage, 'GS')) && strcmp(cond.Feeding, 'F')%    'VpOmFp'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'gregarious\cages\','VpOmFp'])
            elseif (strcmp(cond.Cage, 'G') || strcmp(cond.Cage, 'GS')) && strcmp(cond.Feeding, 'NF')%	'VpOmFm'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'gregarious\cages\','VpOmFm'])
            end%if condition
        else% ---------------- solitarious animals ----------------
            if strcmp(cond.Cage, 'S') && strcmp(cond.Feeding, 'F')%                                     'VpOpFp'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'solitarious\cages\','VpOpFp'])
            elseif strcmp(cond.Cage, 'S') && strcmp(cond.Feeding, 'NF')%                                'VpOpFm'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'solitarious\cages\','VpOpFm'])
            elseif strcmp(cond.Cage, 'N') && strcmp(cond.Feeding, 'F')%                                 'VmOpFp'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'solitarious\cages\','VmOpFp'])
            elseif strcmp(cond.Cage, 'N') && strcmp(cond.Feeding, 'NF')%                                'VmOpFm'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'solitarious\cages\','VmOpFm'])
            elseif (strcmp(cond.Cage, 'G') || strcmp(cond.Cage, 'GS')) && strcmp(cond.Feeding, 'F')%    'VpOmFp'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'solitarious\cages\','VpOmFp'])
            elseif (strcmp(cond.Cage, 'G') || strcmp(cond.Cage, 'GS')) && strcmp(cond.Feeding, 'NF')%	'VpOmFm'
                copyfile(...
                    [SET.BasePath,'unsorted_cages\', currFile],...
                    [SET.BasePath,'solitarious\cages\','VpOmFm'])
            end%if condition
        end%if phase        
    end%if
end%iFile
% Close waitbar
close(hWait)
clear cond Content currFile hWait iCondition idx iFile iPhase



% -------------------------------------------------------------------------
%                      Iterate over content of the folder "unsorted_cotton"
% -------------------------------------------------------------------------
% Waitbar
hWait = waitbar(0,'Please wait...');
% Get content
Content = dir([SET.BasePath,'unsorted_cotton\']);
% Iterate over files
for iFile = 1:size(Content,1)
    % Waitbar
    hWait = waitbar(iFile/size(Content,1),hWait);
    % Check whether it is a file
    if ~Content(iFile).isdir
        % Get current file
        currFile = Content(iFile).name;
        % Find "_" in string to navigate
        idx = strfind(currFile, '_');
        % Get feeding condition
        cond.Feeding = currFile(idx(3)+1:idx(4)-1);
        % Get phase
        cond.Phase = currFile(idx(8)+1:idx(9)-1);
        
        % Differentiate between conditions and copy files accordingly
        if strcmp(cond.Phase, 'Gr')
            if ~isempty(strfind(currFile, SET.Conditions_cotton{1}))%                                     'LoLeavesCo'
                copyfile(...
                    [SET.BasePath,'unsorted_cotton\', currFile],...
                    [SET.BasePath,'gregarious\cotton\',SET.Conditions_cotton{1}])
            elseif ~isempty(strfind(currFile, SET.Conditions_cotton{2}))%                                'LoLeavesZ'
                copyfile(...
                    [SET.BasePath,'unsorted_cotton\', currFile],...
                    [SET.BasePath,'gregarious\cotton\',SET.Conditions_cotton{2}])           
            end%if condition
        else% ---------------- solitarious animals ----------------
            if ~isempty(strfind(currFile, SET.Conditions_cotton{1}))%                                     'LoLeavesCo'
                copyfile(...
                    [SET.BasePath,'unsorted_cotton\', currFile],...
                    [SET.BasePath,'gregarious\cotton\',SET.Conditions_cotton{1}])
            elseif ~isempty(strfind(currFile, SET.Conditions_cotton{2}))%                                'LoLeavesZ'
                copyfile(...
                    [SET.BasePath,'unsorted_cotton\', currFile],...
                    [SET.BasePath,'gregarious\cotton\',SET.Conditions_cotton{2}])           
            end%if condition
        end%if phase        
    end%if
end%iFile
% Close waitbar
close(hWait)





































