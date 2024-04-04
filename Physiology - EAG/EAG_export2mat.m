%% EAG_export2mat
% bla
%
% Version: 04-April-24 (R2023a)

% Clean up
clc; clear all; close all;

%% Settings

% --- where can we find the data
SET.path2data = '...\EAG\';
% --- where can we find the CED scripts tp load spike data
SET.path2CEDS64ML = 'C:\CEDMATLAB\CEDS64ML';
% --- where can we find the annotation
SET.path2annotation = 'ANNOTATION.csv';
% --- infos on the recording settings
meta_info.channel = 1;
meta_info.event_channel = 2;
meta_info.freq = 25e3;
% --- High-pass filter cutoff freq
SET.hpf = 3; % [Hz]
% --- window size to correct drift (signal-movmedian(signal))
SET.drift_window = 4;
% --- crop data to before the first and after the third rep [sec]
SET.crop_all = [4 4];
% --- baseline window [sec]
SET.crop_rep = 2;
SET.stim_dur = 2;

%% CED environment

% Add CED environment
setenv('CEDS64ML', SET.path2CEDS64ML);
cedpath = getenv('CEDS64ML');
addpath(cedpath);
CEDS64LoadLib( cedpath );
clear cedpath

% Load annotation
ANNOTATION = readtable(SET.path2annotation);

%% Iterate over all files
for iFile = 1:size(ANNOTATION,1)

        % Construct path to file
        meta_info.filename = [...
            num2str(ANNOTATION.date(iFile)),'_Animal',sprintf('%02d',ANNOTATION.animal(iFile)),...
            '_ant0',num2str(ANNOTATION.ant(iFile)),'_',ANNOTATION.ant_side{iFile},...
            '_',ANNOTATION.stim{iFile},...
            '_',ANNOTATION.phase{iFile}(1:4)];
        meta_info.path2file = [...
            SET.path2data,...
            ANNOTATION.phase{iFile},'\',...
            num2str(ANNOTATION.date(iFile)),'_Animal',sprintf('%02d',ANNOTATION.animal(iFile)),'\',...
            'ant0',num2str(ANNOTATION.ant(iFile)),'_',ANNOTATION.ant_side{iFile},'\',...
            meta_info.filename,'.smrx'];

        % Try to open the smrx file
        fhand1 = CEDS64Open(meta_info.path2file);
        if (fhand1 <= 0)
            disp('- - - - -')
            disp('Cannot open file:')
            disp(meta_info.filename)
            disp(' ')
        else
            % Get the event channel
            % --- Get number of ticks
            maxTimeTicks = CEDS64ChanMaxTime( fhand1, meta_info.event_channel )+1;
            % --- Get data for all ticks
            [ ~, vi64T ] = CEDS64ReadEvents( fhand1, meta_info.event_channel, maxTimeTicks, 0);
            % --- Events are registered as on/off times. Get the corresponding
            %     times in seconds
            TimeVec_tics = CEDS64TicksToSecs(fhand1,1:maxTimeTicks);
            meta_info.event_data = TimeVec_tics(vi64T);
            if length(meta_info.event_data)>6
                meta_info.event_data = meta_info.event_data(end-5:end);
                disp('- - - - -')
                disp('many stim in:')
                disp(meta_info.filename)
                disp(' ')
            end%
            clear vi64T TimeVec_tics
            % Now, get the waveform of the correct channel
            % --- Get number of ticks
            maxTimeTicks = CEDS64ChanMaxTime( fhand1, meta_info.channel )+1;
            % --- Get data for all ticks
            [~, fVals, ~] = CEDS64ReadWaveF( fhand1, meta_info.channel, maxTimeTicks, 0 );
            % --- Construct the time vector
            data.time_data = linspace(1/meta_info.freq, length(fVals)/meta_info.freq, length(fVals));
            % --- Keep the data with double precision
            data.wave_data = double(fVals);
            clear fVals fTime maxTimeTicks

            % Correct polarity
            if strcmp(ANNOTATION.polarity{iFile}, 'pos')
                data.wave_data = -data.wave_data;
            end

            % Convert to mV
            data.wave_data = (data.wave_data/ANNOTATION.gain(iFile))*1000;

            % Generally crop data
            [~, ind1] = min(abs(data.time_data - (meta_info.event_data(1) - SET.crop_all(1) )));
            [~, ind2] = min(abs(data.time_data - (meta_info.event_data(end) + SET.crop_all(2) )));
            data.time_data = data.time_data(ind1:ind2);
            data.wave_data = data.wave_data(ind1:ind2);
            % Adjust event and time
            meta_info.event_data = meta_info.event_data - data.time_data(1);
            data.time_data = data.time_data - data.time_data(1);
            % Reshape events
            meta_info.event_data = [meta_info.event_data(1:2:end); meta_info.event_data(2:2:end)];

            % Crop out responses to each repetition
            for iRep = 1:3
                % Get the data
                [~, ind1] = min(abs(data.time_data - (meta_info.event_data(1,iRep) - SET.crop_rep )));
                ind2 = ind1 - 1 + (2*SET.crop_rep+SET.stim_dur)*meta_info.freq;
                temp = data.wave_data(ind1:ind2);
                % Get time vector
                data.time_rep = ((0:length(temp)-1)/meta_info.freq)-SET.crop_rep;

                % Correct baseline drift
                BL = temp(1:(SET.crop_rep*meta_info.freq));
                [xData, yData] = prepareCurveData( linspace(-SET.crop_rep,0,length(BL)), BL(:)' );
                ft = fittype( 'poly1' );
                ft = fit( xData, yData, ft );
                trend = data.time_rep;
                trend = ft.p1*trend + ft.p2;
                temp = temp(:)-trend(:);
                % Correct baseline offset
                BL = temp(1:(SET.crop_rep*meta_info.freq));
                temp = temp - median(BL);

                % Filter data
                if ANNOTATION.AC(iFile) == 0
                    temp_original = temp;
                    % Design a high-pass filter with Bessel-like
                    % characteristics
                    d = designfilt('highpassiir', 'FilterOrder', 1, ...
                        'HalfPowerFrequency', SET.hpf, ...
                        'DesignMethod', 'butter', ...
                        'SampleRate', meta_info.freq);
                    % Extract filter coefficients
                    [b, a] = tf(d);
                    % Filter the signal
                    temp = filter(b, a, temp);
                    clear a b d
                elseif ANNOTATION.AC(iFile) == 3
                    temp_original = cumsum(temp) * (SET.hpf*(1/meta_info.freq));
                end% if offline filtering
                temp = EAG_SubFcn.tc_lowpass(temp, 1/(100*(1/meta_info.freq)));

                % Store
                data.(['wave_data_original_rep',num2str(iRep)]) = temp_original;
                data.(['wave_data_rep',num2str(iRep)]) = temp;

                clearvars -except SET ANNOTATION data meta_info iFile iRep fhand1
            end%iRep
            data.time_rep = ((0:length(data.(['wave_data_rep',num2str(iRep)]))-1)/meta_info.freq)-SET.crop_rep;
            data.time_rep = data.time_rep';

            % Save all
            save([SET.path2data,...
                ANNOTATION.phase{iFile},'\',...
                num2str(ANNOTATION.date(iFile)),'_Animal',sprintf('%02d',ANNOTATION.animal(iFile)),'\',...
                'ant0',num2str(ANNOTATION.ant(iFile)),'_',ANNOTATION.ant_side{iFile},'\',...
                meta_info.filename,'.mat'],...
                'data', 'meta_info')

        end%if file
        CEDS64Close(fhand1);
end%iFile