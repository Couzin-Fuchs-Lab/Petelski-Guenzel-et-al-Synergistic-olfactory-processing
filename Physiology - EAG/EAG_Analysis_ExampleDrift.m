%% EAG_Analysis_ExampleDrift
% bla
%
% Version: 25-Mar-24 (R2023a)

% Clean up
clc; clear all; close all;
% Add paths:
% --- Current paths
addpath(genpath(pwd))
% --- Toolbox for exporting publication quality figures
%     https://github.com/altmany/export_fig
addpath(genpath('...\GitHub\export_fig'))
mkdir('FIGs')

%% Settings

SET.path2data = "...\EAG\gregarious\example_drift_240329_Animal24_ant02_R_ZHAECOL_greg.mat";
% --- infos on the recording settings
meta_info.channel = 1;
meta_info.event_channel = 2;
meta_info.freq = 25e3;
meta_info.gain = 2;
% --- High-pass filter cutoff freq
SET.hpf = 3; % [Hz]

%% Process data

% Get data and corresponding time vector
data = load(SET.path2data);
wave_data = data.Ch1.values;
events = data.Ch2.times;
time_vector = 0 : 1/meta_info.freq : (length(wave_data)/meta_info.freq)-(1/meta_info.freq);
time_vector = time_vector - events(1);
event_vector = zeros(1,length(time_vector));
event_vector(time_vector>0 & time_vector<=2) = 1;
event_vector(time_vector>10 & time_vector<=12) = 1;
event_vector(time_vector>20 & time_vector<=22) = 1;

% Convert to mV
wave_data = (wave_data/meta_info.gain)*1000;

% Correct offset
[~,ind] = min(abs(time_vector-(-10)));
wave_data = wave_data-wave_data(ind);

% Design a high-pass filter with Bessel-like
% characteristics
d = designfilt('highpassiir', 'FilterOrder', 1, ...
    'HalfPowerFrequency', SET.hpf, ...
    'DesignMethod', 'butter', ...
    'SampleRate', meta_info.freq);
% Extract filter coefficients
[b, a] = tf(d);
% Filter the signal
wave_data_filter = filter(b, a, wave_data);
wave_data_filter = EAG_SubFcn.tc_lowpass(wave_data_filter, 1/(100*(1/meta_info.freq)));

% Check response magnitudes of fitlered and unfiltered data
[~,ind_1] = min(abs(time_vector-0));
[~,ind_2] = min(abs(time_vector-10));
[~,ind_3] = min(abs(time_vector-20));

% Get response magnitudes
mag_raw = zeros(1,3);
mag_raw(1) = wave_data(ind_1) - min(wave_data(time_vector>0 & time_vector<=2));
mag_raw(2) = wave_data(ind_2) - min(wave_data(time_vector>10 & time_vector<=12));
mag_raw(3) = wave_data(ind_3) - min(wave_data(time_vector>20 & time_vector<=22));
mag_filter= zeros(1,3);
mag_filter(1) = wave_data_filter(ind_1) - min(wave_data_filter(time_vector>0 & time_vector<=2));
mag_filter(2) = wave_data_filter(ind_2) - min(wave_data_filter(time_vector>10 & time_vector<=12));
mag_filter(3) = wave_data_filter(ind_3) - min(wave_data_filter(time_vector>20 & time_vector<=22));

% Plot everything
hFig = figure('Color','w','Units','normalized','Position',[0 0 1 0.5]);
subplot(1,3,1); hold on
plot(time_vector, wave_data, 'k')
xlim([-10, 32])
axis square
xlabel('time (s)')
ylabel('response magnitude (mV)')

subplot(1,3,2); hold on
plot(time_vector, wave_data_filter, 'k')
xlim([-10, 32])
axis square
xlabel('time (s)')
ylabel('response magnitude')

subplot(1,3,3); hold on
plot(time_vector, event_vector, 'r')
xlim([-10, 32])
axis square
xlabel('time (s)')
ylabel('stimulus')

% Indicate baseline
[xData, yData] = prepareCurveData(time_vector(time_vector>-10 & time_vector<=0), wave_data(time_vector>-10 & time_vector<=0) );
ft = fittype( 'poly1' );
ft = fit( xData, yData, ft );
trend = time_vector;
trend = ft.p1*trend + ft.p2;

subplot(1,3,1); hold on
plot(time_vector, trend, 'r')
plot([00 02], [1 1]*min(wave_data(time_vector>00 & time_vector<=02)), 'r')
plot([10 12], [1 1]*min(wave_data(time_vector>10 & time_vector<=12)), 'r')
plot([20 22], [1 1]*min(wave_data(time_vector>20 & time_vector<=22)), 'r')
plot([00 00],[wave_data(ind_1), min(wave_data(time_vector>00 & time_vector<=02))],'r')
plot([10 10],[wave_data(ind_2), min(wave_data(time_vector>10 & time_vector<=12))],'r')
plot([20 20],[wave_data(ind_3), min(wave_data(time_vector>20 & time_vector<=22))],'r')

subplot(1,3,2); hold on
plot([-10 32],[0 0], 'r')
plot([00 02], [1 1]*min(wave_data_filter(time_vector>00 & time_vector<=02)), 'r')
plot([10 12], [1 1]*min(wave_data_filter(time_vector>10 & time_vector<=12)), 'r')
plot([20 22], [1 1]*min(wave_data_filter(time_vector>20 & time_vector<=22)), 'r')
plot([00 00],[wave_data_filter(ind_1), min(wave_data_filter(time_vector>00 & time_vector<=02))],'r')
plot([10 10],[wave_data_filter(ind_2), min(wave_data_filter(time_vector>10 & time_vector<=12))],'r')
plot([20 20],[wave_data_filter(ind_3), min(wave_data_filter(time_vector>20 & time_vector<=22))],'r')

% Save
export_fig('FIGs\EAG_Analysis_ExampleDrift', '-pdf', '-painters')
close(hFig)
