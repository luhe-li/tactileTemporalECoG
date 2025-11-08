
clear;
% tbUse tactileTemporalECoG

projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
subject = 'ny726';

sessions = 'nyuecog01';
outputFolder = 'demean_CAR';
tasks = 'tacttestascending';
runnums = [];

% addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/fileio/'))
% addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/utilities/'))

% %% apply demean CAR (do it one)

% bidsEcogRereference(projectDir, subject, sessions, tasks, runnums, outputFolder)

%% Load data

inputFolder = 'demean_CAR';
description = 'reref';

% Load data
dataPath = fullfile(projectDir, 'derivatives', inputFolder);
[data, channels, events, srate] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description);

%% Phase analysis

% Number of repeats across this run
n_cycle = 5;
n_finger = 5;

% Run start and end time
run_onset = events.onset(1);
trial_duration = events.duration(1) + events.ISI(1);
run_offset = events.onset(end) + trial_duration;

% Duration for one sweep across all fingers
t_sweep = trial_duration * n_finger;
f0 = 1/t_sweep;

% Sample index
ix = round(run_onset*srate)+1:round(run_offset*srate);
run_data = data(:, ix);
n_samples = size(run_data, 2);
t = (0:n_samples-1)/srate;

% Detrend each channel
detrend_data = detrend(run_data', 'linear')';

% % FFT
% ft = fft(detrend_data,[],2);
% ft = ft(:, 1:floor(n_samples/2)+1);
% freqs = (0:floor(n_samples/2)) * srate / n_samples;
% 
% % Index of target frequency
% [~, f0_idx] = min(abs(freqs - f0));
% amp_f0 = abs(ft(:, f0_idx)) * 2/n_samples;
% 
% % Total magnitude across all frequencies
% total_mag = sqrt(sum(abs(ft).^2, 2));
% 
% % Coherence is normalized amplitude
% coh = amp_f0 ./ total_mag;
% [coh_sorted, ch_idx] = sort(coh, 'descend');
% 
% % return top responsive electrodes
% topN = 10;
% responsive_channels = ch_idx(1:topN);


% Window
w = hann(n_samples)';          
X = fft(detrend_data .* w, [], 2);
X = X(:, 1:floor(n_samples/2)+1);

% Frequency axis
freqs = (0:floor(n_samples/2)) * srate / n_samples;

% Index of f0
[~, f0_idx] = min(abs(freqs - f0));

% Window-corrected amplitude scaling
win_scale = 2 / (sum(w)/n_samples);   % correcting for Hann attenuation
amp_f0 = abs(X(:, f0_idx)) * win_scale;

% Total magnitude across spectrum
total_mag = sqrt(sum(abs(X).^2, 2));

% Coherence
coh = amp_f0 ./ total_mag;
