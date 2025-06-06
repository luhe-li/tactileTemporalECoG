
projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile'; 
subject = 'ny726';
makePlot = 1;

addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/fileio/'))
addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/utilities/'))
addpath(genpath(fullfile(pwd, 'func')))

%% PREPROCESSING

% Define paths and BIDS specs %%
patientID   = 726; % Specify patient's raw folder name here
RawDataDir  = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/tactile/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs: assuming defaults for a first session, full visual set:
projectName = 'tactile';
sub_label   = 'ny726'; % Specify patient's code name here;
ses_label   = 'nyuecog01';
ses_labelt1 = 'som3t01';
acq_label   = 'clinical';
task_label  = {'temporalpattern', ...            
               'temporalpattern', ... 
               'temporalpattern', ... 
               'temporalpattern', ... 
               'temporalpattern', ... 
              };              
run_label = {'01','02','03','04','05'};
% NOTE: task and run labels should be noted in the order they were run!

% Define paths
[dataReadDir, dataWriteDir, stimWriteDir, T1WriteDir, preprocDir] = bidsconvert_getpaths(patientID, RawDataDir, ...
    BIDSDataDir, projectName, sub_label, ses_label, ses_labelt1);

% Read ECoG data
[rawdata, hdr] = bidsconvert_readecogdata(dataReadDir, ses_label);

% Load channels from the first session (shared across sessions)
[~, channels, ~, ~, ~] = bidsEcogReadFiles(projectDir, subject, ses_label, task_label{1}, run_label{1});

%%
% Define the trigger channel name (probably a 'DC' channel, see hdr.label).
triggerChannelName = 'DC2';
triggerChannel = find(strcmp(triggerChannelName,hdr.label));
figure;plot(rawdata(triggerChannel,:)); 
title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);

brokenblock_start = 1348280;
brokenblock_end = 1412550;

% Extract broken block data
brokenblock_data = rawdata(:, brokenblock_start:brokenblock_end);

% Plot broken block trigger channel
figure;
plot(brokenblock_data(triggerChannel,:));
title(['Broken Block Trigger Channel: ' hdr.label{triggerChannel}]);
xlabel('Sample Number');
ylabel('Amplitude');

% Save broken block data separately
brokenblock_hdr = hdr;
brokenblock_hdr.nSamples = size(brokenblock_data,2);

% Extract event onset times from broken block using same parameters as in a3_checkNoise
peakOpts.minPeakHeight = 0.8;
peakOpts.minPeakProminence = 0.8;
peakOpts.minPeakDistance = 0.05;

triggers = brokenblock_data(triggerChannel,:);
triggers = triggers / max(triggers);
t = ((0:brokenblock_hdr.nSamples-1)/brokenblock_hdr.Fs);

[~,trigger_onsets, widths] = findpeaks(triggers, brokenblock_hdr.Fs,...
    'MinPeakHeight',peakOpts.minPeakHeight,...
    'MinPeakProminence',peakOpts.minPeakProminence,...
    'MinPeakDistance', peakOpts.minPeakDistance);

[~,trigger_onsets_idx] = findpeaks(triggers,...
    'MinPeakHeight',peakOpts.minPeakHeight,...
    'MinPeakProminence',peakOpts.minPeakProminence,...
    'MinPeakDistance', peakOpts.minPeakDistance);

% Filter triggers by width as in a3_checkNoise
valid_index = widths < 0.02; % each trigger pulse is less than 0.02 s
trigger_onsets = trigger_onsets(valid_index);
trigger_onsets_idx = trigger_onsets_idx(valid_index);

% Plot detected triggers
figure; hold on;
plot(t, triggers);
plot(trigger_onsets, triggers(trigger_onsets_idx), 'r*');
title('Broken Block Triggers with Detected Onsets');
xlabel('Time (s)');
ylabel('Normalized Amplitude');

% Plot last 12 trials in a 3x4 grid
figure('Position', [100 100 1600 900]);

% Calculate samples for epoch
epoch_length = 0.5; % seconds
samples_per_epoch = round(epoch_length * brokenblock_hdr.Fs);
epoch_time = (-10:samples_per_epoch-1)/brokenblock_hdr.Fs;

% Get last 12 trigger indices
last_twelve = trigger_onsets_idx(end-11:end);

% Plot each trial in its own subplot
for trial = 1:12
    subplot(3, 4, trial)
    trigger_sample = last_twelve(trial);
    
    if trigger_sample + samples_per_epoch <= size(brokenblock_data, 2)
        epoch_data = brokenblock_data(triggerChannel, (trigger_sample-10):trigger_sample+samples_per_epoch-1);
        plot(epoch_time, epoch_data, 'k-', 'LineWidth', 1)
        title(sprintf('Trial %d', trial))
        xlabel('Time (s)')
        ylabel('Amplitude')
        grid on
    end
end

sgtitle('Last 12 Trials in Broken Block')

%% do CAR on the broken block data, all trials

% get excluded channels as a3
badchannels = [32    33    81    93   102   103   109   110   111   125   127   128   129   130   131   132   133];
chan_index = ismember(1:size(channels,1), badchannels);
channels(chan_index,:).status = repmat({'bad'}, sum(chan_index), 1);

% do CAR
[data_reref, channels_reref, group_indices, group_names] = ecog_performCAR(brokenblock_data, channels);           

% Select last 3 trials
last_three_onsets = trigger_onsets(end-2:end);
last_three_onsets_idx = trigger_onsets_idx(end-2:end);

% Define epoch window
epoch_t = [-0.2 0.8];
epoch_samples = round(epoch_t * brokenblock_hdr.Fs); % Convert to samples
samples_per_epoch = diff(epoch_samples) + 1;

% Initialize data matrix (time x trial x electrode)
num_timepoints = samples_per_epoch;
num_trials = 3;
num_channels = size(data_reref, 1);
epoched_data = zeros(num_timepoints, num_trials, num_channels);

% Extract epochs
for trial = 1:num_trials
    onset_idx = last_three_onsets_idx(trial);
    start_idx = onset_idx + epoch_samples(1);
    end_idx = onset_idx + epoch_samples(2);
    
    if start_idx > 0 && end_idx <= size(data_reref, 2)
        epoched_data(:,trial,:) = data_reref(:,start_idx:end_idx)';
    end
end

% Create time vector for epochs
t = linspace(epoch_t(1), epoch_t(2), num_timepoints);

% group electrode by regions
regions = {'H','M','P','R','S','V','W','Y','Z'};

% Calculate subplot layout - 3x3 grid
subplot_dims = [3 3];

% Create a figure for each trial
for trial = 1:num_trials
    % Create full screen figure
    figure('Position', get(0, 'ScreenSize')); 
    
    % Plot each region
    for r = 1:length(regions)
        % select electrodes for this region
        chanidx = find(cellfun(@(x) ~isempty(x), strfind(brokenblock_hdr.label, regions{r})));
        
        if ~isempty(chanidx)
            subplot(subplot_dims(1), subplot_dims(2), r)
            hold on
            
            % Plot each electrode's response for this trial
            plot(t, squeeze(epoched_data(:,trial,chanidx)), 'LineWidth', 1);
            
            xline(0); yline(0)
            grid on
            
            xlabel('Time (s)')
            ylabel('Amplitude (μV)')
            title(['Region ' regions{r}])
            set(gca, 'XTick', min(epoch_t):0.2:max(epoch_t))
            set(gca, 'FontSize', 10)
        end
    end
    
    sgtitle(sprintf('Trial %d Response by Region', trial), 'FontSize', 14)
end
%% Plot power spectrum for last 3 trials

% Create full screen figure for power spectra
figure('Position', get(0, 'ScreenSize'));

% Calculate frequency parameters
srate  = hdr.Fs;
nfft = 2^nextpow2(srate); % Length of FFT
freq = linspace(0, srate/2, nfft/2+1); % Frequency vector
freq_idx = freq <= 140; % Index for frequencies up to 140 Hz

% Plot power spectrum for each region
for r = 1:length(regions)
    % Select electrodes for this region
    chanidx = find(cellfun(@(x) ~isempty(x), strfind(brokenblock_hdr.label, regions{r})));
    
    if ~isempty(chanidx)
        subplot(subplot_dims(1), subplot_dims(2), r)
        hold on
        
        % Calculate and plot power spectrum for each electrode
        for ch = 1:length(chanidx)
            % Get data for this electrode across last 3 trials
            signal = squeeze(epoched_data(:,end-2:end,chanidx(ch)));
            signal = signal(:)'; % Concatenate trials
            
            % Calculate power spectrum using Welch's method
            [pxx,~] = pwelch(signal, hanning(srate), [], nfft, srate);
            
            % Plot up to 140 Hz
            plot(freq(freq_idx), 10*log10(pxx(freq_idx)), 'LineWidth', 1);
        end
        
        grid on
        xlabel('Frequency (Hz)')
        ylabel('Power (dB)')
        title(['Region ' regions{r}])
        set(gca, 'XTick', 0:20:140)
        set(gca, 'FontSize', 10)
    end
end

sgtitle('Power Spectrum by Region (Last 3 Trials)', 'FontSize', 14)
