
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

% Select last 5 trials
last_ten_onsets = trigger_onsets(end-9:end);
last_ten_onsets_idx = trigger_onsets_idx(end-9:end);

% Define epoch window
epoch_t = [-0.4 1.8];
epoch_samples = round(epoch_t * brokenblock_hdr.Fs); % Convert to samples
samples_per_epoch = diff(epoch_samples) + 1;

% Initialize data matrix (time x trial x electrode)
num_timepoints = samples_per_epoch;
num_trials = 5;
num_channels = size(brokenblock_data, 1);
epoched_data = zeros(num_timepoints, num_trials, num_channels);

% Extract epochs
for trial = 1:num_trials
    onset_idx = last_ten_onsets_idx(trial);
    start_idx = onset_idx + epoch_samples(1);
    end_idx = onset_idx + epoch_samples(2);
    
    if start_idx > 0 && end_idx <= size(brokenblock_data, 2)
        epoched_data(:,trial,:) = brokenblock_data(:,start_idx:end_idx)';
    end
end

% Create time vector for epochs
t = linspace(epoch_t(1), epoch_t(2), num_timepoints);

% group electrode by regions
regions = {'V','W','Y','Z','H','M','P','R','S'};
% plot each region separately
for r = 1:length(regions)
    % select electrodes for this region
    chanidx = find(cellfun(@(x) ~isempty(x), strfind(brokenblock_hdr.label, regions{r})));
    
    if ~isempty(chanidx)
        % Plot trial-averaged response for all electrodes in this region
        figure('Position', [100 100 1600 1200]); hold on
        set(gca, 'FontSize', 15)
        plot(t, squeeze(mean(epoched_data(:,:,chanidx), 2)), 'LineWidth', 1);
        xline(0); yline(0)
        set(gca, 'XGrid', 'on')
        set(gca, 'XTick', -0.4:0.1:1.8)
        title(['Region ' regions{r} ' - Trial-averaged response'])
        xlabel('Time (s)')
        ylabel('Amplitude')
    end
end