clear

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

%% check visual triggger. Exclude the broken-off 4th run.

% Define the trigger channel name (probably a 'DC' channel, see hdr.label).
triggerChannelName = 'DC2';
triggerChannel = find(strcmp(triggerChannelName,hdr.label));
figure;plot(rawdata(triggerChannel,:)); 
title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);

run_start = 733962; % Manually determined from plot of trigger channel 
t2 = 1319100;
t3 = 1414680;
run_end   = 1993730; 

% clip the data
clip1 = rawdata(:,run_start:t2);
clip2 = rawdata(:,t3:run_end);
data = [clip1, clip2];
hdr.nSamples = size(data,2);
% 
% % check the data
% figure;plot(data(triggerChannel,:)/max(data(triggerChannel,:))); 
% title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);

% extract event onset times
peakOpts.minPeakHeight = 0.8;
peakOpts.minPeakProminence = 0.8;
peakOpts.minPeakDistance = 0.05;

triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);
t = ((0:hdr.nSamples-1)/hdr.Fs);

[~,trigger_onsets, widths] = findpeaks(triggers, hdr.Fs,...
    'MinPeakHeight',peakOpts.minPeakHeight,...
    'MinPeakProminence',peakOpts.minPeakProminence,...
    'MinPeakDistance', peakOpts.minPeakDistance);

[~,trigger_onsets_idx] = findpeaks(triggers,...
    'MinPeakHeight',peakOpts.minPeakHeight,...
    'MinPeakProminence',peakOpts.minPeakProminence,...
    'MinPeakDistance', peakOpts.minPeakDistance);

% use figure; histogram(widths(valid_index)) to check the threshold of trigger duration
valid_index = widths < 0.02; % each trigger pulse is less than 0.02 s
trigger_onsets = trigger_onsets(valid_index);
trigger_onsets_idx = trigger_onsets_idx(valid_index);

% Segment triggers into blocks using the gap threshold (e.g., 5 seconds)
dt = diff(trigger_onsets);                   % Inter-trigger intervals
block_break_indices = find(dt > 5);            % Find large gaps indicating block boundaries
numBlocks = length(block_break_indices) + 1;   % Number of blocks is one more than the number of breaks

% Split trigger_onsets into blocks
[blocks_idx, blocks]= deal(cell(numBlocks, 1));
startIdx = 1;
for i = 1:length(block_break_indices)
    endIdx = block_break_indices(i);
    blocks{i} = trigger_onsets(startIdx:endIdx);
    blocks_idx{i} = trigger_onsets_idx(startIdx:endIdx);
    startIdx = endIdx + 1;
end
blocks{numBlocks} = trigger_onsets(startIdx:end);
blocks_idx{numBlocks} = trigger_onsets_idx(startIdx:end);

% Create full screen figure
figure('Position', get(0, 'ScreenSize'));
% Create full screen figure
figure('Position', get(0, 'ScreenSize'));

% Calculate subplot layout for all channels
num_channels = size(data,1);
subplot_rows = ceil(sqrt(num_channels));
subplot_cols = ceil(num_channels/subplot_rows);

% Calculate samples for epoch window
pre_trigger = 0.3; % seconds before trigger
post_trigger = 0.3; % seconds after trigger
pre_samples = round(pre_trigger * hdr.Fs);
post_samples = round(post_trigger * hdr.Fs);
epoch_time = linspace(-pre_trigger, post_trigger, pre_samples + post_samples + 1);

% Plot each channel
for ch_idx = 1:num_channels
    subplot(subplot_rows, subplot_cols, ch_idx)
    hold on

    % Get all trigger indices
    all_triggers = trigger_onsets_idx;

    % Plot each trial
    for trial = 1:length(all_triggers)
        trigger_sample = all_triggers(trial);

        % Check if we have enough samples before and after trigger
        if trigger_sample - pre_samples > 0 && trigger_sample + post_samples <= size(data, 2)
            epoch_data = data(ch_idx, ...
                trigger_sample-pre_samples:trigger_sample+post_samples);
            plot(epoch_time, epoch_data, 'k-', 'LineWidth', 0.5, 'Color', [0 0 0 0.2])
        end
    end

    % Add vertical line at trigger onset
    xline(0, '--r', 'LineWidth', 1);
    yline(0, '--k');

    % Labels and formatting
    title(['Channel: ' hdr.label{ch_idx}])
    xlabel('Time (s)')
    ylabel('Amplitude (μV)')
    grid on
    set(gca, 'FontSize', 8)
end

sgtitle('All Channels: Trigger-aligned Responses', 'FontSize', 14)

% Save figure
saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-visual_triggers',sub_label, ses_label)), 'epsc');

% Remove first trigger from each block
for i = 1:numBlocks
    blocks{i} = blocks{i}(2:end);  % Remove first trigger
    blocks_idx{i} = blocks_idx{i}(2:end); % Remove first trigger index
end

% Display block counts after removing first triggers
fprintf('\nTrigger counts per block AFTER removing first triggers:\n');
for i = 1:numBlocks
    fprintf('Block %d: %d triggers\n', i, numel(blocks{i}));
end

% Recombine blocks into single vector
trigger_onsets = horzcat(blocks{:});
trigger_onsets_idx = horzcat(blocks_idx{:});
fprintf('\nTotal triggers after removing first triggers: %d\n', numel(trigger_onsets));


%% channel selection

% Define time axis (in seconds). First time point = 0 (this is assumed by
% the function we used to detect triggers below, and also in fieldtrip).
t = ((0:hdr.nSamples-1)/hdr.Fs); 

% Automatic identification
% Identify channels with values exceeding the threshold
threshold = 200;
exclude_inx = [];
for cChan = 1:size(data, 1)
    if mean(data(cChan, :)) > threshold || mean(data(cChan, :)) <- threshold || max(data(cChan,:))>1500
        exclude_inx = [exclude_inx, cChan];
    end
end

% Manually add bad channels
manual_exclude_inx = [32,33,81,93,109,125];
exclude_inx = sort(unique([exclude_inx, manual_exclude_inx]));

% % Check bad channels
% for cChan = 1:1:numel(exclude_inx)
%     exc_ch = exclude_inx(cChan);
%     figure;plot(t,data(exc_ch,:));
%     title([num2str(exc_ch) ': ' hdr.label{exc_ch}]);
%     xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); set(gca,'fontsize',16);
% end

% Trigger channel name (probably a 'DC' channel, see hdr.label)
triggerChannelName = 'DC2'; % visual trigger

% Specify reasons for marked as bad, e.g. spikes, elipeptic,
% outlierspectrum, lowfreqdrift
BADCHANNELS_MANUALTABLE = [num2cell(exclude_inx)' repmat({'spikes'}, [length(exclude_inx) 1])];
badChannels = cell2mat(BADCHANNELS_MANUALTABLE(:,1));
badChannelsDescriptions = BADCHANNELS_MANUALTABLE(:,2);

% Generate electrode files
[electrode_table, channel_table] = bidsconvert_getelectrodefiles(dataReadDir, hdr, triggerChannel, badChannels, badChannelsDescriptions);

%% read in stimulus files

% READ IN relevant data files %%%%%%%%%%%%%%%%%%

%   - stimulus log files generated by stimulus code
stimDir = fullfile(dataReadDir, 'stimdata');
stimMatFiles = [dir(fullfile(stimDir, sprintf('sub-*%d*%s*.mat', patientID, ses_label))) dir(fullfile(stimDir, sprintf('sub-*%d*%s*.mat', patientID, upper(ses_label))))];
stimTSVFiles = [dir(fullfile(stimDir, sprintf('sub-*%d*%s*.tsv', patientID, ses_label))) dir(fullfile(stimDir, sprintf('sub-*%d*%s*.tsv', patientID, upper(ses_label))))];

% CHECK: Do we have all the stimfiles?
fprintf('[%s] Asserting that we have all the stimFiles \n', mfilename);
assert(isequal(length(stimMatFiles), length(run_label)))
assert(isequal(length(stimTSVFiles), length(run_label)))

% NOTE: It's important to read the StimFiles in in the correct order.
% Because the logfiles are copied over, we can't assume that
% stimFiles(ii).date is correct. Instead, we sort them based on the field
% "experimentDateandTime" in params saved out by the stimulus code.

% Read the stimfiles
clear stimData;
fprintf('[%s] Reading in the following stimFiles: \n', mfilename);
runTimes = [];
for ii = 1:length(stimMatFiles)
    stimData(ii) = load([stimDir filesep stimMatFiles(ii).name]);
    disp(stimData(ii).fname)
    runTimes{ii} = stimData(ii).params.experimentDateandTime;
end
% Sort the stimfiles, count how many triggers were requested
[~, runIndex] = sort(runTimes);
fprintf('[%s] Sorting runs in recorded order. New run order is: \n', mfilename);
stimData_sorted = stimData(runIndex);
runTimes_sorted = runTimes(runIndex);
requestedTriggerCount = 0;
for ii = 1:length(stimMatFiles)
    disp([stimData_sorted(ii).fname])
    trigSeq = stimData_sorted(ii).stimulus.trigSeq;
    % Remove first and last trigger (block start and end)
    trigSeq = trigSeq(2:end-1);
    num_triggers = length(find(trigSeq));
    requestedTriggerCount = requestedTriggerCount+num_triggers;
end

% CHECK: Does the number of requested triggers match the number of triggers
% that were detected in the trigger channel?
fprintf('[%s] Asserting that number of triggers match requested number \n', mfilename);
foundTriggerCount = length(trigger_onsets);
if isequal(requestedTriggerCount, foundTriggerCount)
    triggersAreMatched = 1;
else
    triggersAreMatched = 0;
    warning('[%s] Number of triggers does not match requested number!!! \n Triggers in data = %d, Triggers requested = %d. \n Please double check the cause of this discrepancy before continuing.', mfilename, foundTriggerCount, requestedTriggerCount);
end

stimData = stimData_sorted;
runTimes = runTimes_sorted;

if makePlot
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-triggers_requested',sub_label, ses_label)), 'epsc');
end

%% writing files (do it once)

% Write run files
[dataFileNames] = bidsconvert_writerunfiles(dataWriteDir, stimWriteDir, ...
    sub_label, ses_label, task_label, acq_label, run_label, ...
    data, hdr, stimData, channel_table, trigger_onsets);

% Write session files
bidsconvert_writesessionfiles(dataReadDir, dataWriteDir, T1WriteDir, ...
    sub_label, ses_label, acq_label, ses_labelt1, electrode_table, dataFileNames, runTimes);

%% common average reference (do it once)

bidsEcogRereference(projectDir, subject);

%% CHECK CAR data

bidsRootPath = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
dataPath = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR');
subject           = 'ny726';
session           = 'nyuecog01';
task              = 'temporalpattern';

[data, channels, events, srate] = bidsEcogGetPreprocData(dataPath, subject, [], task);

align_onset_epoch_t     = [-0.2 0.8];  % stimulus epoch window
align_offset_epoch_t     = [-0.5 0.5];  % stimulus epoch window

% group electrode by regions
regions = {'H','M','P','R','S','V','W','Y','Z'};

% epoch: time series x trial x electrode
[epochs, t] = ecog_makeEpochs(data, events.onset, align_onset_epoch_t, channels.sampling_frequency(1));

%% plot: align trials by trial onset

% Calculate subplot layout
num_regions = length(regions);
num_rows = ceil(sqrt(num_regions));
num_cols = ceil(num_regions/num_rows);

% Create full screen figure for trial-averaged responses
figure('Position', get(0, 'ScreenSize')); 

for r = 1:length(regions)
    % select electrodes for this region
    chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
    
    if ~isempty(chanidx)
        subplot(num_rows, num_cols, r)
        hold on
        plot(t, squeeze(mean(epochs(:,:,chanidx), 2)),'LineWidth',1);
        xline(0); yline(0)
        grid on
        title(['Region ' regions{r}])
        
        % Only show x label on bottom row
        if r > (num_rows-1)*num_cols
            xlabel('Time (s)')
        end
        
        % Only show y label on first column
        if mod(r-1, num_cols) == 0
            ylabel('Amplitude (μV)')
        end
        
        set(gca, 'XTick', min(align_onset_epoch_t):0.2:max(align_onset_epoch_t))
        set(gca, 'FontSize', 10)
    end
end
sgtitle('Trial-Averaged Response by Region', 'FontSize', 14)

% Save figure
figureDir = fullfile(bidsRootPath, 'derivatives', 'ECoGFigures', ['sub-' subject], 'CAR_figures');
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end
saveas(gcf, fullfile(figureDir, 'all_regions_avg.png'))

%% align trials by trial offset

% Calculate event offsets based on trial type
offsets = zeros(size(events.onset));
for i = 1:length(events.trial_name)
    if contains(events.trial_name{i}, 'TWO-PULSE')
        % For two pulse trials, offset = onset + 2*duration + ISI 
        offsets(i) = events.onset(i) + 2*events.duration(i) + events.ISI(i);
    else
        % For one pulse trials, offset = onset + duration
        offsets(i) = events.onset(i) + events.duration(i);
    end
end

% Make epochs aligned to offset
[epochs_offset, t_offset] = ecog_makeEpochs(data, offsets, align_offset_epoch_t, channels.sampling_frequency(1));

% Create full screen figure for offset-aligned responses
figure('Position', get(0, 'ScreenSize')); 

for r = 1:length(regions)
    % select electrodes for this region
    chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
    
    if ~isempty(chanidx)
        subplot(num_rows, num_cols, r)
        hold on
        plot(t_offset, squeeze(mean(epochs_offset(:,:,chanidx), 2)),'LineWidth',1);
        xline(0); yline(0)
        grid on
        title(['Region ' regions{r}])
        
        % Only show x label on bottom row
        if r > (num_rows-1)*num_cols
            xlabel('Time from offset (s)')
        end
        
        % Only show y label on first column
        if mod(r-1, num_cols) == 0
            ylabel('Amplitude (μV)')
        end
        
        set(gca, 'XTick', min(align_offset_epoch_t):0.2:max(align_offset_epoch_t))
        set(gca, 'FontSize', 10)
    end
end
sgtitle('Trial-Averaged Response by Region (Offset-Aligned)', 'FontSize', 14)

% Save figure
saveas(gcf, fullfile(figureDir, 'all_regions_offset.png'))

%% check one trial aligned by onset and offset
% Loop through all trials
for current_trial = 1:size(epochs,2)
    % Create subplots for onset-aligned data
    figure('Position', get(0, 'ScreenSize')); 
    
    % Calculate subplot layout - 3x3 grid
    subplot_dims = [3 3];
    
    % Extract epochs for current trial
    selected_epoch = epochs(:,current_trial,:);
    selected_epoch_offset = epochs_offset(:,current_trial,:);
    
    % Plot each region for onset-aligned data
    for r = 1:length(regions)
        % select electrodes for this region
        chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
        
        if ~isempty(chanidx)
            subplot(subplot_dims(1), subplot_dims(2), r)
            hold on
            
            % Plot each electrode's response for this trial
            plot(t, squeeze(selected_epoch(:,1,chanidx)), 'LineWidth', 1);
            xline(0); yline(0)
            grid on
            
            xlabel('Time (s)')
            ylabel('Amplitude (μV)')
            title(['Region ' regions{r}])
            set(gca, 'XTick', min(align_onset_epoch_t):0.2:max(align_onset_epoch_t))
            set(gca, 'FontSize', 10)
        end
    end
    
    sgtitle(sprintf('Trial %d - Onset aligned', current_trial), 'FontSize', 14)
    
    % % Create subplots for offset-aligned data
    % figure('Position', get(0, 'ScreenSize')); 
    % 
    % % Plot each region for offset-aligned data
    % for r = 1:length(regions)
    %     % select electrodes for this region
    %     chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
    % 
    %     if ~isempty(chanidx)
    %         subplot(subplot_dims(1), subplot_dims(2), r)
    %         hold on
    % 
    %         % Plot each electrode's response for this trial
    %         plot(t_offset, squeeze(selected_epoch_offset(:,1,chanidx)), 'LineWidth', 1);
    %         xline(0); yline(0)
    %         grid on
    % 
    %         xlabel('Time from offset (s)')
    %         ylabel('Amplitude (μV)')
    %         title(['Region ' regions{r}])
    %         set(gca, 'XTick', min(align_offset_epoch_t):0.2:max(align_offset_epoch_t))
    %         set(gca, 'FontSize', 10)
    %     end
    % end
    % 
    % sgtitle(sprintf('Trial %d - Offset aligned', current_trial), 'FontSize', 14)
    
    % Wait for mouse click
    w = waitforbuttonpress;
    if w == 0 % Mouse click
        close all;
        continue;
    else % Keyboard press
        close all;
        break;
    end
end


%% check all channels with trial-aligned responses

% Create a figure to show all channels with trial-aligned responses
figure('Position', get(0, 'ScreenSize'));

% Calculate number of channels and subplot layout
num_channels = size(data, 1);
subplot_rows = ceil(sqrt(num_channels));
subplot_cols = ceil(num_channels/subplot_rows);

% Calculate epoch window samples
pre_samples = round(align_onset_epoch_t(1) * srate);
post_samples = round(align_onset_epoch_t(2) * srate);
epoch_time = linspace(align_onset_epoch_t(1), align_onset_epoch_t(2), post_samples - pre_samples + 1);

% Plot each channel
for ch = 1:num_channels
    subplot(subplot_rows, subplot_cols, ch)
    hold on
    
    % Plot each trial for this channel
    for trial = 1:length(events.onset)
        % Get sample index for trial onset
        onset_idx = round(events.onset(trial) * srate);
        
        % Calculate epoch indices
        start_idx = onset_idx + pre_samples;
        end_idx = onset_idx + post_samples;
        
        % Plot if we have enough samples
        if start_idx > 0 && end_idx <= size(data, 2)
            trial_data = data(ch, start_idx:end_idx);
            plot(epoch_time, trial_data, 'k-', 'LineWidth', 0.5, 'Color', [0 0 0 0.2])
        end
    end
    
    % Add vertical line at stimulus onset
    xline(0, '--r', 'LineWidth', 1);
    yline(0, '--k');
    
    % Labels and formatting
    title(['Channel: ' channels.name{ch}])
    xlabel('Time (s)')
    ylabel('Amplitude (μV)')
    grid on
    set(gca, 'FontSize', 8)
end

sgtitle('All Channels: Trial-aligned Responses', 'FontSize', 14)

% Save figure
saveas(gcf, fullfile(figureDir, 'all_channels_trial_aligned.jpg'))

% Create full screen figure for offset-aligned responses
figure('Position', get(0, 'ScreenSize'));

% Calculate number of channels and subplot layout
num_channels = size(data, 1);
subplot_rows = ceil(sqrt(num_channels));
subplot_cols = ceil(num_channels/subplot_rows);

% Calculate epoch window samples for offset alignment
pre_samples = round(align_offset_epoch_t(1) * srate);
post_samples = round(align_offset_epoch_t(2) * srate);
epoch_time = linspace(align_offset_epoch_t(1), align_offset_epoch_t(2), post_samples - pre_samples + 1);

% Plot each channel
for ch = 1:num_channels
    subplot(subplot_rows, subplot_cols, ch)
    hold on
    
    % Plot each trial for this channel
    for trial = 1:length(events.offset)
        % Get sample index for trial offset
        offset_idx = round(events.offset(trial) * srate);
        
        % Calculate epoch indices
        start_idx = offset_idx + pre_samples;
        end_idx = offset_idx + post_samples;
        
        % Plot if we have enough samples
        if start_idx > 0 && end_idx <= size(data, 2)
            trial_data = data(ch, start_idx:end_idx);
            plot(epoch_time, trial_data, 'k-', 'LineWidth', 0.5, 'Color', [0 0 0 0.2])
        end
    end
    
    % Add vertical line at stimulus offset
    xline(0, '--r', 'LineWidth', 1);
    yline(0, '--k');
    
    % Labels and formatting
    title(['Channel: ' channels.name{ch}])
    xlabel('Time (s)')
    ylabel('Amplitude (μV)')
    grid on
    set(gca, 'FontSize', 8)
end

sgtitle('All Channels: Offset-aligned Responses', 'FontSize', 14)

% Save figure
saveas(gcf, fullfile(figureDir, 'all_channels_offset_aligned.jpg'))

%% Plot the power spectrum of each electrode of all the data before epoch, each region in a subplot

% Create full screen figure for power spectra
figure('Position', get(0, 'ScreenSize'));

% Calculate subplot dimensions based on number of regions
subplot_dims = [3, 3]; % 3x3 grid for 9 regions

% Calculate frequency parameters
nfft = 2^nextpow2(srate); % Length of FFT
freq = linspace(0, srate/2, nfft/2+1); % Frequency vector
freq_idx = freq <= 140; % Index for frequencies up to 140 Hz

% Plot power spectrum for each region
for r = 1:length(regions)
    % Select electrodes for this region
    chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
    
    if ~isempty(chanidx)
        subplot(subplot_dims(1), subplot_dims(2), r)
        hold on
        
        % Calculate and plot power spectrum for each electrode
        for ch = 1:length(chanidx)
            % Get data for this electrode
            signal = data(chanidx(ch),:);
            
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

% Save figure
saveas(gcf, fullfile(figureDir, 'power_spectrum_by_region.jpg'))


%% check CAR data after filtering out carrier frequence (110 Hz) and line Frequency (60 Hz)

% Apply notch filter at 60 Hz and band-stop filter at 100-130 Hz
fprintf('Applying notch filter at 60 Hz and band-stop filter at 100-130 Hz...\n');

srate = channels.sampling_frequency(1); % Get sampling rate

% Design notch filter at 60 Hz
notchFilt = designfilt('bandstopiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', 59, ...
    'HalfPowerFrequency2', 61, ...
    'SampleRate', srate);

% Design band-stop filter for 100-130 Hz
bandstopFilt = designfilt('bandstopiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', 100, ...
    'HalfPowerFrequency2', 130, ...
    'SampleRate', srate);

% Apply filters to data
data_notch = filtfilt(notchFilt, double(data)')'; % First remove 60 Hz
data_filtered = filtfilt(bandstopFilt, data_notch')'; % Then remove 100-130 Hz

% Make epochs with filtered data aligned to onset
[epochs_filtered, t] = ecog_makeEpochs(data_filtered, events.onset, align_onset_epoch_t, srate);

% Plot filtered onset-aligned data for all regions in one figure
figure('Position', get(0, 'ScreenSize')); 
subplot_dims = [3, 3]; % 3x3 grid for 9 regions

for r = 1:length(regions)
    % Select electrodes for this region
    chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
    
    if ~isempty(chanidx)
        subplot(subplot_dims(1), subplot_dims(2), r)
        hold on
        plot(t, squeeze(mean(epochs_filtered(:,:,chanidx), 2)), 'LineWidth', 1);
        xline(0); yline(0)
        grid on
        title(['Region ' regions{r}])
        xlabel('Time (s)')
        ylabel('Amplitude (μV)')
        set(gca, 'XTick', min(align_onset_epoch_t):0.1:max(align_onset_epoch_t))
        set(gca, 'FontSize', 10)
    end
end
sgtitle('Trial-averaged Response by Region (60 Hz notch, 100-130 Hz band-stop)', 'FontSize', 15)

% Save onset-aligned figure
saveas(gcf, fullfile(figureDir, 'all_regions_filtered_onset.jpg'))

% Make epochs with filtered data aligned to offset
[epochs_filtered_offset, t_offset] = ecog_makeEpochs(data_filtered, offsets, align_offset_epoch_t, srate);

% Plot filtered offset-aligned data for all regions in one figure
figure('Position', get(0, 'ScreenSize')); 

for r = 1:length(regions)
    % Select electrodes for this region
    chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
    
    if ~isempty(chanidx)
        subplot(subplot_dims(1), subplot_dims(2), r)
        hold on
        plot(t_offset, squeeze(mean(epochs_filtered_offset(:,:,chanidx), 2)), 'LineWidth', 1);
        xline(0); yline(0)
        grid on
        title(['Region ' regions{r}])
        xlabel('Time from offset (s)')
        ylabel('Amplitude (μV)')
        set(gca, 'XTick', min(align_offset_epoch_t):0.1:max(align_offset_epoch_t))
        set(gca, 'FontSize', 10)
    end
end
sgtitle('Trial-averaged Response by Region Aligned to Offset (60 Hz notch, 100-130 Hz band-stop)', 'FontSize', 15)

% Save offset-aligned figure
saveas(gcf, fullfile(figureDir, 'all_regions_filtered_offset.jpg'))


% 
% %% EXTRACT BROADBAND (do it once)
% 
% outputFolder      = 'ECoGBroadband_include110Hz';
% bands             = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]; [140 150]; [150 160]; [160 170]];
% % bidsEcogBroadband(projectDir, subject, [], [], [], bands, [], [], outputFolder);
% bidsEcogBroadbandPlotAllchannels(projectDir, subject, [], [], [], bands, [], [], outputFolder);
% 
% outputFolder      = 'ECoGBroadband_exclude110Hz';
% bands             = [[70 80]; [80 90]; [90 100]; [130 140]; [140 150]; [150 160]; [160 170]];
% % bidsEcogBroadband(projectDir, subject, [], [], [], bands, [], [], outputFolder);
% bidsEcogBroadbandPlotAllchannels(projectDir, subject, [], [], [], bands, [], [], outputFolder);
% 
% %% plot broadband timecourses for temporal conditions
% 
% savePlot = 1;
% session           = 'nyuecog01';
% task              = 'temporalpattern';
% 
% % Specify one of the three preprocessed data here:
% inputFolder       = 'ECoGBroadband_exclude110Hz'; 
% description       = 'broadband';
% 
% clear specs;
% specs.epoch_t     = [-0.4 1.8]; % stimulus epoch window
% specs.base_t      = [-0.4 -0.1]; % blank epoch window
% specs.plot_ylim   = [-2 20];
% 
% specs.plot_type   = 'average';
% 
% specs.chan_names  = {'V','W','Y','Z'}; % First half of electrodes
% specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
% bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close
% specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
% bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close
% 
% specs.chan_names  = {'H','M','P','R','S'}; % Second half of electrods
% specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
% bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close
% specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
% bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close
