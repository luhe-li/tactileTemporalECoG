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

% check the data
figure;plot(data(triggerChannel,:)/max(data(triggerChannel,:))); 
title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);

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

% % Plot first 3 and last 3 trials for each block in a 5x6 grid
% figure('Position', [100 100 1500 1000]);
% 
% % Calculate samples for 1.5s epoch
% epoch_length = 0.1; % seconds
% samples_per_epoch = round(epoch_length * hdr.Fs);
% epoch_time = (-10:samples_per_epoch-1)/hdr.Fs;
% 
% for block = 1:numBlocks
%     % Get first 3 and last 3 trigger indices for this block
%     block_triggers = blocks_idx{block};
%     first_three = block_triggers(1:3);
%     last_three = block_triggers(end-2:end);
%     all_triggers = [first_three; last_three];
% 
%     % Plot each trial in its own subplot
%     for trial = 1:6
%         subplot(5, 6, (block-1)*6 + trial)
%         trigger_sample = all_triggers(trial);
% 
%         if trigger_sample + samples_per_epoch <= size(data, 2)
%             epoch_data = data(triggerChannel, (trigger_sample-10):trigger_sample+samples_per_epoch-1);
%             plot(epoch_time, epoch_data, 'k-', 'LineWidth', 1)
% 
%             if trial <= 3
%                 title(sprintf('Block %d\nFirst Trial %d', block, trial))
%             else
%                 title(sprintf('Block %d\nLast Trial %d', block, trial-3))
%             end
% 
%             xlabel('Time (s)')
%             ylabel('Amplitude')
%             grid on
%         end
%     end
% end
% 
% sgtitle('First Channel Activity Around Events (First 3 and Last 3 Trials per Block)')

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

%% check one random trial aligned by onset and offset

% Randomly select one trial
selected_trial = randi([1 720], 1);

% Extract epochs for selected trial
selected_epoch = epochs(:,selected_trial,:);
selected_epoch_offset = epochs_offset(:,selected_trial,:);

% Create subplots for onset-aligned data
figure('Position', get(0, 'ScreenSize')); 

% Calculate subplot layout - 3x3 grid
subplot_dims = [3 3];

% Plot each region
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

sgtitle(sprintf('Random trial %d - Onset aligned', selected_trial), 'FontSize', 14)

% Save figure with high resolution
saveas(gcf, fullfile(figureDir, sprintf('random_trial%d_onset.jpg', selected_trial)))

% Create subplots for offset-aligned data
figure('Position', get(0, 'ScreenSize')); 

% Plot each region
for r = 1:length(regions)
    % select electrodes for this region
    chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
    
    if ~isempty(chanidx)
        subplot(subplot_dims(1), subplot_dims(2), r)
        hold on
        
        % Plot each electrode's response for this trial
        plot(t_offset, squeeze(selected_epoch_offset(:,1,chanidx)), 'LineWidth', 1);
        xline(0); yline(0)
        grid on
        
        xlabel('Time from offset (s)')
        ylabel('Amplitude (μV)')
        title(['Region ' regions{r}])
        set(gca, 'XTick', min(align_offset_epoch_t):0.2:max(align_offset_epoch_t))
        set(gca, 'FontSize', 10)
    end
end

sgtitle(sprintf('Random trial %d - Offset aligned', selected_trial), 'FontSize', 14)

% Save figure
saveas(gcf, fullfile(figureDir, sprintf('random_trial%d_offset.jpg', selected_trial)))


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
