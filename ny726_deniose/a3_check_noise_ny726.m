clear; close all;

tbUse tactileTemporalECoG
bidsRootPath = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
dataPath = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR');
subject           = 'ny726';
session           = 'nyuecog01';
task              = 'temporalpattern';

%% CHECK CAR data

[data, channels, events, srate] = bidsEcogGetPreprocData(dataPath, subject, [], task);

align_onset_epoch_t     = [-0.2 0.2];  % stimulus epoch window
align_offset_epoch_t     = [-0.4 0.4];  % stimulus epoch window

% group electrode by regions
regions = {'H','M','P','R','S','V','W','Y','Z'};

% epoch by onset: time series x trial x electrode
[epochs, t] = ecog_makeEpochs(data, events.onset, align_onset_epoch_t, channels.sampling_frequency(1));

% epoch by offset: time series x trial x electrode
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
[epochs_offset, t_offset] = ecog_makeEpochs(data, offsets, align_offset_epoch_t, channels.sampling_frequency(1));

%% check one trial aligned by onset 

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
            set(gca, 'XTick', min(align_onset_epoch_t):0.01:max(align_onset_epoch_t))
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
    %         set(gca, 'XTick', min(align_offset_epoch_t):0.01:max(align_offset_epoch_t))
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
