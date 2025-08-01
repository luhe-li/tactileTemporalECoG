 clear; close all;

% tbUse tactileTemporalECoG
bidsRootPath = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
dataPath = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR');
subject           = 'ny726';
session           = 'nyuecog01';
task              = 'temporalpattern';

%% CHECK CAR data

[data, channels, events, srate] = bidsEcogGetPreprocData(dataPath, subject, [], task);
long = find(events.duration > 0.25);

align_onset_epoch_t     = [-0.2 0.2];  % stimulus epoch window
align_offset_epoch_t     = [-0.4 0.4];  % stimulus epoch window

% group electrode by regions
regions = {'H','M','P','R','S','V','W','Y','Z'};

% epoch by onset: time series x trial x electrode
[epochs, t] = ecog_makeEpochs(data, events.onset, align_onset_epoch_t, channels.sampling_frequency(1));

%sz = size(epochs);
%epochs = highpass(reshape(epochs, sz(1), []), 70, srate);
%epochs = reshape(epochs, sz);

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

for ii = 1:length(long)
    current_trial = long(ii);
    % Create subplots for onset-aligned data
    figure('Position', get(0, 'ScreenSize')); 
    
    % Calculate subplot layout - 3x3 grid
    subplot_dims = [4 3];
    
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
            %plot(t, squeeze(selected_epoch(:,1,chanidx)), 'LineWidth', 1);
            ts = mean(squeeze(selected_epoch(:,1,chanidx)),2);
            
            tidx = (length(t)-102):length(t);           
            fs1 = (0:length(tidx)-1)./(t(tidx(end))-t(tidx(1)));            
            A1 = abs(fft(ts(tidx)));
            P1  = angle(fft(ts(tidx)));                      
            [~, carrierIDX]=min(abs(fs1-110));
            Pcarrier1(r) = P1(carrierIDX);
            Acarrier1(r) = A1(carrierIDX);
            Pcarrier1_c(r) = P1(carrierIDX+3);
            Acarrier1_c(r) = A1(carrierIDX+3);

            tidx = 1:103;
            fs2 = (0:length(tidx)-1)./(t(tidx(end))-t(tidx(1)));
            A2 = abs(fft(ts(tidx)));
            P2  = angle(fft(ts(tidx)));  
            [~, carrierIDX]=min(abs(fs1-110));
            Pcarrier2(r) = P2(carrierIDX);
            Acarrier2(r) = A2(carrierIDX);
            Pcarrier2_c(r) = P2(carrierIDX+3);
            Acarrier2_c(r) = A2(carrierIDX+3);
            
           % plot(t, ts, 'LineWidth', 2);

            plot(fs1, A1, fs2, A2,  linewidth=2);
            
            % yyaxis right
            % plot(fs1, P1, 'ro', fs2, P2,'bo', linewidth=2);
            xlim([50 200])

            yscale('linear');
            
            xline(110);
            grid on
            
            %xlabel('Time (s)')
            %ylabel('Amplitude (μV)')
            title(['Region ' regions{r}])
            %set(gca, 'XTick', min(align_onset_epoch_t):1/110:max(align_onset_epoch_t))
            set(gca, 'FontSize', 10)
        end
    end
    
    sgtitle(sprintf('Trial %d - Onset aligned', current_trial), 'FontSize', 14)
    
    subplot(subplot_dims(1), subplot_dims(2),10)
    compassplot([Pcarrier1; Pcarrier2]', [Acarrier1; Acarrier2]',linewidth=3);
    title('110 Hz')
    legend('Stim', 'Blank')

    subplot(subplot_dims(1), subplot_dims(2),11)
    compassplot([Pcarrier1_c; Pcarrier2_c]', [Acarrier1_c; Acarrier2_c]',linewidth=3);
     title('125 Hz')
    legend('Stim', 'Blank')

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

% %% Trial-averaged time series of all electrodes, aligned to onset and offset
% 
% % Calculate subplot layout
% num_regions = length(regions);
% num_rows = ceil(sqrt(num_regions));
% num_cols = ceil(num_regions/num_rows);
% 
% % Create full screen figure for trial-averaged responses
% figure('Position', get(0, 'ScreenSize')); 
% 
% for r = 1:length(regions)
%     % select electrodes for this region
%     chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
% 
%     if ~isempty(chanidx)
%         subplot(num_rows, num_cols, r)
%         hold on
%         plot(t, squeeze(mean(epochs(:,:,chanidx), 2)),'LineWidth',1);
%         xline(0); yline(0)
%         grid on
%         title(['Region ' regions{r}])
% 
%         % Only show x label on bottom row
%         if r > (num_rows-1)*num_cols
%             xlabel('Time from onset (s)')
%         end
% 
%         % Only show y label on first column
%         if mod(r-1, num_cols) == 0
%             ylabel('Amplitude (μV)')
%         end
% 
%         set(gca, 'XTick', min(align_onset_epoch_t):0.2:max(align_onset_epoch_t))
%         set(gca, 'FontSize', 10)
%     end
% end
% sgtitle('Trial-averaged response, all electrodes grouped by region', 'FontSize', 14)
% 
% % Save figure
% figureDir = fullfile(bidsRootPath, 'derivatives', 'ECoGFigures', ['sub-' subject], 'CAR_figures');
% if ~exist(figureDir, 'dir')
%     mkdir(figureDir);
% end
% saveas(gcf, fullfile(figureDir, 'all_regions_avg_onset.png'))
% 
% % Create full screen figure for trial-averaged responses aligned to offset
% figure('Position', get(0, 'ScreenSize')); 
% 
% for r = 1:length(regions)
%     % select electrodes for this region
%     chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
% 
%     if ~isempty(chanidx)
%         subplot(num_rows, num_cols, r)
%         hold on
%         plot(t_offset, squeeze(mean(epochs_offset(:,:,chanidx), 2)),'LineWidth',1);
%         xlim(align_offset_epoch_t)
%         xline(0); yline(0)
%         grid on
%         title(['Region ' regions{r}])
% 
%         % Only show x label on bottom row
%         if r > (num_rows-1)*num_cols
%             xlabel('Time from offset (s)')
%         end
% 
%         % Only show y label on first column
%         if mod(r-1, num_cols) == 0
%             ylabel('Amplitude (μV)')
%         end
% 
%         set(gca, 'XTick', min(align_offset_epoch_t):0.2:max(align_offset_epoch_t))
%         set(gca, 'FontSize', 10)
%     end
% end
% sgtitle('Trial-averaged response, all electrodes grouped by region', 'FontSize', 14)
% 
% % Save figure
% saveas(gcf, fullfile(figureDir, 'all_regions_avg_offset.png'))

% %% filter line noise and equipment noise
% 
% srate = channels.sampling_frequency(1); % Get sampling rate
% 
% % Design notch filter at 60 Hz
% notchFilt = designfilt('bandstopiir', ...
%     'FilterOrder', 4, ...
%     'HalfPowerFrequency1', 59, ...
%     'HalfPowerFrequency2', 61, ...
%     'SampleRate', srate);
% 
% % Design band-stop filter for 100-130 Hz
% bandstopFilt = designfilt('bandstopiir', ...
%     'FilterOrder', 4, ...
%     'HalfPowerFrequency1', 100, ...
%     'HalfPowerFrequency2', 130, ...
%     'SampleRate', srate);
% 
% % Apply filters to data
% data_notch = filtfilt(notchFilt, double(data)')'; % First remove 60 Hz
% data_filtered = filtfilt(bandstopFilt, data_notch')'; % Then remove 100-130 Hz
% 
% % Make epochs with filtered data aligned to onset
% [epochs_filtered, t] = ecog_makeEpochs(data_filtered, events.onset, align_onset_epoch_t, srate);
% 
% % Plot filtered onset-aligned data for all regions in one figure
% figure('Position', get(0, 'ScreenSize')); 
% subplot_dims = [3, 3]; % 3x3 grid for 9 regions
% 
% for r = 1:length(regions)
%     % Select electrodes for this region
%     chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
% 
%     if ~isempty(chanidx)
%         subplot(subplot_dims(1), subplot_dims(2), r)
%         hold on
%         plot(t, squeeze(mean(epochs_filtered(:,:,chanidx), 2)), 'LineWidth', 1);
%         xline(0); yline(0)
%         grid on
%         title(['Region ' regions{r}])
%         xlabel('Time from onset (s)')
%         ylabel('Amplitude (μV)')
%         set(gca, 'XTick', min(align_onset_epoch_t):0.1:max(align_onset_epoch_t))
%         set(gca, 'FontSize', 10)
%     end
% end
% sgtitle('Trial-averaged response (60 Hz notch, 100-130 Hz band-stop)', 'FontSize', 15)
% 
% % Save onset-aligned figure
% saveas(gcf, fullfile(figureDir, 'all_regions_filtered_onset.jpg'))
% 
% % Make epochs with filtered data aligned to offset
% [epochs_filtered_offset, t_offset] = ecog_makeEpochs(data_filtered, offsets, align_offset_epoch_t, srate);
% 
% % Plot filtered offset-aligned data for all regions in one figure
% figure('Position', get(0, 'ScreenSize')); 
% 
% for r = 1:length(regions)
%     % Select electrodes for this region
%     chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
% 
%     if ~isempty(chanidx)
%         subplot(subplot_dims(1), subplot_dims(2), r)
%         hold on
%         plot(t_offset, squeeze(mean(epochs_filtered_offset(:,:,chanidx), 2)), 'LineWidth', 1);
%         xline(0); yline(0)
%         grid on
%         title(['Region ' regions{r}])
%         xlabel('Time from offset (s)')
%         ylabel('Amplitude (μV)')
%         xlim(align_offset_epoch_t)
%         set(gca, 'XTick', min(align_offset_epoch_t):0.1:max(align_offset_epoch_t))
%         set(gca, 'FontSize', 10)
%     end
% end
% sgtitle('Trial-averaged response (60 Hz notch, 100-130 Hz band-stop)', 'FontSize', 15)
% 
% % Save offset-aligned figure
% saveas(gcf, fullfile(figureDir, 'all_regions_filtered_offset.jpg'))


%%  Calculate ITPC to check for phase locking

% Calculate ITPC on non-filtered data
% Define frequency range of interest
freqs = 1:170;
n_freqs = length(freqs);
n_times = size(epochs,1);
n_trials = size(epochs,2);
n_chans = size(epochs,3);
% Initialize matrices
itpc = zeros(n_freqs, n_times, n_chans);

% Calculate ITPC for each channel using Morlet wavelets
for chan = 1:n_chans
    % Get data for this channel
    chan_data = squeeze(epochs(:,:,chan));
    
    % Initialize time-frequency matrix
    tf_complex = zeros(n_freqs, n_times, n_trials);
    
    % Calculate wavelet transform for each trial
    for trial = 1:n_trials
        % Get trial data
        trial_data = chan_data(:,trial);
        
        % Calculate wavelet transform for each frequency
        for f = 1:n_freqs
            % Create Morlet wavelet
            freq = freqs(f);
            sigma = 7/(2*pi*freq); % Width of Gaussian
            wavelet = exp(2*1i*pi*freq.*(-(3*sigma):1/srate:(3*sigma))).*exp(-(-3*sigma:1/srate:3*sigma).^2./(2*sigma^2));
            wavelet = wavelet./sum(abs(wavelet));
            
            % Convolve with trial data
            conv_result = conv(trial_data, wavelet, 'same');
            tf_complex(f,:,trial) = conv_result;
            
        end
    end
    
    % Calculate ITPC by normalizing phase angles and averaging across trials
    tf_norm = tf_complex ./ abs(tf_complex); % Normalize to unit circle
    itpc(:,:,chan) = abs(mean(tf_norm,3));
    
end

%
% % Plot ITPC for each channel
% figure('Position', get(0, 'ScreenSize'));
% 
% % Calculate subplot layout
% n_total_chans = size(itpc,3);
% n_rows = ceil(sqrt(n_total_chans));
% n_cols = ceil(n_total_chans/n_rows);
% 
% for chan = 1:n_total_chans
%     subplot(n_rows, n_cols, chan)
% 
%     % Plot time-frequency map for this channel
%     imagesc(t, freqs(1:50), itpc(:,1:50,chan));
%     axis xy
%     colorbar
%     caxis([0 0.5]) % Set color limits for ITPC (0-1)
% 
%     title(['Channel ' channels.name{chan}])
%     xlabel('Time from onset (s)')
%     ylabel('Frequency (Hz)')
% 
%     % Add vertical line at stimulus onset
%     hold on
%     line([0 0], ylim, 'Color', 'w', 'LineStyle', '--')
% end
% 
% 
% % Save ITPC figure
% saveas(gcf, fullfile(figureDir, 'all_channels_ITPC.jpg'))
