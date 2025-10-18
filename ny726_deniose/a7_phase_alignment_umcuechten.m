clear; close all;

tbUse tactileTemporalECoG
bidsRootPath = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
dataPath = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR');

subject = 'umcuechten';
session           = 'umcuiemu01';
task              = 'vtstemporalpattern';
inputFolder       = 'ECoGCAR';
description       = 'reref';

figureDir = fullfile(pwd, mfilename);
if ~exist(figureDir,'dir'); mkdir(figureDir); end

[data, channels, events, srate] = bidsEcogGetPreprocData(dataPath, subject, [], task);
long = find(events.duration > 0.25);
fs = channels.sampling_frequency(1);
%%
% epoch by onset: time series x trial x electrode
align_onset_epoch_t     = [-0.2 0.2];  % stimulus epoch window
[epochs, t] = ecog_makeEpochs(data, events.onset, align_onset_epoch_t, fs);

% group electrode by regions
regions = {'C','D'};

%% plot phase of electrodes within a trial

for ii = 1:length(long)

    current_trial = long(ii);
    % Create subplots for onset-aligned data
    figure('Position', get(0, 'ScreenSize')); 
    
    % Calculate subplot layout - 3x3 grid
    subplot_dims = [2,2];
    
    % Extract epochs for current trial
    selected_epoch = epochs(:,current_trial,:);
    
    % Plot each region for onset-aligned data
    for r = 1:length(regions)
        % select electrodes for this region
        chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
        
        if ~isempty(chanidx)
            subplot(subplot_dims(1), subplot_dims(2), r)
            hold on
            
            % Initialize arrays to store phase and amplitude for each channel
            Pcarrier1 = zeros(length(chanidx), 1);
            Acarrier1 = zeros(length(chanidx), 1);
            Pcarrier1_c = zeros(length(chanidx), 1);
            Acarrier1_c = zeros(length(chanidx), 1);
            Pcarrier2 = zeros(length(chanidx), 1);
            Acarrier2 = zeros(length(chanidx), 1);
            Pcarrier2_c = zeros(length(chanidx), 1);
            Acarrier2_c = zeros(length(chanidx), 1);
            
            % Process each channel separately
            for ch = 1:length(chanidx)
                ts = squeeze(selected_epoch(:,1,chanidx(ch)));
                
                tidx = (length(t)-102):length(t);           
                fs1 = (0:length(tidx)-1)./(t(tidx(end))-t(tidx(1)));
                A1 = abs(fft(ts(tidx)));
                P1 = angle(fft(ts(tidx)));                      
                [~, carrierIDX]=min(abs(fs1-110));
                Pcarrier1(r, ch) = P1(carrierIDX);
                Acarrier1(r, ch) = A1(carrierIDX);
                Pcarrier1_c(r, ch) = P1(carrierIDX+3);
                Acarrier1_c(r, ch) = A1(carrierIDX+3);

                tidx = 1:103;
                fs2 = (0:length(tidx)-1)./(t(tidx(end))-t(tidx(1)));
                A2 = abs(fft(ts(tidx)));
                P2 = angle(fft(ts(tidx)));  
                [~, carrierIDX]=min(abs(fs1-110));
                Pcarrier2(r, ch) = P2(carrierIDX);
                Acarrier2(r, ch) = A2(carrierIDX);
                Pcarrier2_c(r, ch) = P2(carrierIDX+3);
                Acarrier2_c(r, ch) = A2(carrierIDX+3);
                
                % Plot frequency spectrum for each channel
                plot(fs1, A1, fs2, A2, 'LineWidth', 1);
            end
            
            xlabel('Frequency')
            ylabel('Amplitude')
            xlim([1 150])
            yscale('linear');
            xline(110);
            grid on
            title(['Region ' regions{r}])
            set(gca, 'FontSize', 10)
        end
    end
    
    sgtitle(sprintf('Trial %d - Onset aligned', current_trial), 'FontSize', 14)
    
    % Plot compass plots with all channels from all regions
    subplot(subplot_dims(1), subplot_dims(2),3)
    compassplot([Pcarrier1(:)'; Pcarrier2(:)']', [Acarrier1(:)'; Acarrier2(:)']',linewidth=3);
    title('110 Hz')
    legend('Stim', 'Blank')

    subplot(subplot_dims(1), subplot_dims(2),4)
    compassplot([Pcarrier1_c(:)'; Pcarrier2_c(:)']', [Acarrier1_c(:)'; Acarrier2_c(:)']',linewidth=3);
     title('125 Hz')
    legend('Stim', 'Blank')

    % Wait for mouse click
    saveas(gcf, fullfile(figureDir, sprintf('trila%i.jpg',ii)))
    w = waitforbuttonpress;
    if w == 0 % Mouse click
        close all;
        continue;
    else % Keyboard press
        close all;
        break;
    end

end
