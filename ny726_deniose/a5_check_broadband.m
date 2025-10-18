
clear;
tbUse tactileTemporalECoG

projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
subject = 'ny726';
inputFolder = 'ECoGCAR';
outputFolder = 'ECoGBroadband_exclude110Hz';
bandopts = [[70 80]; [80 90]; [90 100]; [130 140]; [140 150]; [150 160]; [160 170]];

%% specifiy data

[sessions] = bidsSpecifySessions(projectDir, subject);
idx = find(contains(lower(sessions), {'ecog', 'iemu'}));
sessions = sessions(idx);
tasks = [];
runnums = 5;
description = 'reref';
[session, tasks, runnums] = bidsSpecifyData(projectDir, subject, sessions{1}, tasks, runnums);
dataPath = fullfile(projectDir, 'derivatives', inputFolder);

% choose the first run
task = tasks{1};
runnum = '01';

%% read data
[data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, subject, session, task, runnum, description);
chan_index = find(contains(lower(channels.type), {'ecog', 'seeg'}));
data_bb = data;
temp_data = data(chan_index,:)'; % inputs to ecog_extractBroadband should be time x channels

%% apply notch filter at 110 hz
srate = hdr.Fs;
x = ecog_notch(temp_data, srate, 110, 2);

% format the bands input
if isa(bandopts, 'cell')

    % Entire range for broadband
    band_rg  = bandopts{1};

    % Bin width
    band_w   = bandopts{2};

    % All bins
    lb       = band_rg(1):band_w:band_rg(2)-band_w;
    ub       = lb+band_w;
    bands   = [lb; ub]';
else
    bands = bandopts;
end

% band pass filter each sub-band
bp  = repmat({zeros(size(x), 'double')},size(bands,1),1);

% replace any nans in the data with zeros
if any(isnan(x(:)))
    warning('[%s] Found nans in the data! Replacing with zeros. Please check data for extensive nans \n', mfilename);
    x(isnan(x)) = 0;
end

% bandpass filter
for ii = 1:size(bands,1)
    fprintf('[%s] Filtering signal in band %d-%d\n', mfilename, bands(ii,1),bands(ii,2));
    bp{ii} = butterpass_eeglabdata_nyu(x,bands(ii,:),srate);
end
bp = cat(ndims(x)+1,bp{:});

%% check frequency removal

% average across bands
mbp = mean(bp,3); 
for i = 1:size(mbp,2)

    electrode_data = mbp(:,i);
    [pxx, f] = pwelch(electrode_data, [], [], [], srate);
    idx_range = f >= 50 & f <= 180;

    freqs_in_range(i,:) = f(idx_range);
    power_in_range(i,:) = pxx(idx_range);

end
ch_slc = 1:60;
figure; plot(freqs_in_range(1,:),power_in_range(ch_slc,:));

ch_slc = 60:120;
figure; plot(freqs_in_range(1,:),power_in_range(ch_slc,:));

%% plot time series of different bands, of electrodes within one region, at one trial

regions = {'H','M','P','R','S','V','W','Y','Z'};
r = 4;
chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
epoch_time = [-0.2, 1.8]; % second

clear epochs
for bb = 1:size(bandopts,1)
    % epochs: (samples x trials x channels x bands)   
    [epochs(:,:,:,bb), t] = ecog_makeEpochs(squeeze(bp(:,chanidx,bb))', events.onsets, epoch_time, srate);
end

cmap = parula(size(bandopts,1));
for tt = 1:size(epochs,2)
    figure;set(gcf, 'Position', get(0, 'Screensize'));
    hold on
    h = []; 
    for bb = 1:size(bandopts,1)
        plot(t, squeeze(epochs(:,tt,:,bb)), 'Color', cmap(bb,:));
        h(bb) = plot(t, squeeze(epochs(:,tt,1,bb)), 'Color', cmap(bb,:));
    end
    xlabel('Time to trial onset (s)')
    xline(0)
    xlim(epoch_time)
    ylabel('uV')
    legend(h(:), arrayfun(@(x,y) sprintf('%d-%dHz', x, y), bandopts(:,1), bandopts(:,2), 'UniformOutput', false))
    trial_name = events.trial_name{tt};
    title(trial_name)
    if contains(trial_name, 'ONE-PULSE')
        stim_dur = events.duration(tt);
        plot([0, stim_dur],[10, 10],'r','LineWidth', 2,'DisplayName','Stimulus')
    elseif contains(trial_name, 'TWO-PULSE')
        plot([0, 0.2],[10, 10],'r','LineWidth', 2,'DisplayName','Stimulus')
        ISI = events.ISI(tt);
        second_start = 0.2 + ISI;
        plot([second_start, second_start+0.2],[10, 10],'r','LineWidth', 2,'HandleVisibility','off')
    end

    waitforbuttonpress

end

%% plot hilbert transformed envelop of different bands, of electrodes within one region, at one trial

pw = abs(hilbert(bp)).^2;

regions = {'H','M','P','R','S','V','W','Y','Z'};
r = 4;
chanidx = find(cellfun(@(x) ~isempty(x), strfind(channels.name, regions{r})));
epoch_time = [-0.2, 1.8]; % second

clear epochs
for bb = 1:size(bandopts,1)
    % epochs: (samples x trials x channels x bands)   
    [epochs(:,:,:,bb), t] = ecog_makeEpochs(squeeze(pw(:,chanidx,bb))', events.onsets, epoch_time, srate);
end

cmap = parula(size(bandopts,1));
for tt = 1:size(epochs,2)
    figure;set(gcf, 'Position', get(0, 'Screensize'));
    hold on
    h = []; 
    for bb = 1:size(bandopts,1)
        plot(t, squeeze(epochs(:,tt,:,bb)), 'Color', cmap(bb,:),'LineWidth',1.5);
        h(bb) = plot(t, squeeze(epochs(:,tt,1,bb)), 'Color', cmap(bb,:));
    end
    xlabel('Time to trial onset (s)')
    xline(0)
    xlim(epoch_time)
    ylabel('uV')
    legend(h(:), arrayfun(@(x,y) sprintf('%d-%dHz', x, y), bandopts(:,1), bandopts(:,2), 'UniformOutput', false))
    trial_name = events.trial_name{tt};
    title(trial_name)
    if contains(trial_name, 'ONE-PULSE')
        stim_dur = events.duration(tt);
        plot([0, stim_dur],[10, 10],'r','LineWidth', 2,'DisplayName','Stimulus')
    elseif contains(trial_name, 'TWO-PULSE')
        plot([0, 0.2],[10, 10],'r','LineWidth', 2,'DisplayName','Stimulus')
        ISI = events.ISI(tt);
        second_start = 0.2 + ISI;
        plot([second_start, second_start+0.2],[10, 10],'r','LineWidth', 2,'HandleVisibility','off')
    end

    waitforbuttonpress

end
