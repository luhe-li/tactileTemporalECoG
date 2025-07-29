
clear;
tbUse tactileTemporalECoG

bidsRootPath = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
dataPath = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR');
subject           = 'ny726';
task              = 'temporalpattern';

%% extract all session data

% data: channel x time series across 5 runs
[data, channels, events, srate] = bidsEcogGetPreprocData(dataPath, subject, [], task);

% epoch by onset: time series x trial x channel
align_onset_epoch_t = [-1, 1.5];  % stimulus epoch window
[epochs, t] = ecog_makeEpochs(data, events.onset, align_onset_epoch_t, channels.sampling_frequency(1));

save('ny726_data','data','epochs','t','channels','events','srate','subject')