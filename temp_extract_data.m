
clear;
tbUse tactileTemporalECoG

bidsRootPath = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
dataPath = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR');
subject           = 'ny726';
task              = 'temporalpattern';

%% extract all session data

[data, channels, events, srate] = bidsEcogGetPreprocData(dataPath, subject, [], task);