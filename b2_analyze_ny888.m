% sub-ny888 is recoded as sub-p15 when converting to BIDS

clear;
tbUse tactileTemporalECoG

projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
subject = 'ny888';

sessions = 'nyuecog01';
outputFolder = 'demean_CAR';
tasks = 'temporalpattern';
runnums = [];

addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/fileio/'))
addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/utilities/'))

%% apply demean CAR (do it one)

bidsEcogRereference(projectDir, subject, sessions, tasks, runnums, outputFolder)

%% extract/check broadband (do it once)

inputFolder       = 'demean_CAR';
outputFolder      = 'ECoGBroadband_exclude110Hz';
bands             = [[70 80]; [80 90]; [90 100]; [130 140]; [140 150]; [150 160]; [160 170]];
bidsEcogBroadband(projectDir, subject, sessions, tasks, [], bands, [], inputFolder, outputFolder);

%% plot broadband timecourses for temporal conditions

savePlot = 1;
session           = 'nyuecog01';
task              = 'temporalpattern';

% Specify one of the three preprocessed data here:
inputFolder       = 'ECoGBroadband_exclude110Hz'; 
description       = 'broadband';

%% plot broadband timecourses for temporal conditions

clear specs;
specs.epoch_t     = [-0.4 1.8]; % stimulus epoch window
specs.base_t      = [-0.4 -0.1]; % blank epoch window

specs.plot_type   = 'average';

specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); 
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); 
