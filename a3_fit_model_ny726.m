
clear;
tbUse tactileTemporalECoG

addpath(genpath(pwd, 'codes'))
compute     = false;
bidsDir     = tt_bidsRootPath;
subject     = 'ny726';
session     = 'nyuecog01';
task        = 'temporalpattern';
numruns     = [];
inputFolder = 'ECoGBroadband_exclude110Hz';
description = 'broadband';

%% 1. Load broadband data

% Select epochs and channels, average trials within stimulus condition
specs.plot_data    = false;
specs.plot_smooth  = 1; % could try some other values
specs.epoch_t      = [-0.4, 1.8]; % stimulus epoch window
specs.base_t       = [-0.4 -0.1]; % blank epoch window
specs.chan_names   = {'C03', 'C04', 'C11', 'C12'};
specs.stim_names   = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6',...
    'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
[data, channel, t, srate]   = tt_prepareData(bidsDir, subject, session, task, numruns, inputFolder, description, specs);

% Generate stimulus timecourses
[stim_ts, stim_info] = tt_generateStimulusTimecourses(specs.stim_names, t);

%% 2. Fit model to the s1 and m1 data


modelfun = '@DN';

% Define parameters
options.doplots   = true;
options.xvalmode  = 0;      % 0 = none, 1 = stimulus leave-one-out
options.display   = 'iter';  % 'iter' 'final' 'off'
options.algorithm = 'bads';
options.average_elecs = true;
options.nfits     = 100; % if fit average

% Compute model fit(s); data and fits will be saved to 'derivative/modelFit/results' folder
tt_doModelFits(modelfun, stim_ts, s1_data, channel, srate, t, stim_info, options, [], 's1');
tt_doModelFits(modelfun, stim_ts, m1_data, channel, srate, t, stim_info, options, [], 'm1');

%% 3. Evaluate model fits

xvalmode = 1;
datatype = 'electrodeaverages';

[D{1}] = tt_loadDataForFigure(modelfun, xvalmode, datatype, [], 's1');
[D{2}] = tt_loadDataForFigure(modelfun, xvalmode, datatype, [], 'm1');

saveDir = fullfile(bidsDir, 'derivatives', 'modelFit', 'figure', subject);
tt_plotTactileParams(D, modelfun, saveDir);