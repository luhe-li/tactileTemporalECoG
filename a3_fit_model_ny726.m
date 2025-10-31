
clear;
tbUse tactileTemporalECoG

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
specs.chan_names   = {'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'W4','W5','W6','W9','W10','W11','W12','P3','P5','Z9','Z10','Z13','Z14','Z15','S9','S10'};
specs.stim_names   = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6',...
    'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
[data, channel, t, srate]   = tt_prepareData(bidsDir, subject, session, task, numruns, inputFolder, description, specs);

% Generate stimulus timecourses
[stim_ts, stim_info] = tt_generateStimulusTimecourses(specs.stim_names, t);

%% 2. Model fitting to average electrodes, both models, both full fit and crossvalidation

modelfun          = {@DN,@LINEAR};

% Define parameters
options.doplots   = false;
options.xvalmode  = 0;      % 0 = none, 1 = stimulus leave-one-out
options.display   = 'iter';  % 'iter' 'final' 'off'
options.algorithm = 'bads';
options.average_elecs = true;

% Compute model fit(s); data and fits will be saved to 'derivative/modelFit/results' folder
tt_doModelFits(modelfun, stim_ts, data, channel, srate, t, stim_info, options, [], subject);

% Do crossvalidation
options.xvalmode  = 1;
tt_doModelFits(modelfun, stim_ts, data, channel, srate, t, stim_info, options, [], subject);

%% 3. Model fitting to each electrodes, DN model and full fit only to get parameter estimates confidence interval

modelfun          = @DN;
options.xvalmode  = 0;      % 0 = none, 1 = stimulus leave-one-out
options.average_elecs = false;
tt_doModelFits(modelfun, stim_ts, data, channel, srate, t, stim_info, options, [], subject);
