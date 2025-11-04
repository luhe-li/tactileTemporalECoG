
clear;
% tbUse tactileTemporalECoG

compute     = false;
bidsDir     = tt_bidsRootPath;

%% 1. Load broadband data of two patients

% ny726
subject     = 'ny726';
session     = 'nyuecog01';
task        = 'temporalpattern';
numruns     = [];
inputFolder = 'ECoGBroadband_exclude110Hz';
description = 'broadband';

% Select epochs and channels, average trials within stimulus condition
specs.plot_data    = false;
specs.plot_smooth  = 1; % could try some other values
specs.epoch_t      = [-0.4, 1.8]; % stimulus epoch window
specs.base_t       = [-0.4 -0.1]; % blank epoch window
specs.chan_names   = {'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'W4','W5','W6','W9','W10','W11','W12','P3','P5','Z13','Z14','Z15','S9','S10'};
specs.stim_names   = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6',...
    'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
nyu_stim_names     = specs.stim_names;
[data1, channel1, t, srate]   = tt_prepareData(bidsDir, subject, session, task, numruns, inputFolder, description, specs);

% umcudrouwen
subject     = 'umcudrouwen';
session     = 'umcuiemu01';
task        = 'vtstemporalpattern';
numruns     = '1';
inputFolder = 'ECoGBroadband_exclude110Hz';
description = 'broadband';

% Select epochs and channels, average trials within stimulus condition
specs.plot_data    = false;
specs.plot_smooth  = 1; % could try some other values
specs.epoch_t      = [-0.4, 1.8]; % stimulus epoch window
specs.base_t       = [-0.4 -0.1]; % blank epoch window
specs.chan_names   = {'C03', 'C04', 'C11', 'C12'};
specs.stim_names   = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6',...
    'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};

% Resample data using sampling frequency from the first patient
[data2, channel2, t, srate]   = tt_prepareData(bidsDir, subject, session, task, numruns, inputFolder, description, specs, srate);

% Make sure the samples align between two patients
assert(size(data1,1) == size(data2,1))

%% 2. Concatenate data from patients

data = cat(3, data1, data2);
channel = vertcat(channel1, channel2);

%% 3. Fit model to the data averaged between patients, across electrodes

% Generate stimulus timecourses, using the stimulus info file of nyu patients
[stim_ts, stim_info] = tt_generateStimulusTimecourses(nyu_stim_names, t);

% modelfun = {@DN, @LINEAR};
modelfun = @LINEAR;

% Define parameters
options.doplots   = true;
options.xvalmode  = 0;      % 0 = none, 1 = stimulus leave-one-out
options.display   = 'off';  % 'iter' 'final' 'off'
options.algorithm = 'bads';
options.average_elecs = true;

% Compute model fit(s); data and fits will be saved to 'derivative/modelFit/results' folder
tt_doModelFits(modelfun, stim_ts, data, channel, srate, t, stim_info, options, [], 'group_average');

% Crossvalidate
options.xvalmode  = 1;      % 0 = none, 1 = stimulus leave-one-out
tt_doModelFits(modelfun, stim_ts, data, channel, srate, t, stim_info, options, [], 'group_average');

%% Fit individual electrodes

modelfun = @DN;
options.xvalmode = 0;
options.average_elecs = false;
tt_doModelFits(modelfun, stim_ts, data, channel, srate, t, stim_info, options, [], 'group_average');


% %% 4. Model evaluation
% 
% % Load data and fits
% modelfun = @DN;
% xvalmode = 1;
% datatype = 'electrodeaverages';
% [D] = tt_loadDataForFigure(modelfun, xvalmode, datatype);
% 
% % Compute R2 and derived parameters
% objFunction = modelfun;
% includeDerivedParams = false;
% [results] = tt_evaluateModelFit(D,includeDerivedParams);