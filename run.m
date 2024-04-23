

%% 1: load ECoG data

% Load or (re)compute the processed data
compute     = false;
bidsDir     = tt_bidsRootPath;
subject     = 'umcudrouwen';
session     = 'umcuiemu01';
task        = 'vtstemporalpattern';
numruns     = 1;
inputFolder = 'ECoGBroadband_exclude110Hz';
description = 'broadband';

% Select epochs and channels, average trials within stimulus condition
specs.plot_data    = false;
specs.plot_smooth  = 1;
specs.epoch_t      = [-0.4 1.8]; % stimulus epoch window
specs.base_t       = [-0.4 -0.1]; % blank epoch window
specs.chan_names   = {'C03', 'C04', 'C11', 'C12'};
specs.stim_names   = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6',...
                     'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
[data, channel, t, srate]   = tt_prepareData(bidsDir, subject, session, task, numruns, inputFolder, description, specs);

% Generate stimulus timecourses
[stim_ts, stim_info] = tt_generateStimulusTimecourses(specs.stim_names, t);

%% 2: Model fitting
% try fitting DN model by reusing Iris' codes, removing probabilistic
% resample step

fname             = 'DN';
modelfun          = str2func(fname);

% Define options
options           = [];
options.doplots   = false;
options.xvalmode  = 0;      % 0 = none, 1 = stimulus leave-one-out
options.display   = 'off';  % 'iter' 'final' 'off'
options.algorithm = 'bads';
options.average_elecs = true;
options.fitaverage = false;
options.nfits     = 1000; % if fit average 

% Compute model fit(s); data and fits will be saved to 'derivative/analysis/results' folder
tt_doModelFits(modelfun, stim_ts, data, channel, srate, t, stim_info, options);

%% 3: Model evaluation

% Load data and fits
modelfun = @DN;
xvalmode = 0;
datatype = 'electrodeaverages';
[D] = tt_loadDataForFigure(modelfun, xvalmode, datatype);

% Compute R2 and derived parameters
objFunction = modelfun;
includeDerivedParams = false;
[results] = tt_evaluateModelFit(D,includeDerivedParams);

%% 4. Plot timecourses and fits

% Provide a directory to save figures (optional)
saveDir = fullfile(bidsDir, 'derivatives', 'modelFit', 'figure', subject);

% Plot multiple model predictions (superimposed)
tt_plotDataAndFits(results, D.data, D.channels, D.stim, D.stim_info, D.t, D.options, saveDir, {'ONEPULSE', 'TWOPULSE'})

%% 5. Plot derived and fitted parameters

