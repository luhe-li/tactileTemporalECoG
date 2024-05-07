

%% knobs and basic info

fitIndivElecs = true;
options           = [];
if fitIndivElecs
    options.average_elecs = false;
    datatype = 'individualelecs';
else
    options.average_elecs = true;
    datatype = 'electrodeaverages';
end

% Load or (re)compute the processed data
compute     = false;
bidsDir     = tt_bidsRootPath;
subject     = 'umcudrouwen';
session     = 'umcuiemu01';
task        = 'vtstemporalpattern';
numruns     = 1;
inputFolder = 'ECoGBroadband_exclude110Hz';
description = 'broadband';

%% 1: load ECoG data

% Select epochs and channels, average trials within stimulus condition
specs.plot_data    = false;
specs.plot_smooth  = 1; % could try some other values
specs.epoch_t      = [-0.4, 2.0];%[-0.4 1.8]; % stimulus epoch window
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

fname             = 'LINEAR';
modelfun          = str2func(fname);

% Define options
options.doplots   = true;
options.xvalmode  = 0;      % 0 = none, 1 = stimulus leave-one-out
options.display   = 'off';  % 'iter' 'final' 'off'
options.algorithm = 'bads';
% options.average_elecs = false;
options.fitaverage = false;
options.nfits     = 1000; % if fit average

% Compute model fit(s); data and fits will be saved to 'derivative/modelFit/results' folder
tt_doModelFits(modelfun, stim_ts, data, channel, srate, t, stim_info, options);

%% 3: Model evaluation

% Load data and fits
modelfun = @LINEAR;
xvalmode = 0;
% datatype = 'electrodeaverages';
datatype = 'individualelecs';
[D] = tt_loadDataForFigure(modelfun, xvalmode, datatype);

% Compute R2 and derived parameters
objFunction = modelfun;
includeDerivedParams = false;
[results] = tt_evaluateModelFit(D,includeDerivedParams);

%% 4. Plot timecourses of data and fits

% Provide a directory to save figures (optional)
saveDir = fullfile(bidsDir, 'derivatives', 'modelFit', 'figure', subject);

% Plot multiple model predictions (superimposed)
tt_plotDataAndFits(results, D.data, D.channels, D.stim, D.stim_info, D.t, D.options, saveDir, {'ONEPULSE', 'TWOPULSE'})

%% 5. Plot summed response of data and fits

% choose to plot recovery from adapatation for TWOPULSE
saveDir = fullfile(bidsDir, 'derivatives', 'modelFit', 'figure', subject);
timepointsOfInterest = [0, 2];
tt_plotSumDataAndFits(results, D.data, D.channels, D.stim, D.stim_info, D.t, D.options, saveDir, {'ONEPULSE', 'TWOPULSE'},timepointsOfInterest)

%% 6. Plot derived and fitted parameters

% Provide a directory to save figures (optional)
saveDir = fullfile(bidsDir, 'derivatives', 'modelFit', 'figure', subject);

% plot fitted parameters (see compareVSParam.m for comparison between visual and tactile datasets)
tt_plotParams(results, D.channels, D.options, saveDir);%close;
