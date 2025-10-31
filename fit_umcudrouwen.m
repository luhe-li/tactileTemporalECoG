
%% knobs and basic info
options       = [];

% Load or (re)compute the processed data
compute     = false;
bidsDir     = tt_bidsRootPath;
subject     = 'umcudrouwen';
session     = 'umcuiemu01';
task        = 'vtstemporalpattern';
numruns     = '1';
inputFolder = 'ECoGBroadband_exclude110Hz';
description = 'broadband';

%% 1: load ECoG data

% Select epochs and channels, average trials within stimulus condition
specs.plot_data    = false;
specs.plot_smooth  = 1; % could try some other values
specs.epoch_t      = [-0.4, 2.0]; % stimulus epoch window
specs.base_t       = [-0.4 -0.1]; % blank epoch window
specs.chan_names   = {'C03', 'C04', 'C11', 'C12'};
specs.stim_names   = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6',...
    'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
[data, channel, t, srate]   = tt_prepareData(bidsDir, subject, session, task, numruns, inputFolder, description, specs);

% Fix inconsistency between experiments, use the 'ONE-PULSE-1' as sample
% condition names
if ~strcmp(specs.stim_names, 'ONE-PULSE-1')
   specs.stim_names = regexprep(specs.stim_names, 'PULSE', '-PULSE');
end

% Generate stimulus timecourses
[stim_ts, stim_info] = tt_generateStimulusTimecourses(specs.stim_names, t);

%% 2: Model fitting to average electrodes, both models, both full fit and crossvalidation

% Fitting DN/Linear model by reusing Iris' codes, removing probabilistic
% resample step
modelfun          = {@DN,@LINEAR};

% Define options
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

%% 3. Model fitting to each electrodes, DN model only to get parameter estimates confidence interval

modelfun          = @DN;
options.xvalmode  = 0;      % 0 = none, 1 = stimulus leave-one-out
options.average_elecs = false;
tt_doModelFits(modelfun, stim_ts, data, channel, srate, t, stim_info, options, [], subject);

% %% 4: Model evaluation
% 
% % Load data and fits
% modelfun = @DN;
% xvalmode = 0;
% [D] = tt_loadDataForFigure(modelfun, xvalmode, datatype);
% 
% % Compute R2 and derived parameters
% objFunction = modelfun;
% includeDerivedParams = false;
% [results] = tt_evaluateModelFit(D,includeDerivedParams);
% 
% %% 4. Plot timecourses of data and fits
% 
% % Provide a directory to save figures (optional)
% saveDir = fullfile(bidsDir, 'derivatives', 'modelFit', 'figure', subject);
% 
% % Plot multiple model predictions (superimposed)
% tt_plotDataAndFits(results, D.data, D.channels, D.stim, D.stim_info, D.t, D.options, saveDir, {'ONEPULSE', 'TWOPULSE'})

% %% 5. Plot summed response of data and fits
% 
% % choose to plot recovery from adapatation for TWOPULSE
% saveDir = fullfile(bidsDir, 'derivatives', 'modelFit', 'figure', subject);
% timepointsOfInterest = [0, 2];
% tt_plotSumDataAndFits(results, D.data, D.channels, D.stim, D.stim_info, D.t, D.options, saveDir, {'ONEPULSE', 'TWOPULSE'},timepointsOfInterest)
% 
% %% 6. Plot derived and fitted parameters
% 
% % Provide a directory to save figures (optional)
% saveDir = fullfile(bidsDir, 'derivatives', 'modelFit', 'figure', subject);
% 
% % plot fitted parameters (see s2_compareVTparam.m for comparison between visual and tactile datasets)
% tt_plotParams(results, D.channels, D.options, saveDir); %close;
