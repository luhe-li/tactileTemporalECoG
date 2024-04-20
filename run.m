

%% 1: load ECoG data

compute     = false;
bidsDir     = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
subject     = 'umcudrouwen';
session     = 'umcuiemu01';
task        = 'vtstemporalpattern';
numruns     = 1;
inputFolder = 'ECoGBroadband_exclude110Hz';
description = 'broadband';

specs.plot_data    = true;
specs.plot_smooth  = 1;
specs.epoch_t     = [-0.4 1.8]; % stimulus epoch window
specs.base_t      = [-0.4 -0.1]; % blank epoch window
specs.chan_names  = {'C03', 'C04', 'C11', 'C12'};
specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6',...
                     'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};

[data, channel, epochs]   = tt_prepareData(bidsDir, subject, session, task, numruns, inputFolder, description, specs);

%% 2: Model fitting


%% 3: Model evaluation


%% 4. Plot timecourses and fits


%% 5. Plot derived and fitted parameters