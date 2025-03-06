
projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile'; 
subject = 'ny726';

%% common average reference
bidsEcogRereference(projectDir, subject);
 
%% extract broadband 
outputFolder      = 'ECoGBroadband_exclude110Hz';
bands             = [[60 70]; [70 80]; [80 90]; [120 130]; [130 140]; [160 170]; [170 180]; [180 190]];
bidsEcogBroadband(projectDir, subject, [], [], [], bands, [], [], outputFolder);

%% plot broadband timecourses for temporal conditions
savePlot = 0;

session           = 'nyuecog01';
task              = 'temporalpattern';

% Specify one of the three preprocessed data here:
inputFolder       = 'ECoGBroadband_exclude110Hz'; 
description       = 'broadband';

clear specs;

specs.epoch_t     = [-0.4 1.8]; % stimulus epoch window
specs.base_t      = [-0.4 -0.1]; % blank epoch windown
specs.plot_ylim   = [-2 10];

% all channels on the grid
specs.chan_names  = {'P','M','S'};
% specs.subplotdims = [4 8];
% specs.subplotidx  = 1:32;
specs.plot_type   = 'average';

specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close
