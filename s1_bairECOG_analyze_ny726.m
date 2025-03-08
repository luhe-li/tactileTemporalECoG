
projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile'; 
subject = 'ny726';

%% common average reference (do it once)
% bidsEcogRereference(projectDir, subject);
 
%% extract broadband (do it once)
% outputFolder      = 'ECoGBroadband_exclude110Hz';
% bands             = [[60 70]; [70 80]; [80 90]; [120 130]; [130 140]; [160 170]; [170 180]; [180 190]];
% bidsEcogBroadband(projectDir, subject, [], [], [], bands, [], [], outputFolder);

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
specs.plot_ylim   = [-2 20];

% all channels on the grid
% specs.subplotdims = [4 8];
% specs.subplotidx  = 1:32;
specs.plot_type   = 'average';

specs.chan_names  = {'V','W','Y','Z'}; % First half of electrodes
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close

specs.chan_names  = {'H','M','P','R','S'}; % Second half of electrods
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close

% subeset of channels
specs.chan_names  = {'P','M','S'}; % Region of interests based on figure
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close

%% plot spectra of the re-referenced data

clear specs;
specs.epoch_t     = [-1 2]; % stimulus epoch window
% specs.fft_blank_t = [-0.8 -0.1]; % fft blank epoch window
% specs.fft_w       = 0.4;
% specs.fft_ov      = 0.1;
% specs.plot_ylim   = [10^-3 10^3];

specs.chan_names  = {'V','W','Y','Z'}; % First half of electrodes 
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
specs.fft_stim_t  = [0 1.2];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], [], [], specs, savePlot);
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
specs.fft_stim_t  = [0 1.8];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], [], [], specs, savePlot);

specs.chan_names  = {'H','M','P','R','S'}; % Second half of electrods
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
specs.fft_stim_t  = [0 1.2];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], [], [], specs, savePlot);
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
specs.fft_stim_t  = [0 1.8];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], [], [], specs, savePlot);

%% Plot the time series for selected electrodes

clear specs;
specs.epoch_t     = [-0.4 1.8]; % stimulus epoch window
specs.base_t      = [-0.4 -0.1]; % blank epoch window
specs.plot_ylim   = [-2 30];
specs.plot_type   = 'averageSE';

specs.chan_names  = {'M11', 'M12', 'P8', 'P9'};
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close

specs.plot_ylim   = [-2 20];
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close

%% Average time series for selected electrodes, each panel per condition

bidsDir     = tt_bidsRootPath;

clear specs;
specs.plot_data    = true;
specs.plot_smooth  = 1; % could try some other values
specs.epoch_t      = [-0.4, 1.8]; % stimulus epoch window
specs.base_t       = [-0.4 -0.1]; % blank epoch window
specs.chan_names  = {'M'};
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
[data, channels, t, srate]   = tt_prepareData(bidsDir, subject, session, task, [], inputFolder, description, specs);

specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
[data, channels, t, srate]   = tt_prepareData(bidsDir, subject, session, task, [], inputFolder, description, specs);

%%
% Generate stimulus timecourses
[stim_ts, stim_info] = tt_generateStimulusTimecourses(specs.stim_names, t);

specs.average_elecs = true;
%%
tt_plotData(data, channels, t, specs, savePlot)