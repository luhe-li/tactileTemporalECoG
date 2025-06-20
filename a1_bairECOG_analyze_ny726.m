
clear;
tbUse tactileTemporalECoG

projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile'; 
subject = 'ny726';

% addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/fileio/'))
% addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/utilities/'))

%% common average reference (do it once)
% bidsEcogRereference(projectDir, subject);

%% extract/check broadband (do it once)

outputFolder      = 'ECoGBroadband_exclude110Hz';
bands             = [[70 80]; [80 90]; [90 100]; [130 140]; [140 150]; [150 160]; [160 170]];
bidsEcogBroadband(projectDir, subject, [], [], [], bands, [], [], outputFolder);

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

% all channels on the gri
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

%% Plot the overlapped time series for selected electrodes

specs.plot_type   = 'averageSE';

specs.chan_names  = {'M11', 'M12', 'P8', 'P9'};
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close

specs.plot_ylim   = [-2 20];
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close

%% plot broadband spectra

clear specs;
specs.epoch_t     = [-1 2]; % stimulus epoch window
specs.fft_blank_t = [-0.8 -0.1]; % fft blank epoch window
% specs.fft_w       = 0.4;
% specs.fft_ov      = 0.1;
% specs.plot_ylim   = [10^-3 10^3];

specs.chan_names  = {'V','W','Y','Z'}; % First half of electrodes 
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
specs.fft_stim_t  = [0 1.2];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot);
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
specs.fft_stim_t  = [0 1.8];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot);

specs.chan_names  = {'H','M','P','R','S'}; % Second half of electrods
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
specs.fft_stim_t  = [0 1.2];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot);
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
specs.fft_stim_t  = [0 1.8];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot);

%% plot broadband spectra of CAR

inputFolder       = 'ECoGCAR';
description       = 'reref';

clear specs;
specs.epoch_t     = [-1 2]; % stimulus epoch window
specs.fft_blank_t = [-0.8 -0.1]; % fft blank epoch window
% specs.fft_w       = 0.4;
% specs.fft_ov      = 0.1;
% specs.plot_ylim   = [10^-3 10^3];

specs.chan_names  = {'V','W','Y','Z'}; % First half of electrodes 
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
specs.fft_stim_t  = [0 1.2];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot);
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
specs.fft_stim_t  = [0 1.8];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot);

specs.chan_names  = {'H','M','P','R','S'}; % Second half of electrods
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
specs.fft_stim_t  = [0 1.2];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot);
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
specs.fft_stim_t  = [0 1.8];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot);


%% plot ERP time courses to check *onsets* 

inputFolder       = 'ECoGCAR';
description       = 'reref';

clear specs;

specs.plot_type = 'averageSE';
specs.epoch_t     = [-0.2 0.2];

% specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};

% selected electrodes
specs.chan_names  = {'M5', 'M6', 'M11', 'M12'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

specs.chan_names  = {'P3', 'P5', 'P8', 'P9'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

specs.chan_names  = {'S9','S10','S13','S14'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

% all electrodes
specs.chan_names  = {'H','M','P','R','S'}; % Second half of electrods
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

specs.chan_names  = {'V','W','Y','Z'}; % First half of electrodes 
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

%% plot ERP time courses to check *offsets*

inputFolder       = 'ECoGCAR';
description       = 'reref';

clear specs;
specs.plot_offset = 1;
specs.plot_type   = 'averageSE';
specs.epoch_t     = [-0.2 0.2];

% specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};

% selected electrodes
specs.chan_names  = {'M5', 'M6', 'M11', 'M12'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

specs.chan_names  = {'P3', 'P5', 'P8', 'P9'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

specs.chan_names  = {'S9','S10','S13','S14'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

% all electrodes
specs.chan_names  = {'H','M','P','R','S'}; % Second half of electrods
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

specs.chan_names  = {'V','W','Y','Z'}; % First half of electrodes 
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

% %% Average time series for selected electrodes, each panel per condition
% 
% bidsDir     = tt_bidsRootPath;
% 
% clear specs;
% specs.plot_data    = true;
% specs.plot_smooth  = 1; % could try some other values
% specs.epoch_t      = [-0.4, 1.8]; % stimulus epoch window
% specs.base_t       = [-0.4 -0.1]; % blank epoch window
% specs.chan_names  = {'M'};
% specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
% [data, channels, t, srate]   = tt_prepareData(bidsDir, subject, session, task, [], inputFolder, description, specs);
% 
% specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
% [data, channels, t, srate]   = tt_prepareData(bidsDir, subject, session, task, [], inputFolder, description, specs);
% 
% %%
% % Generate stimulus timecourses
% [stim_ts, stim_info] = tt_generateStimulusTimecourses(specs.stim_names, t);
% 
% specs.average_elecs = true;
% %%
% tt_plotData(data, channels, t, specs, savePlot)






