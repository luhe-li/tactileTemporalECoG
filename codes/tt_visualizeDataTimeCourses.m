
% replicating Iris' script 'bairECOG_analyze_umcudrouwen.m'

projectDir      = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile'; 
subject         = 'umcudrouwen';
recompute       = 0;
savePlot        = 1;
session          = 'umcuiemu01';
task             = 'vtstemporalpattern';

%% common average reference
if recompute

%     bidsEcogRereference(projectDir, subject);

%     outputFolder      = 'ECoGBroadband';
%     bidsEcogBroadband(projectDir, subject, [], task, [], [], [], [], outputFolder, savePlot);

    outputFolder      = 'ECoGBroadband_include110Hz';
    bands             = [[60 70]; [70 80]; [80 90]; [110 120]; [120 130]; [130 140]; [160 170]; [170 180]; [180 190]];
    bidsEcogBroadband(projectDir, subject, [], [], [], bands, [], [], outputFolder);

    outputFolder      = 'ECoGBroadband_exclude110Hz';
    bands             = [[60 70]; [70 80]; [80 90]; [120 130]; [130 140]; [160 170]; [170 180]; [180 190]];
    bidsEcogBroadband(projectDir, subject, [], [], [], bands, [], [], outputFolder);

    outputFolder      = 'ECoGBroadband_110Hz';
    bands             = [105 115];
    method            = @(bp) abs(hilbert(bp)).^2; % geomean should be eliminated when filtering in single band
    bidsEcogBroadband(projectDir, subject, [], [], [], bands, method, [], outputFolder);

end

% plot spectra for temporal conditions
clear specs;

% Specify one of the three preprocessed data here:
inputFolder       = 'ECoGCAR'; 
description       = 'reref';

specs.epoch_t     = [-0.4 1.8];%[-1 2]; % stimulus epoch window
specs.fft_blank_t = [-0.4 -0.1]; %[-0.8 -0.1]; % fft blank epoch window
specs.fft_w       = 0.4;
specs.fft_ov      = 0.1;
specs.plot_ylim   = [10^-3 10^3];
specs.fft_stim_t  = [0 1.2];

% all channels on the grid
specs.chan_names  = {'C'};
specs.subplotdims = [4 8];
specs.subplotidx  = [32:-1:1];
specs.plot_type   = 'average';

specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
bidsEcogPlotSpectra(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); 

specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
bidsEcogPlotSpectra(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); 

% subsets of channels
specs.subplotdims = [2 2];
specs.subplotidx  = [4:-1:1];

specs.chan_names  = {'C03', 'C04', 'C11', 'C12'};

specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
specs.fft_stim_t  = [0 1.2];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], [], [], specs, savePlot); %close;

specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
specs.fft_stim_t  = [0 1.8];
bidsEcogPlotSpectra(projectDir, subject, session, task, [], [], [], specs, savePlot); %close;

%% plot broadband timecourses for temporal conditions

% Specify one of the three preprocessed data here:
inputFolder       = 'ECoGBroadband_exclude110Hz'; 
description       = 'broadband';

clear specs;

specs.epoch_t     = [-0.4 1.8]; % stimulus epoch window
specs.base_t      = [-0.4 -0.1]; % blank epoch window
specs.plot_ylim   = [-2 10];

% all channels on the grid
specs.chan_names  = {'C'};
specs.subplotdims = [4 8];
specs.subplotidx  = [32:-1:1];
specs.plot_type   = 'average';

specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close

specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close

% subsets of channels
specs.subplotdims = [2 2];
specs.subplotidx  = [4:-1:1];

specs.chan_names  = {'C03', 'C04', 'C11', 'C12'};

specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close;

specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot);