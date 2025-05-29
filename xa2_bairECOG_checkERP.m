
projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile'; 
subject = 'ny726';
addpath(fullfile(pwd, 'func'))

%% plot ERP time courses to check *onsets* 

savePlot = 1;
session           = 'nyuecog01';
task              = 'temporalpattern';

% Specify one of the three preprocessed data here:
inputFolder       = 'ECoGCAR'; 
description       = 'reref';

clear specs;
specs.epoch_t     = [-0.5 0.5];
specs.plot_ylim   = [-250, 250];
specs.runnums     = 5;
% specs.plot_single_channel = true;
specs.stim_names  = {'ONE-PULSE-6'};

% H13-H16
specs.chan_names  = {'H13','H14','H15','H16'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% M10-M12
specs.chan_names  = {'M10','M11','M12'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% P1-P3
specs.chan_names  = {'P1','P2','P3'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% P8-P10
specs.chan_names  = {'P8','P9','P10'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% S1-S5
specs.chan_names  = {'S1','S2','S3','S4','S5'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% S8-S12
specs.chan_names  = {'S8','S9','S10','S11','S12'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% S13-S15
specs.chan_names  = {'S13','S14','S15'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% Y1-Y4
specs.chan_names  = {'Y1','Y2','Y3','Y4'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% Z12-Z15
specs.chan_names  = {'Z12','Z13','Z14','Z15'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

%% plot offset

specs.plot_offset = true;

% H13-H16
specs.chan_names  = {'H13','H14','H15','H16'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% M10-M12
specs.chan_names  = {'M10','M11','M12'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% P1-P3
specs.chan_names  = {'P1','P2','P3'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% P8-P10
specs.chan_names  = {'P8','P9','P10'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% S1-S5
specs.chan_names  = {'S1','S2','S3','S4','S5'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% S8-S12
specs.chan_names  = {'S8','S9','S10','S11','S12'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% S13-S15
specs.chan_names  = {'S13','S14','S15'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% Y1-Y4
specs.chan_names  = {'Y1','Y2','Y3','Y4'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)

% Z12-Z15
specs.chan_names  = {'Z12','Z13','Z14','Z15'};
bidsEcogPlotChannelAveragedTrials(projectDir, subject, session, task, [], ...
    inputFolder, description, specs, savePlot)
