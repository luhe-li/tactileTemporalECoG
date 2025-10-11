% This script follows a0 and takes in preprocessed data (remoevd bad channel, epoched, converted to BIDS)

clear;
% tbUse tactileTemporalECoG

projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
subject = 'ny726';
sessions = 'nyuecog01';
outputFolder = 'group_ref';
tasks = 'temporalpattern';
runnums = [];

%% apply group-specific reference
[session, tasks, runnums] = bidsSpecifyData(projectDir, subject, sessions, tasks, runnums);
fprintf('[%s] Starting group-averaging for subject: %s session: %s\n', mfilename, subject, session);

% define paths
writePath = fullfile(projectDir, 'derivatives', outputFolder);


for jj = 1:length(tasks)
    for kk = 1:length(runnums{jj})

        task = tasks{jj};
        runnum = runnums{jj}{kk};
        fprintf('[%s] Task = %s, Run = %s \n', mfilename, task, runnum);

        [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(projectDir, subject, session, task, runnum);

        % Overwrite channels sampling frequency with header info
        channels.sampling_frequency = repmat(hdr.Fs, [height(channels) 1]);

        %% Apply group-specific reference

        % Get unique electrode groups from first letter of channel names
        groups = unique(cellfun(@(x) x(1), channels.name, 'UniformOutput', false));

        % Initialize
        data_reref = data;
        channels_reref = channels;
        group_indices = cell(size(groups)-2); % skip DC and EKG channels

        % Loop through each electrode group
        for g = 3:length(groups) % skip DC and EKG channels

            % Find channels belonging to this group
            chan_index = startsWith(channels.name, groups{g});

            % Only include good channels in the reference
            good_channels = find(contains(channels(chan_index,:).status, 'good'));

            % The mean of the good channels is regressed out of all channels
            % (i.e. including the bad ones)
            temp = ecog_carRegress(data(chan_index,:), good_channels);
            data_reref(chan_index,:) = temp;
            channels_reref.reference(chan_index,:) = {'group_ref'};
        end

        % In the ieeg.json file, update the iEEGreference field
        ieeg_json.iEEGReference = 'A group-specific reference was computed separately for each group of good channels and regressed out of each channel.';

        % Update the description and save out the data
        [fname_out] = bidsEcogWriteFiles(writePath, subject, session, task, runnum, 'reref', ...
            data_reref, channels_reref, events, ieeg_json, hdr);

    end
end

%%
% %% extract/check broadband (do it once)
% inputFolder = 'group_ref';
% outputFolder      = 'groupRef_ECoGBroadband_exclude110Hz';
% bands             = [[70 80]; [80 90]; [90 100]; [130 140]; [140 150]; [150 160]; [160 170]];
% bidsEcogBroadband(projectDir, subject, [], [], [], bands, [], inputFolder, outputFolder);

%% plot broadband timecourses for temporal conditions

savePlot = 1;
session           = 'nyuecog01';
task              = 'temporalpattern';

% Specify one of the three preprocessed data here:
inputFolder       = 'groupRef_ECoGBroadband_exclude110Hz'; 
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
