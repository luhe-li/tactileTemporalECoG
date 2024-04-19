function [data, channel, epochs] = prepareData(projectDir, subject, sessions, tasks, runnums, ...
    inputFolder, description, specs)

% <projectDir>
if ~exist('projectDir', 'var') || isempty(projectDir)
    error('projectDir not defined');
end

% <subject>
if ~exist('subject', 'var') || isempty(subject)
    error('subject not defined');
end

if ~exist('sessions', 'var'), sessions = []; end
if ~exist('tasks', 'var'), tasks = []; end
if ~exist('runnums', 'var'), runnums = []; end

% <inputFolder>
if ~exist('inputFolder', 'var') || isempty(inputFolder)
    inputFolder = 'ECoGBroadband';
end

% <description>
if ~exist('description', 'var') || isempty(description)
    description = 'broadband';
end

% <specs>
if ~exist('specs', 'var') || isempty(specs), specs = struct(); end
if ~isfield(specs,'epoch_t') || isempty(specs.epoch_t), specs.epoch_t = [-0.2 1.2];end
if ~isfield(specs,'base_t') || isempty(specs.base_t), specs.base_t = [min(specs.epoch_t) 0];end
if ~isfield(specs,'chan_names'), specs.chan_names = []; end
if ~isfield(specs,'stim_names'), specs.stim_names = []; end
if ~isfield(specs,'plot_data'), specs.plot_data = true; end
if ~isfield(specs,'plot_cmap') || isempty(specs.plot_cmap), specs.plot_cmap = 'parula'; end
if ~isfield(specs,'average_channels') || isempty(specs.average_channels), specs.average_channels = false; end

%% load data as bidsEcogPlotTrials

% Load data

dataPath = fullfile(projectDir, 'derivatives', inputFolder);
writePath = fullfile(projectDir, 'derivatives', 'ECoGFigures', subject, 'DNfit');

[data, channels, events] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description, 512);

% Select channels
chan_names = specs.chan_names;

if ~iscell(chan_names), chan_names = {chan_names}; end
if isempty([chan_names{:}])
    chan_idx = contains(lower(channels.type), {'ecog', 'seeg'});

elseif any(contains(chan_names, string(0:10)))
    % specs.chan_names contains numbers, so the user is (probably)
    % referencing individual channels. Match each individual channel:
    chan_idx = ecog_matchChannels(chan_names, channels.name);

else
    % specs.chan_names does not contain any numbers, so the user is
    % (probably) trying to plot a group of channels based on a common
    % character (e.g. G): match all channels with this character:
    % ALTERNATIVE: use channel_group as second way to select channels?
    chan_idx = contains(channels.name, chan_names);
    chan_idx = chan_idx | matches(channels.group, chan_names);

end
if ~any(chan_idx), error('Did not find any matching channels! Please check channel names.'), end

data = data(chan_idx,:);
channels = channels(chan_idx,:);

% Epoch the data
[epochs, t] = ecog_makeEpochs(data, events.onset, specs.epoch_t, channels.sampling_frequency(1));
fprintf('[%s] Found %d epochs across %d runs and %d sessions \n', ...
    mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));

% Baseline correct the epochs
switch description
    case 'reref'
        baseType = 'subtractwithintrial';
    case 'broadband'
        baseType = 'percentsignalchange';
end
fprintf('[%s] Baseline correcting epochs using %s \n', mfilename, baseType);
[epochs] = ecog_normalizeEpochs(epochs, t, specs.base_t, baseType);

% Select trials
fprintf('[%s] Selecting stimulus conditions... \n', mfilename);
stim_names = specs.stim_names;
if ~iscell(stim_names), stim_names = {stim_names}; end
if isempty([stim_names{:}])
    stim_names = unique(events.trial_name);
end
stim_idx = cell(length(stim_names),1);

for ii = 1:length(stim_names)
    if ~isnumeric(stim_names)
        %         stim_idx{ii} = find(contains(events.trial_name, stim_names{ii}));
        stim_idx{ii} = find(strcmp(events.trial_name, stim_names{ii}));
    else
        stim_idx{ii} = find(events.trial_type == stim_names(ii));
    end
    fprintf('[%s] Found %d trials for condition %s \n', mfilename, length(stim_idx{ii}), stim_names{ii});
end
if specs.average_stims, stim_idx = {vertcat(stim_idx{:})}; stim_names = {'all stim averaged'}; end%{[stim_names{:}]};


%% plot data

if specs.plot_data

    % Set plot settings
    if ischar(specs.plot_cmap)
        cmap = eval(specs.plot_cmap);
        colors = cmap(1:round((length(cmap)/length(stim_idx))):length(cmap),:);
    else
        colors = specs.plot_cmap(1:length(stim_idx),:);
    end



end
end