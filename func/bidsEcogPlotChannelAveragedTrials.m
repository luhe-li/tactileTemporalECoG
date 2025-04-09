function [out] = bidsEcogPlotChannelAveragedTrials(projectDir, subject, sessions, tasks, runnums, ...
    inputFolder, description, specs, savePlot)

% Plots all single-trial ERPs for a single condition across all channels.
%
% Input
%     projectDir:       path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     sessions:         BIDS session name (string, all lower case)
%     tasks:            one or more BIDS tasks (string or cell array of strings)
%     runnums:          BIDS run numbers (vector or cell array of vectors)
%     inputFolder:      Name of a derivatives folder where broadband is
%                       to be computed on (e.g., rereferenced data)
%     description:      String stating the 'desc-' label in the name of
%                       the input data files
%     specs:            A struct specifying instructions for what to plot,
%                       with the following possible fields:
%                        - epoch_t: Epoch window (s)
%                                           default: [-0.2 1.2]
%                        - base_t:  Baseline time window (s)
%                                           default: all t<0 in epochTime.
%                        - chan_names: Cell array with channel names 
%                                           default: all channels in channels table                           
%                        - stim_name: Single stimulus name (string)
%                                           default: first unique name in events table                         
%                        - plot_ylim: limits of y axis e.g. [-1 10]
%                                           default: automatic by matlab
%                        - plot_smooth: extent of smoothing of time course
%                                           default: 0
%                        - plot_single_channel: Flag to plot each channel
%                          as a subplot with single trials
%                                           default: false
%     savePlot:         Flag indicating whether to save the plots in 
%                       derivatives/ECoGfigures
%                           default: true

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
if ~isfield(specs,'stim_name') || isempty(specs.stim_name), specs.stim_name = []; end
if ~isfield(specs,'plot_ylim'), specs.plot_ylim = []; end
if ~isfield(specs,'plot_smooth'), specs.plot_smooth = 0; end
if ~isfield(specs,'plot_single_channel'), specs.plot_single_channel = false; end
if ~isfield(specs,'plot_offset'), specs.plot_offset = false; end

% <plot save>
if ~exist('savePlot', 'var') || isempty(savePlot), savePlot = true; end

%% PREPARE DATA

% Load data
dataPath = fullfile(projectDir, 'derivatives', inputFolder);
writePath = fullfile(projectDir, 'derivatives', 'ECoGFigures');
[data, channels, events] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description, 512);
% Select channels
chan_names = specs.chan_names;
if ~iscell(chan_names), chan_names = {chan_names}; end
if isempty([chan_names{:}]) 
    chan_idx = contains(lower(channels.type), {'ecog', 'seeg'}); 
else
    chan_idx = false(height(channels), 1);
    for i = 1:length(chan_names)
        exact_match = strcmp(channels.name, chan_names{i});
        chan_idx = chan_idx | exact_match;
    end
end
if ~any(chan_idx), error('Did not find any matching channels! Please check channel names.'), end

data = data(chan_idx,:);
channels = channels(chan_idx,:);

% Epoch the data
if specs.plot_offset % epoch by trial offset

    event_offset = zeros(height(events), 1);

    % Convert trial_name to string if needed
    if iscell(events.trial_name)
        trial_names = string(events.trial_name);
    else
        trial_names = events.trial_name;
    end

    % Identify 'ONE-PULSE' and 'TWO-PULSE' trials
    is_one_pulse = contains(trial_names, "ONE-PULSE");
    is_two_pulse = contains(trial_names, "TWO-PULSE");

    % Compute offsets
    event_offset(is_one_pulse) = events.onset(is_one_pulse) + events.duration(is_one_pulse);
    event_offset(is_two_pulse) = events.onset(is_two_pulse) + 2 * events.duration(is_two_pulse) + events.ISI(is_two_pulse);

    % Add to the table
    events.event_offset = event_offset;
    [epochs, t] = ecog_makeEpochs(data, events.event_offset, specs.epoch_t, channels.sampling_frequency(1));
    fprintf('[%s] Found %d epochs across %d runs and %d sessions \n', ...
        mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));

    % Add offset to figure filename and figure title
    specs.str = 'offset';

else % epoch by trial onset
    [epochs, t] = ecog_makeEpochs(data, events.onset, specs.epoch_t, channels.sampling_frequency(1));
    fprintf('[%s] Found %d epochs across %d runs and %d sessions \n', ...
        mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));
    specs.str = 'onset';
end

% Baseline correct the epochs
switch description
    case 'reref'
        baseType = 'subtractwithintrial';
    case 'broadband'
        baseType = 'percentsignalchange';
end
[epochs] = ecog_normalizeEpochs(epochs, t, specs.base_t, baseType);
% Select trials for the specified condition
stim_name = specs.stim_names;
if isempty(stim_name)
    stim_name = unique(events.trial_name{1});
end
stim_idx = find(contains(events.trial_name, stim_name));
if isempty(stim_idx), error('No trials found for the specified condition.'); end
if iscell(stim_name)
    str_stim_name = stim_name{1}; % Take the first element if it's a cell
end

%% PLOT DATA

if specs.plot_single_channel
    % Plot each channel as a subplot
    num_channels = size(epochs, 3);
    num_rows = ceil(sqrt(num_channels));
    num_cols = ceil(num_channels / num_rows);
    
    figure('Name', sprintf('Single-trial ERPs for %s (Each Channel)', str_stim_name), ...
        'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    
    for ch = 1:num_channels
        subplot(num_rows, num_cols, ch);
        hold on;
        
        % Plot each trial for the current channel
        colors = parula(length(stim_idx)); 
        for trial = 1:length(stim_idx)
            this_trial = epochs(:, stim_idx(trial), ch);
            if specs.plot_smooth > 0
                this_trial = smooth(this_trial, specs.plot_smooth);
            end
            plot(t, this_trial, 'Color', colors(trial,:),'LineWidth', 1);
        end
        
        title(channels.name{ch}, 'Interpreter', 'none');
        if ~isempty(specs.plot_ylim)
            ylim(specs.plot_ylim);
        end
        if ch > (num_rows - 1) * num_cols
            xlabel('Time (s)');
        else
            set(gca, 'XTickLabel', []);
        end
        if mod(ch, num_cols) == 1
            ylabel(sprintf('%s (mV)', description));
        else
            set(gca, 'YTickLabel', []);
        end
        set(gca, 'FontSize', 10);
    end
else
    % Average time series across channels
    avg_epochs = mean(epochs, 3);

    figure('Name', sprintf('Single-trial ERPs for %s (Channels: %s)', str_stim_name, strjoin(channels.name, ', ')), ...
        'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    hold on;

    % Add dashed lines to mark x = 0 and y = 0
    xline(0, '--k', 'LineWidth', 1); % Dashed line at x = 0
    yline(0, '--k', 'LineWidth', 1); % Dashed line at y = 0

    % Plot each trial with a colormap from parula
    colors = parula(length(stim_idx)); % Generate a colormap using parula
    for trial = 1:length(stim_idx)
        this_trial = avg_epochs(:, stim_idx(trial));
        if specs.plot_smooth > 0
            this_trial = smooth(this_trial, specs.plot_smooth);
        end
        plot(t, this_trial, 'Color', colors(trial, :), 'LineWidth', 1, ...
            'DisplayName', sprintf('Trial %d', trial)); % Use parula colors for each trial
    end
    % legend()

    xlabel(['Time (s), 0 = ' specs.str]);
    ylabel(sprintf('%s (mV)', description));
    title(sprintf('Single-trial ERPs for %s (averaged over %d channels: %s)', ...
        str_stim_name, size(epochs, 3), strjoin(channels.name, ', ')));
    if ~isempty(specs.plot_ylim)
        ylim(specs.plot_ylim);
    end
    set(gca, 'FontSize', 12);
end

%% Save Plot
if savePlot
    figSaveDir = fullfile(writePath, sprintf('sub-%s', subject), 'figures', 'singletrial');
    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir); 
    end
    if specs.plot_single_channel
        saveas(gcf, fullfile(figSaveDir, sprintf('sub-%s_%s_indiv_%s_channels-%s_singletrial.png', ...
            subject, specs.str, str_stim_name, strjoin(channels.name, '-'))), 'png');
    else
        saveas(gcf, fullfile(figSaveDir, sprintf('sub-%s_%s_%s_channels-%s_singletrial.png', ...
        subject, specs.str, str_stim_name, strjoin(channels.name, '-'))), 'png');
    end
end

out.t = t;
out.epochs = epochs;
out.channels = channels;
out.stim_name = str_stim_name;

end
