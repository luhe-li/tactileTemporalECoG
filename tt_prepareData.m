function [data, channels, epochs, t, data_slc, ci_slc] = tt_prepareData(projectDir, subject, sessions, tasks, runnums, ...
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
if ~isfield(specs,'save_plot'), specs.save_plot = true; end
if ~isfield(specs,'epoch_t') || isempty(specs.epoch_t), specs.epoch_t = [-0.2 1.2];end
if ~isfield(specs,'base_t') || isempty(specs.base_t), specs.base_t = [min(specs.epoch_t) 0];end
if ~isfield(specs,'chan_names'), specs.chan_names = []; end
if ~isfield(specs,'stim_names'), specs.stim_names = []; end
if ~isfield(specs,'plot_smooth'), specs.plot_smooth = 0; end
if ~isfield(specs,'plot_data'), specs.plot_data = true; end
if ~isfield(specs,'plot_cmap') || isempty(specs.plot_cmap), specs.plot_cmap = 'parula'; end
if ~isfield(specs,'average_stims') || isempty(specs.average_stims), specs.average_stims = false; end

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
%- check if need to convert run_name into numbers
if iscell(events.run_name), run_name = cell2mat(events.run_name); end
fprintf('[%s] Found %d epochs across %d runs and %d sessions \n', ...
    mfilename, size(epochs,2), length(unique(run_name)), length(unique(events.session_name)));

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

data_slc = nan(length(t), length(stim_names), numel(chan_idx));
ci_slc   = nan(length(t), 2, length(stim_names), numel(chan_idx));

% Loop over electrodes
for ee = 1:numel(chan_idx)

    if ~isnan(chan_idx(ee))

        % Loop over trial types
        for ss = 1:length(stim_idx)

            this_epoch = epochs(:,stim_idx{ss}, ee);
            this_trial = mean(this_epoch,2, 'omitnan');
            llim = this_trial - std(this_epoch,0,2, 'omitnan')/sqrt(length(stim_idx{ss}));
            ulim = this_trial + std(this_epoch,0,2, 'omitnan')/sqrt(length(stim_idx{ss}));
            CI = [llim ulim];

            % Plot
            if specs.plot_smooth > 0
                this_trial = smooth(this_trial,specs.plot_smooth);
                if ~isempty(CI)
                    CI(:,1) = smooth(CI(:,1),specs.plot_smooth);
                    CI(:,2) = smooth(CI(:,2),specs.plot_smooth);
                end
            end

            % Collect data in output
            data_slc(:,ss,ee) = this_trial;
            ci_slc(:,:,ss,ee) = CI; % CI between trials
        end

    end
end

if specs.plot_data

    figureName = sprintf('selecteddata_bystimulus_allelectrodes');
    FontSz = 20;
    FigSz = [400 200 2000 1200];

    %     % Plot each condition as a separate timecourse, all electrodes superimposed
    figure('Name', figureName); hold on;
    colors = parula(height(channels)+1);
    for ii = 1:height(channels), plot(flatten(data_slc(:,:,ii)), 'LineWidth', 2, 'Color', colors(ii,:));end
    set(gca, 'xtick', (0:length(stim_names))*size(data_slc,1)+1, 'xgrid', 'on', 'xticklabel', stim_names, 'xticklabelrotation', 45);
    axis tight
    set(gcf, 'Position', FigSz);
    set(gca, 'FontSize', FontSz);
    legend(channels.name);
    ylabel('x-fold increase in broadband power');

    if specs.save_plot 
        saveDir = fullfile(pwd, 'figures', 'data');
        if ~exist(saveDir, 'dir'), mkdir(saveDir);end
        fprintf('[%s] Saving figures to %s \n',mfilename, saveDir);
        saveas(gcf, fullfile(saveDir, figureName), 'png'); close;
    end

    % Plot each condition as a separate timecourse, all electrodes
    % averaged
    figureName = sprintf('selecteddata_bystimulus_electrodeaverage');
    figure('Name', figureName); hold on;
    plot(flatten(mean(data_slc, 3, 'omitnan')), 'LineWidth', 2, 'Color', 'k');
    set(gca, 'xtick', (0:length(stim_names))*size(data_slc,1)+1, 'xgrid', 'on', 'xticklabel', stim_names, 'xticklabelrotation', 45);
    axis tight
    set(gcf, 'Position', FigSz);
    set(gca, 'FontSize', FontSz);
    ylabel('x-fold increase in broadband power');

    % save Plot?
    if specs.save_plot 
        saveDir = fullfile(pwd, 'figures', 'data');
        if ~exist(saveDir, 'dir'), mkdir(saveDir);end
        fprintf('[%s] Saving figures to %s \n',mfilename, saveDir);
        saveas(gcf, fullfile(saveDir, figureName), 'png'); close;
    end

end

