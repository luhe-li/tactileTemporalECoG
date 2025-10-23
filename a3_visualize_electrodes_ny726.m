
clear; close all;
% tbUse tactileTemporalECoG;
addpath(fullfile(pwd, 'utils'));
gitDir = fileparts(pwd);
ft_path = fullfile(gitDir, 'fieldtrip');

addpath(ft_path);
ft_defaults;

projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
subject           = 'ny726';
tasks             = 'temporalpattern';
sessions          = 'nyuecog01';
dataDir           = fullfile(projectDir, 'derivatives');

figureDir = fullfile(projectDir, 'derivatives', 'ECoGFigures', 'sub-ny726','electrodes');
if ~exist(figureDir, 'dir'), mkdir(figureDir); end

%% visualize by electrode positions and default freesurfer template

% Load electrode positions from tsv file
fid = fopen('sub-NY726_electrodes.tsv');
txt = textscan(fid, '%s%f%f%f%s', 'delimiter', '\t', 'Headerlines', 1);
fclose(fid);

% Extract electrode positions
epos = [txt{2} txt{3} txt{4}];
channel_names = txt{1};
first_letters = cellfun(@(x) x(1), channel_names, 'UniformOutput', false);
unique_letters = unique(first_letters);
n_groups = length(unique_letters);

% Create plot vector for electrode groups
plot_vec = zeros(size(epos,1), 1);
for i = 1:length(channel_names)
    letter = channel_names{i}(1);
    group_idx = find(strcmp(unique_letters, letter));
    plot_vec(i) = group_idx;
end

figure('Position', get(0, 'Screensize')); hold on
subplot(131); hold on
plot_ecog(plot_vec, ...
    fullfile(ft_path, 'template/anatomy'), ...
    epos, ...
    [1 n_groups], ... % empty limits for default
    0.1, ... % transparency
    [90 0], ... % view angle
    30, ... % electrode size
    1, ... % plot mesh
    1, ... % plot right hemisphere only
    parula, ... % default colormap
    0); % no colorbar

% Plot labels for first electrode of each group
for i = 1:n_groups
    group_letter = unique_letters{i};
    first_elec_idx = find(strcmp(first_letters, group_letter), 1);
    pos = epos(first_elec_idx,:);
    text(pos(1), pos(2), pos(3), group_letter, 'FontSize', 12, 'FontWeight', 'bold');
end

% Sagittal view
subplot(132); hold on
plot_ecog(plot_vec, ...
    fullfile(ft_path, 'template/anatomy'), ...
    epos, ...
    [1 n_groups], ... % empty limits for default
    0.1, ... % transparency
    [0 0], ... % view angle - sagittal
    30, ... % electrode size
    1, ... % plot mesh
    1, ... % plot right hemisphere only
    parula, ... % default colormap
    0); % no colorbar

% Plot labels for first electrode of each group
for i = 1:n_groups
    group_letter = unique_letters{i};
    first_elec_idx = find(strcmp(first_letters, group_letter), 1);
    pos = epos(first_elec_idx,:);
    % Slightly offset the text above the electrode position
    text(pos(1), pos(2), pos(3)+10, group_letter, 'FontSize', 12, 'FontWeight', 'bold');
end

% Vertical view (top-down)
subplot(133); hold on
plot_ecog(plot_vec, ...
    fullfile(ft_path, 'template/anatomy'), ...
    epos, ...
    [1 n_groups], ... % empty limits for default
    0.1, ... % transparency
    [0 90], ... % view angle - vertical
    30, ... % electrode size
    1, ... % plot mesh
    1, ... % plot right hemisphere only
    parula, ... % default colormap
    0); % no colorbar

% Plot labels for first electrode of each group
for i = 1:n_groups
    % Find first electrode in this group
    group_letter = unique_letters{i};
    first_elec_idx = find(strcmp(first_letters, group_letter), 1);
    
    % Get position of first electrode
    pos = epos(first_elec_idx,:);
    
    % Add text label slightly above the electrode position
    text(pos(1), pos(2), pos(3)+10, group_letter, 'FontSize', 14, 'FontWeight', 'bold');
end

saveas(gcf, fullfile(figureDir, 'electrodes_all'), 'epsc');
