% Plot DN and LIN model fits of group-averaged, electrode-averaged tactile data

% Model fits of average-electrode are under bidsRootPath/derivatives/modelFit/results

%% 
server_available = true;
if server_available 
    bids_dir     = tt_bidsRootPath;
    data_dir = fullfile(bids_dir, 'derivatives', 'modelFit', 'results');

    fig_dir     = fullfile(bids_dir, 'derivatives', 'modelFit', 'figure', 'check_tradeoff');
    if ~exist(fullfile(fig_dir), 'dir'), mkdir(fullfile(fig_dir)); end
else
    git_dir = fileparts(pwd);
    data_dir = fullfile(git_dir, 'tactileECoGdata');
    fig_dir = fullfile(pwd, 'figures');
end

modelColors = [218, 62, 82; 45, 125, 210]./255; % red and blue

%% Load all data

models_to_plot = {'DN', 'LINEAR'};
str = 'group_average'; %umcudrouwen/ny726/group_average

%% parameter estimates

% Load individual electrode data for parameter comparison
tactile_indiv = load(fullfile(data_dir,sprintf('DN_xvalmode0_individualelecs_%s.mat', str)), ...
                     'data', 'pred', 'stim_info', 'stim', 't', 'channels', 'params'); % Always the full fit to individual electrodes
visual_indiv = load(fullfile(pwd, 'model_results/visual', 'DN_xvalmode0_individualelecs.mat'), ...
                    'data', 'pred', 'stim_info', 'stim', 't', 'channels', 'params');

% Parameter names from DN model
paramNames = {'$\tau_1$', 'w', '$\tau_2$', '$n$', '$\sigma$', 'shift','scale'};
paramSlc = [1:5,7]; % exclude 'shift' parameter
clut = {[1, 0, 0], [0.7, 0.7, 0.7], [0, 0, 0]};

% Process visual electrodes
[~, v_channels, group_prob] = groupElecsByVisualArea(visual_indiv.channels, 'probabilisticresample');
v_nChans = height(v_channels);

% Compute visual area averages using the probabilities
visual_area_means = nan(length(paramSlc), v_nChans);
visual_area_se = nan(length(paramSlc), v_nChans, 2);

for p_idx = 1:length(paramSlc)
    p = paramSlc(p_idx);
    visual_params = visual_indiv.params(p,:);
    [m, se] = averageWithinArea(visual_params, group_prob);
    visual_area_means(p_idx, :) = m;
    visual_area_se(p_idx, :, :) = se;
end

% Compute tactile summary statistics with bootstrap CI
nboot = 10000;
tactile_means = mean(tactile_indiv.params(paramSlc,:), 2);
tactile_ci = nan(length(paramSlc), 2);

% Bootstrap confidence intervals
for n = 1:length(paramSlc)
    p = paramSlc(n);
    params = tactile_indiv.params(p,:);
    boot_means = bootstrp(nboot, @mean, params);
    tactile_ci(n,:) = prctile(boot_means, [15.87 84.13]);
end

%% plot
paramFig = figure('Color', [1 1 1], 'Position', [100 100 1800 600]);
set(paramFig, 'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[1200 600])

% Add tactile electrode groups to x-axis and scatter their parameters
first_letters = cellfun(@(x) x(1), tactile_indiv.channels.name);
unique_groups = unique(first_letters);
n_tactile_groups = numel(unique_groups);

% For positions:
tactile_pos = v_nChans + 1; % Group-averaged tactile
tactile_group_pos = (v_nChans + 2):(v_nChans + 1 + n_tactile_groups);
all_labels = [v_channels.name; {'Tactile'}; cellstr(unique_groups)];

% Gather color for each group from clut
group_colors = parula(n_tactile_groups+1);

% Prepare for use in plotting loop
tactile_group_inds = cell(n_tactile_groups,1);
for g = 1:n_tactile_groups
    tactile_group_inds{g} = find(first_letters == unique_groups(g));
end

for n = 1:numel(paramSlc)
    p = paramSlc(n);

    subplot(2, ceil(numel(paramSlc)/2), n); 
    hold on
    set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out','TickLength', [0.05 0.05]);
    
    % Plot visual areas with error bars
    errorbar(1:v_nChans, visual_area_means(n,:), ...
             visual_area_means(n,:) - visual_area_se(n,:,1), ...
             visual_area_se(n,:,2) - visual_area_means(n,:), ...
             '.', 'Color', 'k', 'MarkerSize', 20, 'LineWidth', 1, ...
             'LineStyle', 'none', 'CapSize', 0);
    
    % Plot tactile average with error bars
    errorbar(tactile_pos, tactile_means(n), ...
             tactile_means(n) - tactile_ci(n,1), ...
             tactile_ci(n,2) - tactile_means(n), ...
             '.', 'Color', 'k', 'MarkerSize', 20, 'LineWidth', 1, ...
             'LineStyle', 'none', 'CapSize', 0);

    errorbar(tactile_pos, tactile_means(n), 0, ...
             '.', 'Color', 'k', 'MarkerSize', 20, 'LineWidth', 1, ...
             'LineStyle', 'none', 'CapSize', 0);

    % Plot tactile groups (right of tactile)
    for g = 1:n_tactile_groups
        these_inds = tactile_group_inds{g};
        these_params = tactile_indiv.params(paramSlc(n), these_inds);
        x_jitter = (rand(size(these_params))-0.5)*0.15;
        scatter(tactile_group_pos(g) + x_jitter, these_params, 55, ...
            'MarkerFaceColor', group_colors(g,:), ...
            'MarkerEdgeColor', 'none', 'LineWidth', 1.2, 'MarkerFaceAlpha', 0.85);
    end

    xlim([0, v_nChans + n_tactile_groups + 2]);
    xticks(1:(v_nChans + 1 + n_tactile_groups));
    xticklabels(all_labels);
    xtickangle(45);

    title(paramNames{p}, 'FontSize', 20, 'Interpreter', 'latex');
    % ylabel('Parameter value', 'FontSize', 20);

    box off

    if p == 1; ylim([0 0.1]); end
        if p == 2; ylim([0 0.8]); end
        if p == 3; ylim([0 0.4]); end
        if p == 4; ylim([0.8 2]); end
        if p == 5; ylim([0 0.15]); end
        if p == 6; ylim([0 0.15]); end
        if p == 7; ylim([0 4]); end
end

% Save figures
saveas(paramFig, fullfile(fig_dir, sprintf('%s_parameters_comparison',str)), 'png');

%% scatter estimates of w vs. tau2

figure;
hold on

w_idx = 2;
tau2_idx = 3;
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out','TickLength', [0.05 0.05]);
scatter(tactile_indiv.params(w_idx,:), tactile_indiv.params(tau2_idx,:), 80, 'k', 'filled');
xlabel(paramNames{w_idx}, 'FontSize', 20, 'Interpreter', 'latex');
ylabel(paramNames{tau2_idx}, 'FontSize', 20, 'Interpreter', 'latex');
box off

saveas(gcf, fullfile(fig_dir, sprintf('%s_scatter_w_vs_tau2',str)), 'png');