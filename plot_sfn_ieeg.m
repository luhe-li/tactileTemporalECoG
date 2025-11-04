% SFN ecog figures
% 

% Plot DN and LIN model fits of group-averaged, electrode-averaged tactile data

% Model fits of average-electrode are under bidsRootPath/derivatives/modelFit/results

close all
%% 
server_available = true;
if server_available 
    bids_dir     = tt_bidsRootPath;
    data_dir = fullfile(bids_dir, 'derivatives', 'modelFit', 'results');

    fig_dir     = fullfile(bids_dir, 'derivatives', 'modelFit', 'figure', 'sfn');
    if ~exist(fullfile(fig_dir), 'dir'), mkdir(fullfile(fig_dir)); end
else
    git_dir = fileparts(pwd);
    data_dir = fullfile(git_dir, 'tactileECoGdata');
    fig_dir = fullfile(pwd, 'figures');
end

%% Load all data

epoch_t = [-0.4, 1.6];
models_to_plot = {'DN', 'LINEAR'};
str = 'group_average'; %umcudrouwen/ny726/group_average

% colors = parula(4);
% modelColors = colors(2:3,:);

modelColors = [218, 62, 82; 45, 125, 210]./255; % red and blue

%% Evaluate cross-valiedated R2

includeDerivedParams = false;
for ii = 1:length(models_to_plot)

    model = models_to_plot{ii};
    cv_D = load(fullfile(data_dir, sprintf('%s_xvalmode1_electrodeaverages_%s.mat', model, str)));
    cv_results{ii} = tt_evaluateModelFit(cv_D,includeDerivedParams);
    CV_R2(ii) = cv_results{ii}.R2.concat_all;
end

%% Plot timecourses

for dd = 1:length(models_to_plot)

    model = models_to_plot{dd};
    D = load(fullfile(data_dir, sprintf('%s_xvalmode0_electrodeaverages_%s.mat', model, str)));
    data = D.data;
    pred = D.pred;
    stim_info = D.stim_info;
    stim_ts = D.stim;
    t = D.t;
    channels = D.channels;

    % Setup indices for different conditions
    one_idx = find(contains(stim_info.name, 'ONE-PULSE'));
    xDur = stim_info.duration(one_idx);
    pair_idx = find(contains(stim_info.name, {'ONE-PULSE-4', 'TWO-PULSE'}));
    xISI = stim_info.ISI(pair_idx);

    %% Plot timecourses

    tcourseFig = figure('Color', [1 1 1], 'Position', [30 + (dd-1)*750, 300, 900, 400]);
    set(tcourseFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[1000 400])
    T1 = tiledlayout(2, 1,'TileIndexing','rowmajor');

    % Combine pulse conditions
    onePulseIndx = find(contains(stim_info.name, 'ONE-PULSE'));
    twoPulseIndx = cat(1, find(contains(stim_info.name, 'ONE-PULSE-4')), ...
        find(contains(stim_info.name, 'TWO-PULSE')));
    allPulsesIndx = cat(1, onePulseIndx, twoPulseIndx);

    % Set y-bounds
    bounds = @(x) [floor(min(x(:))*10) ceil(max(x(:))*10)]/10;
    ybounds = bounds(prctile(data, [1 99], 'all'));
    if mod(ybounds(1), 0.2) > 0; ybounds(1) = ybounds(1)-0.2; end
    if mod(ybounds(2), 0.2) > 0; ybounds(2) = ybounds(2)+0.5; end

    % Calculate summed response across time series (first dimension)
    sum_data = sum(data, 1);
    sum_pred{dd} = sum(pred, 1);
    
    % Panel for single pulse conditions, with one empty subplot at the beginning
    nSingle = numel(onePulseIndx);
    tpanel = tiledlayout(T1,1,nSingle+1,'TileIndexing','columnmajor');
    tpanel.Layout.Tile = 1;
    tpanel.XLabel.String = 'Time (s)';
    tpanel.Title.String = sprintf('Cross-validated r² = %.3f', CV_R2(dd));

    % Add empty subplot as the first tile
    nexttile(tpanel,1)
    axis off

    for ii = 1:nSingle
        nexttile(tpanel,ii+1)
        hold on
        plot(t, stim_ts(:, onePulseIndx(ii)) * ybounds(end), 'Color', [.5 .5 .5], 'HandleVisibility', 'off')
        plot(t([1 end]), [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
        plot(t, data(:,onePulseIndx(ii)), 'k-', 'LineWidth', 1, 'HandleVisibility', 'off')
        plot(t, pred(:,onePulseIndx(ii)), 'Color', modelColors(dd,:), 'LineWidth', 1)
        set(gca,'TickDir', 'out', 'FontSize', 14, 'XColor', 'k', 'YColor', 'k', ...
            'LineWidth', 1, 'TickLength', [0.05 0.05])
        title(sprintf('Dur %.2fs', xDur(ii)), 'FontSize', 10)
        ylim(ybounds)
        box off
        % set(gca, 'XTickLabel', [], 'YTickLabel', [], 'XLabel', [], 'YLabel', []);
    end

    % Panel for paired pulse conditions
    tpanel2 = tiledlayout(T1,1,numel(twoPulseIndx),'TileIndexing','columnmajor');
    tpanel2.Layout.Tile = 2;
    % tpanel2.XLabel.String = 'Time (s)';
    % tpanel2.Title.String = sprintf('%s - Paired Pulse Conditions', upper(dataset));

    for ii = 1:numel(twoPulseIndx)
        nexttile(tpanel2,[1 1])
        hold on
        plot(t, stim_ts(:, twoPulseIndx(ii)) * ybounds(end), 'Color', [.5 .5 .5], 'HandleVisibility', 'off')
        plot(t([1 end]), [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
        plot(t, data(:,twoPulseIndx(ii)), 'k-', 'LineWidth', 1, 'HandleVisibility', 'off')
        plot(t, pred(:,twoPulseIndx(ii)), 'Color', modelColors(dd,:), 'LineWidth', 1)
        set(gca,'TickDir', 'out', 'FontSize', 14, 'XColor', 'k', 'YColor', 'k', ...
            'LineWidth', 1, 'TickLength', [0.03 0.03])
        title(sprintf('ISI %.2fs', xISI(ii)), 'FontSize', 10)
        ylim(ybounds)
        box off
        % set(gca, 'XTickLabel', [], 'YTickLabel', [], 'XLabel', [], 'YLabel', []);
    end

    % Save figures
    saveas(tcourseFig, fullfile(fig_dir, sprintf('%s_%s_timecourse', str, model)), 'pdf');


end

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

    % Fits of individual electrode
    for g = 1:n_tactile_groups
        these_inds = tactile_group_inds{g};
        these_params = tactile_indiv.params(p, these_inds);
        scatter(tactile_group_pos(g), these_params, 55, ...
            'MarkerFaceColor', group_colors(g,:), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'MarkerFaceAlpha', 0.85);
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
saveas(paramFig, fullfile(fig_dir, sprintf('%s_parameters_comparison',str)), 'pdf');


%% Plot summed responses

summaryFig      = figure('Color', [1 1 1], 'Position', [30 300 800 260]);
set(summaryFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[550 260])
figure(summaryFig)
T = tiledlayout(1, 3,'TileIndexing','rowmajor');

% Load individual electrode data to get confidence interval
tactile_indiv = load(fullfile(data_dir,sprintf('DN_xvalmode0_individualelecs_%s.mat', str)), ...
                     'data', 'pred', 'stim_info', 'stim', 't', 'channels', 'params');

cResponses68 = prctile(squeeze(sum(tactile_indiv.data, 1)), [15.87 84.13], 2);
cResponses95 = prctile(squeeze(sum(tactile_indiv.data, 1)), [2.5, 97.5], 2);

num_one_pulse = numel(onePulseIndx);
num_two_pulse = numel(twoPulseIndx);

xvalues = [0.05, 0.1, 0.2, 0.4, 0.8, 1.2];

    % D = load(fullfile(data_dir, sprintf('%s_xvalmode0_electrodeaverages_%s.mat', model, str)));
    % data = D.data;
    % pred = D.pred;
    % stim_info = D.stim_info;
    % stim_ts = D.stim;
    
    % One pulse conditions
    t1 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
    t1.Layout.Tile = 1;
    % ax1 = nexttile(t1,[1 1]);
    % hold on,
    % plot([-0.1 0.1], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    % plot(0 .* [1; 1], cResponses95(1,:), 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
    % plot(0 .* [1; 1], cResponses68(1,:), 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
    % plot(0,  sum_data(1), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

    ax2 = nexttile(t1,[1 2]);
    hold on
    plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    plot(xvalues .* [1; 1], cResponses95(1:num_one_pulse,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
    plot(xvalues .* [1; 1], cResponses68(1:num_one_pulse,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
    plot(xvalues,  sum_data(1:num_one_pulse), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

    % Paired pulse conditions
    t2 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
    t2.Layout.Tile = 2;
    ax3 = nexttile(t2,[1 1]);
    hold on
    plot([-0.1 0.1], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    plot(0 .* [1; 1], cResponses95(twoPulseIndx(1),:), 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
    plot(0 .* [1; 1], cResponses68(twoPulseIndx(1),:), 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
    plot(0,  sum_data(twoPulseIndx(1)), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')
    
    ax4 = nexttile(t2,[1 2]);
    hold on
    plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    plot(xvalues .* [1; 1], cResponses95(twoPulseIndx(2:end),:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
    plot(xvalues .* [1; 1], cResponses68(twoPulseIndx(2:end),:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
    plot(xvalues,  sum_data(twoPulseIndx(2:end)), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

    % Linear xaxis
    t3 = tiledlayout(T,2,1,'TileIndexing','columnmajor');
    t3.Layout.Tile = 3;
    ax5 = nexttile(t3,[1 1]);
    hold on
    plot([-0.05 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    plot(xvalues .* [1; 1], cResponses95(onePulseIndx,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
    plot(xvalues .* [1; 1], cResponses68(onePulseIndx,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
    plot(xvalues,  sum_data(onePulseIndx), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

    ax6 = nexttile(t3,[1 1]);
    hold on
    plot([-0.05 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    plot([0 xvalues].* [1; 1], cResponses95(twoPulseIndx,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
    plot([0 xvalues] .* [1; 1], cResponses68(twoPulseIndx,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
    plot([0 xvalues],  sum_data(twoPulseIndx), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')


    for dd = 1:length(models_to_plot)

        plot(ax2, xvalues, sum_pred{dd}(1:num_one_pulse), 'Color', modelColors(dd,:), 'LineWidth', 1.5)
        plot(ax3, 0, sum_pred{dd}(twoPulseIndx(1)), 'Color', modelColors(dd,:), 'LineWidth', 1.5)
        plot(ax4, xvalues, sum_pred{dd}(twoPulseIndx(2:end)), 'Color', modelColors(dd,:), 'LineWidth', 1.5)
        plot(ax5, xvalues, sum_pred{dd}(onePulseIndx), 'Color', modelColors(dd,:), 'LineWidth', 1.5 )
        plot(ax6, [0 xvalues], sum_pred{dd}(twoPulseIndx), 'Color', modelColors(dd,:), 'LineWidth', 1.5 )

    end

    ybounds = get(ax2, 'YLim');
    
    ax2.Box = 'off';
    t1.Title.String = 'Single-pulse conditions';
    t1.XLabel.String = 'Stimulus duration (s)';
    t1.YLabel.String = {'Summed broadband time-series';'(X-fold)'};
    

    set(ax2, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
        'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
        'YTick', ybounds(1):200:ybounds(2), 'YLim', ybounds, 'XLim', [0.035 1.3]);
    ax2.Box = 'off';
    % ax2.YAxis.Visible = 'off';
  
    set(ax3, 'TickDir', 'out', 'XTick', 0, 'XTickLabel', 0, ...
        'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
        'YTick', ybounds(1):200:ybounds(2), 'YLim', ybounds, 'XLim', [-0.005 0.01]);
    ax3.Box = 'off';
    t2.Title.String = 'Paired-pulse conditions';
    t2.XLabel.String = 'Interstimulus interval (s)';
    t2.YLabel.String = {'Summed broadband time-series';'(X-fold)'};
    
    set(ax4, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
        'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
        'YTick', ybounds(1):200:ybounds(2), 'YLim', ybounds, 'XLim', [0.04 1.3]);
    ax4.Box = 'off';
    ax4.YAxis.Visible = 'off';
     

    set(ax5, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
        'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
        'YTick', ybounds(1):200:ybounds(2), 'YLim', ybounds, 'XLim', [-0.1 1.3]);
    ax5.Box = 'off';
    ax5.Title.String = 'Single-pulse conditions';
    ax5.XLabel.String = 'Duration (s)';
    ax5.DataAspectRatio = [1 1200 1];
    t3.YLabel.String = {'Summed broadband time-series';'(X-fold)'};

    set(ax6, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
        'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
        'YTick', ybounds(1):200:ybounds(2), 'YLim', ybounds, 'XLim', [-0.1 1.3]);
    ax6.Box = 'off';
    ax6.Title.String = 'Paired-pulse conditions';
    ax6.DataAspectRatio = [1 1200 1];
    ax6.XLabel.String = 'ISI (s)';

    saveas(summaryFig, fullfile(fig_dir, sprintf('%s_summed_responses', str)), 'pdf');

 