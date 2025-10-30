% Plot DN and LIN model fits of group-averaged, electrode-averaged tactile data

% Model fits of average-electrode are under bidsRootPath/derivatives/modelFit/results

%% 
server_available = false;
if server_available 
    bids_dir     = tt_bidsRootPath;
    data_dir = fullfile(bids_dir, 'derivatives', 'modelFit', 'results');
    str         = 'group_average';

    figRoot     = fullfile(bids_dir, 'derivatives', 'modelFit', 'figure', str);
    if ~exist(fullfile(figRoot), 'dir'), mkdir(fullfile(figRoot)); end
else
    git_dir = fileparts(pwd);
    data_dir = fullfile(git_dir, 'tactileECoGdata');
    fig_dir = fullfile(pwd, 'figures');
end

modelColors = [218, 62, 82; 45, 125, 210]./255; % red and blue


%% Load all data

models_to_plot = {'DN', 'LINEAR'};
str = 'ny726';

%% Evaluate cross-valiedated R2

includeDerivedParams = false;
for ii = 1:length(models_to_plot)

    model = models_to_plot{ii};
    cv_D = load(fullfile(data_dir, sprintf('%s_xvalmode1_electrodeaverages_%s.mat', model, str)));
    cv_results{ii} = tt_evaluateModelFit(cv_D,includeDerivedParams);
    CV_R2(ii) = cv_results{1}.R2.concat_all;
end

%% Plot full-fit data

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
    set(tcourseFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[700 400])
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

    % if dd == 2; ybounds = [-0.5, 20]; end

    Responses = NaN(numel(allPulsesIndx), numel(t));
    % No bootstrapped CI in this version, so just use NaNs
    ciResponses = NaN(numel(allPulsesIndx), numel(t), 2);

    % Panel for single pulse conditions, with one empty subplot at the beginning
    nSingle = numel(onePulseIndx);
    tpanel = tiledlayout(T1,1,nSingle+1,'TileIndexing','columnmajor');
    tpanel.Layout.Tile = 1;
    tpanel.XLabel.String = 'Time (s)';
    tpanel.Title.String = sprintf('R² = %.3f', CV_R2(dd));

    % Add empty subplot as the first tile
    nexttile(tpanel,1)
    axis off

    for ii = 1:nSingle
        nexttile(tpanel,ii+1)
        hold on
        plot(t, stim_ts(:, onePulseIndx(ii)) * ybounds(end), 'Color', [.5 .5 .5], 'HandleVisibility', 'off')
        plot(t([1 end]), [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
        plot(t, data(:,onePulseIndx(ii)), 'k-', 'LineWidth', 1, 'HandleVisibility', 'off')
        plot(t, pred(:,onePulseIndx(ii)), 'Color', modelColors(2,:), 'LineWidth', 1)
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
        plot(t, pred(:,twoPulseIndx(ii)), 'Color', modelColors(2,:), 'LineWidth', 1)
        set(gca,'TickDir', 'out', 'FontSize', 14, 'XColor', 'k', 'YColor', 'k', ...
            'LineWidth', 1, 'TickLength', [0.03 0.03])
        title(sprintf('ISI %.2fs', xISI(ii)), 'FontSize', 10)
        ylim(ybounds)
        box off
        % set(gca, 'XTickLabel', [], 'YTickLabel', [], 'XLabel', [], 'YLabel', []);
    end

    % Save figures
    % saveas(tcourseFig, fullfile(fig_dir, sprintf('timecourse_%s_%s', model, str)), 'pdf');

end