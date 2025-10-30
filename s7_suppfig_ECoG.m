
% Dependencies

% Model fits of both average-electrodes and individual-electrodes under dataRootDir/ECoG/visual and dataRootDir/ECoG/tactile
% ECoG_utils: groupElecsByVisualArea, averageWithinArea, bootstrWithinArea
% Functions under Code/ECoG/

%% A-B. Plot cross-validated average-electrodes data with DN prediction

% Path

datasets = {'tactile','visual'};
dataRootDir = fullfile(pwd, '/Code/ECoG/data/');
dataStr = [];

% datasets = {'tactile', 'visual'};
% dataStr         = 'ECoG';
% [~, dataRootDir]= rootPath(false);
figRoot         = fullfile(dataRootDir, '..', 'Figures');
figDir          = fullfile(figRoot, dataStr);
if ~exist(fullfile(figDir), 'dir'), mkdir(fullfile(figDir)); end

% modelColors = [30, 120, 180; 227, 27, 27]./255;
modelColors = [0,0,0;77, 118, 255]./255;
allSumData = [];
allSumPred = [];


% Collect all summed data for y-axis
for dd = 1:length(datasets)
    dataset = datasets{dd};
    fits = load(fullfile(dataRootDir, dataStr, dataset, 'DN_xvalmode1_electrodeaverages.mat'), ...
                'data', 'pred', 'stim_info', 'stim', 't', 'channels');
    sumResp_data = squeeze(sum(fits.data,1));
    sumResp_pred = squeeze(sum(fits.pred,1));
    allSumData = [allSumData, sumResp_data(:)'];
    allSumPred = [allSumPred, sumResp_pred(:)'];
end
allData = [allSumData, allSumPred];
globalSumBounds = [floor(min(allData))-0.5, ceil(max(allData))+0.5];

for dd = 1:length(datasets)
    dataset = datasets{dd};
    
    % Load electrode averages data
    fits = load(fullfile(dataRootDir, dataStr, dataset, 'DN_xvalmode1_electrodeaverages.mat'), ...
                'data', 'pred', 'stim_info', 'stim', 't', 'channels');
    
    data = fits.data;
    pred = fits.pred;
    stim_info = fits.stim_info;
    stim_ts = fits.stim;
    t = fits.t;
    
    % Compute concatenated R2
    DATA = data(:);
    MODEL = pred(:);
    R2concat = 1 - sum((DATA-MODEL).^2) ./ sum((DATA-mean(DATA)).^2);
    fprintf('%s R2 = %.3f\n', dataset, R2concat);
    
    % Setup indices for different conditions
    one_idx = find(contains(stim_info.name, 'ONEPULSE'));
    xDur = stim_info.duration(one_idx);
    pair_idx = find(contains(stim_info.name, {'ONEPULSE-4', 'TWOPULSE'}));
    xISI = stim_info.ISI(pair_idx);

    %% Time course figure

    tcourseFig = figure('Color', [1 1 1], 'Position', [30 + (dd-1)*750, 300, 900, 400]);
    set(tcourseFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[700 400])
    T1 = tiledlayout(2, 1,'TileIndexing','rowmajor');

    % Combine pulse conditions
    onePulseIndx = find(contains(stim_info.name, 'ONEPULSE'));
    twoPulseIndx = cat(1, find(contains(stim_info.name, 'ONEPULSE-4')), ...
        find(contains(stim_info.name, 'TWOPULSE')));
    allPulsesIndx = cat(1, onePulseIndx, twoPulseIndx);

    % Set y-bounds
    bounds = @(x) [floor(min(x(:))*10) ceil(max(x(:))*10)]/10;
    ybounds = bounds(prctile(data, [1 99], 'all'));
    if mod(ybounds(1), 0.2) > 0; ybounds(1) = ybounds(1)-0.2; end
    if mod(ybounds(2), 0.2) > 0; ybounds(2) = ybounds(2)+0.5; end

    if dd == 2; ybounds = [-0.5, 20]; end

    Responses = NaN(numel(allPulsesIndx), numel(t));
    % No bootstrapped CI in this version, so just use NaNs
    ciResponses = NaN(numel(allPulsesIndx), numel(t), 2);

    % Panel for single pulse conditions, with one empty subplot at the beginning
    nSingle = numel(onePulseIndx);
    tpanel = tiledlayout(T1,1,nSingle+1,'TileIndexing','columnmajor');
    tpanel.Layout.Tile = 1;
    tpanel.XLabel.String = 'Time (s)';
    tpanel.Title.String = sprintf('%s - Single Pulse Conditions (R² = %.3f)', upper(dataset), R2concat);

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
    tpanel2.Title.String = sprintf('%s - Paired Pulse Conditions', upper(dataset));

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
    saveas(tcourseFig, fullfile(figDir, sprintf('timecourse_%s', dataset)), 'pdf');
end

%% D. Parameter estimates comparison
% Load individual electrode data for parameter comparison
tactile_indiv = load(fullfile(dataRootDir, dataStr, 'tactile', 'DN_xvalmode0_individualelecs.mat'), ...
                     'data', 'pred', 'stim_info', 'stim', 't', 'channels', 'params');
visual_indiv = load(fullfile(dataRootDir, dataStr, 'visual', 'DN_xvalmode0_individualelecs.mat'), ...
                    'data', 'pred', 'stim_info', 'stim', 't', 'channels', 'params');

% Parameter names from DN model
% paramNames = {'$\tau_1$', 'w', '$\tau_2$', '$n$', '$\sigma$', 'shift', 'scale'};
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

paramFig = figure('Color', [1 1 1], 'Position', [100 100 1200 600]);
set(paramFig, 'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[1200 600])


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
    tactile_pos = v_nChans + 1;
    errorbar(tactile_pos, tactile_means(n), ...
             tactile_means(n) - tactile_ci(n,1), ...
             tactile_ci(n,2) - tactile_means(n), ...
             '.', 'Color', 'k', 'MarkerSize', 20, 'LineWidth', 1, ...
             'LineStyle', 'none', 'CapSize', 0);
    
    all_labels = [v_channels.name; {'Tactile'};];
    
    xlim([0, v_nChans + 2]);
    xticks(1:(v_nChans + 1));
    xticklabels(all_labels);
    xtickangle(45);
    
    % title(paramNames{p}, 'FontSize', 20, 'Interpreter', 'latex');
    % ylabel('Parameter value', 'FontSize', 20);
    
    box off

    if p == 1; ylim([0 0.1]); end
        if p == 2; ylim([0 0.6]); end
        if p == 3; ylim([0 0.4]); end
        if p == 4; ylim([0.8 1.8]); end
        if p == 5; ylim([0 0.1]); end
        if p == 6; ylim([0 0.15]); end
        if p == 7; ylim([0 4]); end

end

% Save figures
saveas(paramFig, fullfile(figDir, 'parameters_comparison'), 'pdf');
