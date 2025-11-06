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

    fig_dir     = fullfile(bids_dir, 'derivatives', 'modelFit', 'figure', 'sfn_fixw');
    if ~exist(fullfile(fig_dir), 'dir'), mkdir(fullfile(fig_dir)); end
else
    git_dir = fileparts(pwd);
    data_dir = fullfile(git_dir, 'tactileECoGdata');
    fig_dir = fullfile(pwd, 'figures');
end

%% Load all data

epoch_t = [-0.4, 1.6];
models_to_plot = {'DN_fixw','LINEAR'};
str = 'umcudrouwen'; %umcudrouwen/ny726/group_average

modelColors = [218, 62, 82; 45, 125, 210]./255; % red and blue

%% Evaluate cross-validated R2

includeDerivedParams = false;
CV_R2 = nan(1, numel(models_to_plot));
cv_results = cell(1, numel(models_to_plot));
for ii = 1:numel(models_to_plot)
    model = models_to_plot{ii};
    cv_D = load(fullfile(data_dir, sprintf('%s_xvalmode1_electrodeaverages_%s.mat', model, str)));
    cv_results{ii} = tt_evaluateModelFit(cv_D,includeDerivedParams);
    CV_R2(ii) = cv_results{ii}.R2.concat_all;
end

%% Plot timecourses

tcourseFig = figure('Color', [1 1 1], 'Position', [0, 300, 1400, 400]);
set(tcourseFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[1000 400])
T1 = tiledlayout(2, 1,'TileIndexing','rowmajor');

% Load first model for metadata and share data for all plots
Dbase = load(fullfile(data_dir, sprintf('%s_xvalmode0_electrodeaverages_%s.mat', models_to_plot{1}, str)));
data = Dbase.data;
stim_info = Dbase.stim_info;
stim_ts = Dbase.stim;
t = Dbase.t;
channels = Dbase.channels;

is_onepulse = contains(stim_info.name, 'ONE-PULSE');
is_twopulse = contains(stim_info.name, 'ONE-PULSE-4') | contains(stim_info.name, 'TWO-PULSE');

onePulseIndx = find(is_onepulse);
twoPulseIndx = find(is_twopulse);
allPulsesIndx = [onePulseIndx; twoPulseIndx];
n_one = numel(onePulseIndx);
n_two = numel(twoPulseIndx);

xDur = stim_info.duration(onePulseIndx)';
xISI = stim_info.ISI(twoPulseIndx)';

sum_data  = sum(data,1);

% Set y-bounds based on actual data
bounds = @(x) [floor(min(x(:))*10) ceil(max(x(:))*10)]/10;
ybounds = bounds(prctile(data, [1 99], 'all'));
ybounds(1) = ybounds(1)-0.5;
ybounds(2) = ybounds(2)+0.5;

% Load predictions from both models
preds = cell(1, numel(models_to_plot));
sum_pred = cell(1, numel(models_to_plot));
for dd = 1:numel(models_to_plot)
    D = load(fullfile(data_dir, sprintf('%s_xvalmode0_electrodeaverages_%s.mat', models_to_plot{dd}, str)));
    preds{dd} = D.pred;
    sum_pred{dd} = sum(D.pred,1);
end

% Only plot model prediction of DN model
models_to_plot = {'DN_fixw'};
modelNames = models_to_plot;

% --- Panel for single pulse conditions, with one empty subplot at the beginning ---
nSingle = numel(onePulseIndx);
tpanel = tiledlayout(T1,1,nSingle+1,'TileIndexing','columnmajor');
tpanel.Layout.Tile = 1;
tpanel.XLabel.String = 'Time (s)';
tpanel.Title.String = sprintf('Cross-validated r^2: %s = %.3f, sub-%s', ...
    modelNames{1}, CV_R2(1),str);

% Add empty subplot as the first tile
nexttile(tpanel,1)
axis off

for ii = 1:nSingle
    nexttile(tpanel,ii+1)
    hold on
    % Stimulus and data once
    plot(t, stim_ts(:, onePulseIndx(ii)) * ybounds(end), 'Color', [.5 .5 .5], 'HandleVisibility', 'off')
    plot(t([1 end]), [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    plot(t, data(:,onePulseIndx(ii)), 'k-', 'LineWidth', 1, 'DisplayName', 'Data')
    % Loop for model predictions
    for dd = 1:numel(models_to_plot)
        plot(t, preds{dd}(:,onePulseIndx(ii)), 'Color', modelColors(dd,:), 'LineWidth', 1, ...
            'DisplayName', modelNames{dd})
    end
    set(gca,'TickDir', 'out', 'FontSize', 14, 'XColor', 'k', 'YColor', 'k', ...
        'LineWidth', 1, 'TickLength', [0.05 0.05])
    title(sprintf('Dur %.2fs', xDur(ii)), 'FontSize', 10)
    ylim(ybounds)
    box off

end

% --- Panel for paired pulse conditions ---
tpanel2 = tiledlayout(T1,1,n_two,'TileIndexing','columnmajor');
tpanel2.Layout.Tile = 2;

for ii = 1:n_two
    nexttile(tpanel2,ii)
    hold on
    plot(t, stim_ts(:, twoPulseIndx(ii)) * ybounds(end), 'Color', [.5 .5 .5], 'HandleVisibility', 'off')
    plot(t([1 end]), [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    plot(t, data(:,twoPulseIndx(ii)), 'k-', 'LineWidth', 1, 'DisplayName', 'Data')
    % Loop for model predictions
    for dd = 1:numel(models_to_plot)
        plot(t, preds{dd}(:,twoPulseIndx(ii)), 'Color', modelColors(dd,:), 'LineWidth', 1, ...
            'DisplayName', modelNames{dd})
    end
    set(gca,'TickDir', 'out', 'FontSize', 14, 'XColor', 'k', 'YColor', 'k', ...
        'LineWidth', 1, 'TickLength', [0.03 0.03])
    title(sprintf('ISI %.2fs', xISI(ii)), 'FontSize', 10)
    ylim(ybounds)
    box off

end

% Save figure only once with both models included
saveas(tcourseFig, fullfile(fig_dir, sprintf('%s_timecourse', str)), 'pdf');

%% Plot summed responses

models_to_plot = {'DN_fixw','LINEAR'};

summaryFig      = figure('Color', [1 1 1], 'Position', [30 300 800 260]);
set(summaryFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[550 260])
figure(summaryFig)
T = tiledlayout(1, 3,'TileIndexing','rowmajor');

% Load individual electrode data to get confidence interval
tactile_indiv = load(fullfile(data_dir,sprintf('DN_xvalmode0_individualelecs_%s.mat', str)), ...
    'data', 'pred', 'stim_info', 'stim', 't', 'channels', 'params');

cResponses68 = prctile(squeeze(sum(tactile_indiv.data, 1)), [15.87 84.13], 2);
cResponses95 = prctile(squeeze(sum(tactile_indiv.data, 1)), [2.5, 97.5], 2);

num_one_pulse = n_one;
num_two_pulse = n_two;

xvalues = xDur; % [0.05, 0.1, 0.2, 0.4, 0.8, 1.2];

% One pulse conditions
t1 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
t1.Layout.Tile = 1;

ax2 = nexttile(t1,[1 2]);
hold on
plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xvalues .* [1; 1], cResponses95(onePulseIndx,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xvalues .* [1; 1], cResponses68(onePulseIndx,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xvalues,  sum_data(onePulseIndx), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

% Paired pulse conditions
t2 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
t2.Layout.Tile = 2;
ax3 = nexttile(t2,[1 1]);
hold on
plot([-0.1 0.1], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xISI(1).* [1; 1], cResponses95(twoPulseIndx(1),:), 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(1).* [1; 1], cResponses68(twoPulseIndx(1),:), 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(1),  sum_data(twoPulseIndx(1)), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

ax4 = nexttile(t2,[1 2]);
hold on
plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xISI(2:end) .* [1; 1], cResponses95(twoPulseIndx(2:end),:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(2:end) .* [1; 1], cResponses68(twoPulseIndx(2:end),:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(2:end),  sum_data(twoPulseIndx(2:end)), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

% Linear xaxis
t3 = tiledlayout(T,2,1,'TileIndexing','columnmajor');
t3.Layout.Tile = 3;
ax5 = nexttile(t3,[1 1]);
hold on
plot([-0.05 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xDur .* [1; 1], cResponses95(onePulseIndx,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xDur .* [1; 1], cResponses68(onePulseIndx,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xDur,  sum_data(onePulseIndx), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

ax6 = nexttile(t3,[1 1]);
hold on
plot([-0.05 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xISI .* [1; 1], cResponses95(twoPulseIndx,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI.* [1; 1], cResponses68(twoPulseIndx,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI, sum_data(twoPulseIndx), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')


for dd = 1:numel(models_to_plot)
    plot(ax2, xDur, sum_pred{dd}(onePulseIndx), 'Color', modelColors(dd,:), 'LineWidth', 1.5)
    plot(ax3, xISI(1), sum_pred{dd}(twoPulseIndx(1)), 'Color', modelColors(dd,:), 'LineWidth', 1.5)
    plot(ax4, xISI(2:end), sum_pred{dd}(twoPulseIndx(2:end)), 'Color', modelColors(dd,:), 'LineWidth', 1.5)
    plot(ax5, xDur, sum_pred{dd}(onePulseIndx), 'Color', modelColors(dd,:), 'LineWidth', 1.5 )
    plot(ax6, xISI, sum_pred{dd}(twoPulseIndx), 'Color', modelColors(dd,:), 'LineWidth', 1.5 )
end

ybounds = get(ax2, 'YLim');

ax2.Box = 'off';
t1.Title.String = 'Single-pulse conditions';
t1.XLabel.String = 'Stimulus duration (s)';
t1.YLabel.String = {'Summed broadband time-series';'(X-fold)'};


if strcmp(str, 'ny726')
    tickspace = 200;
elseif strcmp(str, 'umcudrouwen')
    tickspace = 1000;
end

set(ax2, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
    'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds(1):tickspace:ybounds(2), 'YLim', ybounds, 'XLim', [0.035 1.3]);
ax2.Box = 'off';
% ax2.YAxis.Visible = 'off';

set(ax3, 'TickDir', 'out', 'XTick', 0, 'XTickLabel', 0, ...
    'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds(1):tickspace:ybounds(2), 'YLim', ybounds, 'XLim', [-0.005 0.01]);
ax3.Box = 'off';
t2.Title.String = 'Paired-pulse conditions';
t2.XLabel.String = 'Interstimulus interval (s)';
t2.YLabel.String = {'Summed broadband time-series';'(X-fold)'};

set(ax4, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
    'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds(1):tickspace:ybounds(2), 'YLim', ybounds, 'XLim', [0.04 1.3]);
ax4.Box = 'off';
ax4.YAxis.Visible = 'off';

set(ax5, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
    'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds(1):tickspace:ybounds(2), 'YLim', ybounds, 'XLim', [-0.1 1.3]);
ax5.Box = 'off';
ax5.Title.String = 'Single-pulse conditions';
ax5.XLabel.String = 'Duration (s)';
t3.YLabel.String = {'Summed broadband time-series';'(X-fold)'};

set(ax6, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
    'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds(1):tickspace:ybounds(2), 'YLim', ybounds, 'XLim', [-0.1 1.3]);
ax6.Box = 'off';
ax6.Title.String = 'Paired-pulse conditions';
ax6.XLabel.String = 'ISI (s)';

saveas(summaryFig, fullfile(fig_dir, sprintf('%s_summed_responses', str)), 'pdf');

%% Load fMRI data, see if we can scale model prediction to fMRI data

fmri_data_dir = '/Volumes/server/Projects/TemporalTactileCounting/Data/modelOutput';

% Select one of the model outputs since we only need the summed data
fmri_results(1) = load(fullfile(fmri_data_dir, 'HRF', 'sub-group_localizerROI-S1_model-HRF_crossval-noCross_optimizer-bads_modelOutput.mat'));
fmri_results(2) = load(fullfile(fmri_data_dir, 'HRF', 'sub-group_localizerROI-S1_model-HRF_crossval-noCross_optimizer-bads_btsmodelOutput.mat'));

% combine one-pulse conditions
is_fmri_onepulse = contains(fmri_results(1).condOrder, 'BLANK') | contains(fmri_results(1).condOrder, 'ONE');
is_fmri_twopulse = contains(fmri_results(1).condOrder, 'ONE-PULSE-4') | contains(fmri_results(1).condOrder, 'TWO');

onePulseIndx    = find(is_fmri_onepulse);
twoPulseIndx    = find(is_fmri_twopulse);
allPulsesIndx   = [onePulseIndx; twoPulseIndx];
n_fmri_one = numel(onePulseIndx);
n_fmri_two = numel(twoPulseIndx);

% Corresponding ecog conditions
predOnePulseIndx = 1:6; % There's no blank prediction from ecog
predTwoPulseIndx = [4,7:12];

fmri_summed_data = zeros(numel(allPulsesIndx), 1);
[fmri_ci_95, fmri_ci_68] = deal(zeros(numel(allPulsesIndx), 2));

for ii = 1:numel(allPulsesIndx)
    fmri_summed_data(ii) = sum(fmri_results(1).y_data(:,allPulsesIndx(ii)));
    fmri_btst = sum(squeeze(fmri_results(2).y_data(:,allPulsesIndx(ii),:)),1);
    fmri_ci_95(ii,:) = prctile(fmri_btst, [2.5, 97.5]);
    fmri_ci_68(ii,:) = prctile(fmri_btst, [15.87, 84.13]);
end

% Predict fMRI summed data by scaling model predictions of ecog data
scaler = mean(fmri_summed_data) / mean(sum_data);
fmri_pred = scaler * sum_pred{1};

%% make plot of fmri

fmriSummaryFig = figure('Color', [1 1 1], 'Position', [30 300 800 260]);
set(fmriSummaryFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[550 260])
figure(fmriSummaryFig)
T = tiledlayout(1, 3,'TileIndexing','rowmajor');

ybounds_fmri = [-1 5];

%-- one pulse
t1 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
t1.Layout.Tile = 1;
ax1 = nexttile(t1,[1 1]);
hold on,
plot([-0.1 0.1], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(0 .* [1; 1],fmri_ci_95(1,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(0 .* [1; 1], fmri_ci_68(1,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(0,  fmri_summed_data(1), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

ax2 = nexttile(t1,[1 2]);
hold on
plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xDur .* [1; 1], fmri_ci_95(2:n_fmri_one,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xDur.*[1; 1], fmri_ci_68(2:n_fmri_one,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xDur, fmri_summed_data(2:n_fmri_one), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

%-- paired pulse
t2 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
t2.Layout.Tile = 2;
ax3 = nexttile(t2,[1 1]);
hold on
plot([-0.1 0.1], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xISI(1) .* [1; 1], fmri_ci_95((1)+n_fmri_one,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(1) .* [1; 1], fmri_ci_68((1)+n_fmri_one,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(1),  fmri_summed_data((1)+n_fmri_one), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

ax4 = nexttile(t2,[1 2]);
hold on
plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xISI(2:end) .* [1; 1], fmri_ci_95((2+n_fmri_one):end,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(2:end).* [1; 1], fmri_ci_68((2+n_fmri_one):end,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(2:end),  fmri_summed_data((2+n_fmri_one):end), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

%-- linear xaxis
t3 = tiledlayout(T,2,1,'TileIndexing','columnmajor');
t3.Layout.Tile = 3;
ax5 = nexttile(t3, [1 1]);
hold on
plot([-0.05 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot([0 xDur].* [1; 1], fmri_ci_95(1:n_fmri_one,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot([0 xDur].* [1; 1], fmri_ci_68(1:n_fmri_one,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot([0 xDur],  fmri_summed_data(1:n_fmri_one), '.k', 'MarkerSize', 10,  'HandleVisibility', 'off')

ax6 = nexttile(t3,[1 1]);
hold on
plot([-0.05 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xISI .* [1; 1], fmri_ci_95((1:n_fmri_one)+n_fmri_one,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI .* [1; 1], fmri_ci_68((1:n_fmri_one)+n_fmri_one,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI,  fmri_summed_data((1:n_fmri_one)+n_fmri_one), '.k', 'MarkerSize', 10,  'HandleVisibility', 'off')

% Model predictions
plot(ax2, xDur, fmri_pred(predOnePulseIndx), 'Color', modelColors(1,:), 'LineWidth', 1.5)
scatter(ax3, xISI(1)+0.1, fmri_pred(predTwoPulseIndx(1)), 30,'filled','MarkerFaceColor', modelColors(1,:),'MarkerEdgeColor','none')
plot(ax4, xISI(2:end), fmri_pred(predTwoPulseIndx(2:end)), 'Color', modelColors(1,:), 'LineWidth', 1.5)
plot(ax5, xDur, fmri_pred(predOnePulseIndx), 'Color', modelColors(1,:), 'LineWidth', 1.5 )
plot(ax6, xISI, fmri_pred(predTwoPulseIndx), 'Color', modelColors(1,:), 'LineWidth', 1.5 )

% Linear prediction
lin_pred = sum(fmri_results(1).y_est,1);
plot(ax2, xDur, lin_pred(onePulseIndx(2:end)), 'Color', modelColors(2,:), 'LineWidth', 1.5)
plot(ax3, xISI(1)+0.2, lin_pred(twoPulseIndx(1)), 'Color', modelColors(2,:), 'LineWidth', 1.5)
plot(ax4, xISI(2:end), lin_pred(twoPulseIndx(2:end)), 'Color', modelColors(2,:), 'LineWidth', 1.5)
plot(ax5, xDur, lin_pred(onePulseIndx(2:end)), 'Color', modelColors(2,:), 'LineWidth', 1.5 )
plot(ax6, xISI, lin_pred(twoPulseIndx), 'Color', modelColors(2,:), 'LineWidth', 1.5 )

tickspace = 1;
ax2.Box = 'off';
t1.Title.String = 'Single-pulse conditions';
t1.XLabel.String = 'Stimulus duration (s)';
t1.YLabel.String = {'Summed BOLD time series';'(%SC)'};

set(ax1, 'TickDir', 'out', 'XTick', 0, 'XTickLabel', 0, ...
    'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):1:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [-0.005 0.01]);
ax2.Box = 'off';
t1.Title.String = 'Single pulse conditions';
t1.XLabel.String = 'Stimulus duration (s)';
t1.YLabel.String = {'Summed BOLD time series';'(%SC)'};

set(ax2, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
    'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):tickspace:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [0.035 1.3]);
ax2.Box = 'off';
ax2.YAxis.Visible = 'off';

set(ax3, 'TickDir', 'out', 'XTick', 0, 'XTickLabel', 0, ...
    'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):tickspace:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [-0.005 0.01]);
ax3.Box = 'off';
t2.Title.String = 'Paired-pulse conditions';
t2.XLabel.String = 'Interstimulus interval (s)';
t2.YLabel.String = {'Summed BOLD time series';'(%SC)'};

set(ax4, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
    'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):tickspace:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [0.04 1.3]);
ax4.Box = 'off';
ax4.YAxis.Visible = 'off';

set(ax5, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
    'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):tickspace:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [-0.1 1.3]);
ax5.Box = 'off';
ax5.Title.String = 'Single-pulse conditions';
ax5.XLabel.String = 'Duration (s)';
t3.YLabel.String = {'Summed BOLD time series';'(%SC)'};

set(ax6, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
    'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):tickspace:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [-0.1 1.3]);
ax6.Box = 'off';
ax6.Title.String = 'Paired-pulse conditions';
ax6.XLabel.String = 'ISI (s)';

saveas(fmriSummaryFig, fullfile(fig_dir, sprintf('%s_fmri_summed_responses', str)), 'pdf');