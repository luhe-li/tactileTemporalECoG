% SFN ecog figures

% Plot DN and LIN model fits of group-averaged, electrode-averaged tactile data
% Model fits of average-electrode are under bidsRootPath/derivatives/modelFit/results

close all

%% Setup paths and parameters
server_available = true;
if server_available
    bids_dir = tt_bidsRootPath;
    data_dir = fullfile(bids_dir, 'derivatives', 'modelFit', 'results');
    fig_dir = fullfile(bids_dir, 'derivatives', 'modelFit', 'figure', 'sfn_fixw');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
else
    git_dir = fileparts(pwd);
    data_dir = fullfile(git_dir, 'tactileECoGdata');
    fig_dir = fullfile(pwd, 'figures');
end

epoch_t = [-0.4, 1.6];
models_to_plot = {'DN_fixw','LINEAR'};
modelLabels = {'DN', 'LIN'};
nModels = numel(models_to_plot);
str = 'umcudrouwen'; %umcudrouwen/ny726/group_average

modelColors = [218, 62, 82; 45, 125, 210]./255; % red and blue

%% Evaluate cross-validated R2 for each model

includeDerivedParams = false;
CV_R2 = nan(1, numel(models_to_plot));
cv_results = cell(1, numel(models_to_plot));
for ii = 1:numel(models_to_plot)
    model = models_to_plot{ii};
    cv_D = load(fullfile(data_dir, sprintf('%s_xvalmode1_electrodeaverages_%s.mat', model, str)));
    cv_results{ii} = tt_evaluateModelFit(cv_D,includeDerivedParams);
    CV_R2(ii) = cv_results{ii}.R2.concat_all;
end


%% Load first model for metadata and common variables

D = load(fullfile(data_dir, sprintf('%s_xvalmode0_electrodeaverages_%s.mat', models_to_plot{1}, str)));
data = D.data;
stim_info = D.stim_info;
stim_ts = D.stim;
t = D.t;
channels = D.channels;

is_onepulse = contains(stim_info.name, 'ONE-PULSE');
is_twopulse = contains(stim_info.name, 'ONE-PULSE-4') | contains(stim_info.name, 'TWO-PULSE');
idx_onepulse = find(is_onepulse);
idx_twopulse = find(is_twopulse); 

onePulseIndx = find(is_onepulse);
twoPulseIndx = find(is_twopulse);
allPulsesIndx = [onePulseIndx; twoPulseIndx];
n_one = numel(onePulseIndx);
n_two = numel(twoPulseIndx);

xDur = stim_info.duration(onePulseIndx)';
xISI = stim_info.ISI(twoPulseIndx)';

sum_data = sum(data,1);

% For both model predictions
preds = cell(1, numel(models_to_plot));
sum_pred = cell(1, numel(models_to_plot));
for ii = 1:numel(models_to_plot)
    D = load(fullfile(data_dir, sprintf('%s_xvalmode0_electrodeaverages_%s.mat', models_to_plot{ii}, str)));
    preds{ii} = D.pred;
    sum_pred{ii} = sum(D.pred,1);
end

%% Plot timecourses for each condition

tcourseFig = figure('Color', [1 1 1], 'Position', [0, 0, 1000, 300]);
set(tcourseFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[1000 300])
T1 = tiledlayout(2, 1,'TileIndexing','rowmajor');

bounds = @(x) [floor(min(x(:))*10) ceil(max(x(:))*10)]/10;
ybounds = bounds(prctile(data, [1 99], 'all'));
ybounds = [ybounds(1)-0.5, ybounds(2)+0.5];

nSingle = numel(onePulseIndx);
tpanel = tiledlayout(T1,1,nSingle+1,'TileIndexing','columnmajor');
tpanel.Layout.Tile = 1;
tpanel.XLabel.String = 'Time (s)';
tpanel.Title.String = sprintf('Cross-validated r^2: %s = %.3f, sub-%s', ...
    modelLabels{1}, CV_R2(1),str);

nexttile(tpanel,1); axis off;

for i = 1:nSingle
    nexttile(tpanel,i+1)
    hold on
    plot(t, stim_ts(:, onePulseIndx(i)) * ybounds(end), 'Color', [.5 .5 .5], 'HandleVisibility', 'off','LineWidth', 1.5)
    plot(t([1 end]), [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    plot(t, data(:,onePulseIndx(i)), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Data')
    for j = numel(models_to_plot):-1:1
        plot(t, preds{j}(:,onePulseIndx(i)), 'Color', modelColors(j,:), 'LineWidth', 2, 'DisplayName', modelLabels{j})
    end
    set(gca,'TickDir', 'out', 'FontSize', 14, 'XColor', 'k', 'YColor', 'k', ...
        'LineWidth', 1, 'TickLength', [0.05 0.05])
    title(sprintf('Dur %.2fs', xDur(i)), 'FontSize', 10)
    ylim(ybounds)
    if i ~= 1
        yticks([])
    end
    box off
    
end

tpanel2 = tiledlayout(T1,1,n_two,'TileIndexing','columnmajor');
tpanel2.Layout.Tile = 2;
for i = 1:n_two
    nexttile(tpanel2,i)
    hold on
    plot(t, stim_ts(:, twoPulseIndx(i)) * ybounds(end), 'Color', [.5 .5 .5], 'HandleVisibility', 'off', 'LineWidth',1.5)
    plot(t([1 end]), [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    plot(t, data(:,twoPulseIndx(i)), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Data')
    for j = numel(models_to_plot):-1:1
        plot(t, preds{j}(:,twoPulseIndx(i)), 'Color', modelColors(j,:), 'LineWidth', 2, 'DisplayName', modelLabels{j})
    end
    set(gca,'TickDir', 'out', 'FontSize', 14, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, 'TickLength', [0.03 0.03])
    title(sprintf('ISI %.2fs', xISI(i)), 'FontSize', 10)
    if i ~= 1
        yticks([])
    end
    ylim(ybounds)
    box off
end
saveas(tcourseFig, fullfile(fig_dir, sprintf('%s_%s_timecourse_%s', str, models_to_plot{1})), 'pdf');

%% Plot summed responses

summaryFig = figure('Color', [1 1 1], 'Position', [0 0 1000 250]);
set(summaryFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[800 250])
T = tiledlayout(1, 3,'TileIndexing','rowmajor');

tactile_indiv = load(fullfile(data_dir,sprintf('DN_fixw_xvalmode0_individualelecs_%s.mat', str)), ...
    'data', 'pred', 'stim_info', 'stim', 't', 'channels', 'params');
cResponses68 = prctile(squeeze(sum(tactile_indiv.data, 1)), [15.87 84.13], 2);
cResponses95 = prctile(squeeze(sum(tactile_indiv.data, 1)), [2.5, 97.5], 2);

x = stim_info.duration(onePulseIndx)';
if numel(x) == 0; x = [0.05, 0.1, 0.2, 0.4, 0.8, 1.2]; end 

% One pulse: main axis (ax2)
t1 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
t1.Layout.Tile = 1;
ax2 = nexttile(t1,[1 2]);
hold on
plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xDur .* [1; 1], cResponses95(1:n_one,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xDur .* [1; 1], cResponses68(1:n_one,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xDur,  sum_data(onePulseIndx), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

% Paired pulse: first (0) and rest (ax3/ax4)
t2 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
t2.Layout.Tile = 2;
ax3 = nexttile(t2,[1 1]);
hold on
plot([-0.1 0.4], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xISI(1).* [1; 1], cResponses95(twoPulseIndx(1),:), 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(1).* [1; 1], cResponses68(twoPulseIndx(1),:), 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(1),  sum_data(twoPulseIndx(1)), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

ax4 = nexttile(t2,[1 2]);
hold on
plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(xISI(2:end) .* [1; 1], cResponses95(twoPulseIndx(2:end),:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(2:end) .* [1; 1], cResponses68(twoPulseIndx(2:end),:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(xISI(2:end),  sum_data(twoPulseIndx(2:end)), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

% Linear xaxis for one- and two-pulse (ax5/ax6)
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

% Model predictions of one-pulse
nStimNew = 100;
[stim2, stim_info2] = tt_simulateNewStimuli(t,nStimNew);
stim_idx2 = find(contains(stim_info2.name, 'ONEPULSE'));
one_pulse_x2 = stim_info2.duration(stim_idx2)';
one_pulse_pred2 = nan([nModels size(stim2(:,stim_idx2))]);
srate = D.channels.sampling_frequency(1);
for kk = 1:nModels
    DD = load(fullfile(data_dir, sprintf('%s_xvalmode0_electrodeaverages_%s.mat', models_to_plot{kk}, str)));
    prm = DD.params;

    % pred2: model x timepoint x trial/condition
    [~, one_pulse_pred2(kk,:,:)] = DD.objFunction(prm, [], stim2(:,stim_idx2), srate);

    % Plot
    plot(ax2, one_pulse_x2, sum(squeeze(one_pulse_pred2(kk,:,:)),1), 'Color', modelColors(kk,:), 'LineWidth', 2)
    plot(ax5, one_pulse_x2, sum(squeeze(one_pulse_pred2(kk,:,:)),1), 'Color', modelColors(kk,:), 'LineWidth', 2)
end

% Model predictions of two-pulse
nStimNew = 100;
[stim2, stim_info2] = tt_simulateNewStimuli(t,nStimNew);
stim_idx2 = find(contains(stim_info2.name, 'TWOPULSE'));
stim2 = stim2(:,stim_idx2);
stim_info2 = stim_info2(stim_idx2,:);
two_pulse_x2 = stim_info2.ISI';
two_pulse_pred2 = nan([nModels size(stim2)]);

% Overplot model predictions
for kk = 1:numel(models_to_plot)
    DD = load(fullfile(data_dir, sprintf('%s_xvalmode0_electrodeaverages_%s.mat', models_to_plot{kk}, str)));
    prm = DD.params;

    % pred2: model x timepoint x trial/condition
    [~, two_pulse_pred2(kk,:,:)] = DD.objFunction(prm, [], stim2, srate);

    x_ax2 = two_pulse_x2<0.03;
    x_ax3 = two_pulse_x2>=0.03;
    plot(ax3, two_pulse_x2(x_ax2), sum(squeeze(two_pulse_pred2(kk,:,x_ax2)),1), 'Color', modelColors(kk,:), 'LineWidth',2)
    plot(ax4, two_pulse_x2(x_ax3), sum(squeeze(two_pulse_pred2(kk,:,x_ax3)),1), 'Color', modelColors(kk,:), 'LineWidth',2)
    plot(ax6, two_pulse_x2, sum(squeeze(two_pulse_pred2(kk,:,:)),1), 'Color', modelColors(kk,:), 'LineWidth', 2)
end

ybounds = get(ax2, 'YLim');
ax2.Box = 'off';
t1.Title.String = 'Single-pulse conditions';
t1.XLabel.String = 'Stimulus duration (s)';
t1.YLabel.String = {'Summed broadband time-series';'(X-fold)'};

tickspace = 500;
if strcmp(str, 'ny726'), tickspace = 200; end

set(ax2, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
    'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds(1):tickspace:ybounds(2), 'YLim', ybounds, 'XLim', [0.035 1.3]);
set(ax3, 'TickDir', 'out', 'XTick', 0, 'XTickLabel', 0, ...
    'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds(1):tickspace:ybounds(2), 'YLim', ybounds, 'XLim', [-0.005 0.01]);
set(ax4, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
    'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds(1):tickspace:ybounds(2), 'YLim', ybounds, 'XLim', [0.04 1.3]); ax4.YAxis.Visible = 'off';
set(ax5, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
    'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds(1):tickspace:ybounds(2), 'YLim', ybounds, 'XLim', [-0.1 1.3]);
set(ax6, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
    'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds(1):tickspace:ybounds(2), 'YLim', ybounds, 'XLim', [-0.1 1.3]);

ax5.Title.String = 'Single-pulse conditions'; ax5.XLabel.String = 'Duration (s)';
ax6.Title.String = 'Paired-pulse conditions'; ax6.XLabel.String = 'ISI (s)';
t3.YLabel.String = {'Summed broadband time-series';'(X-fold)'};

saveas(summaryFig, fullfile(fig_dir, sprintf('%s_summed_responses', str)), 'pdf');

%% Load fMRI data

clearvars -except xDur xISI stim_info data_dir modelColors fig_dir models_to_plot one_pulse_pred2 two_pulse_pred2 one_pulse_x2 two_pulse_x2

% ----- load fMRI data -----
fmri_data_dir = '/Volumes/server/Projects/TemporalTactileCounting/Data/modelOutput';
fmri_results(1) = load(fullfile(fmri_data_dir, 'HRF', 'sub-group_localizerROI-S1_model-HRF_crossval-noCross_optimizer-bads_modelOutput.mat'));
fmri_results(2) = load(fullfile(fmri_data_dir, 'HRF', 'sub-group_localizerROI-S1_model-HRF_crossval-noCross_optimizer-bads_btsmodelOutput.mat'));

% Indices for fMRI and mapping to ECoG model
f_one = find(contains(fmri_results(1).condOrder, 'ONE'));
f_blank = find(contains(fmri_results(1).condOrder, 'BLANK'));
f_onepulse = [f_blank; f_one]';
f_twopulse = [find(contains(fmri_results(1).condOrder, 'ONE-PULSE-4')); ...
    find(contains(fmri_results(1).condOrder, 'TWO'))]';
f_allpulse = [f_onepulse, f_twopulse];

predOnePulseIndx = 1:numel(xDur);
predTwoPulseIndx = [find(contains(stim_info.name, 'ONE-PULSE-4'),1,'first'); find(contains(stim_info.name, 'TWO-PULSE'))]';

fmri_summed_data = nan(numel(f_allpulse),1);
fmri_ci_95 = nan(numel(f_allpulse),2);
fmri_ci_68 = nan(numel(f_allpulse),2);

for i = 1:numel(f_allpulse)
    fmri_summed_data(i) = sum(fmri_results(1).y_data(:,f_allpulse(i)));

    % sum the time points first, take confidence interval next
    % fmri_btst = sum(squeeze(fmri_results(2).y_data(:,f_allpulse(i),:)), 1);
    % fmri_ci_95(i,:) = prctile(fmri_btst, [2.5, 97.5]);
    % fmri_ci_68(i,:) = prctile(fmri_btst, [15.87, 84.13]);

    % take confidence interval for each time point across bootstrapped
    % trials, and then sum
    fmri_ci_95(i,:) = sum(prctile(squeeze(fmri_results(2).y_data(:,f_allpulse(i),:)), [2.5, 97.5],2),1);
    fmri_ci_68(i,:) = sum(prctile(squeeze(fmri_results(2).y_data(:,f_allpulse(i),:)), [15.87, 84.13],2),1);
end

% ----- predict fmri from group averaged ECoG data -----
subjects = {'umcudrouwen','ny726'};
num_subs = numel(subjects);

% Load and extract DN model predictions for each subject
for ss = 1%:num_subs
    sub_str = subjects{ss};

    for mm = 1:numel(models_to_plot)
        dn_path = fullfile(data_dir, sprintf('%s_xvalmode0_electrodeaverages_%s.mat', models_to_plot{mm}, sub_str));
        
        if exist(dn_path,'file')
            Dsub = load(dn_path);
            sum_preds{mm}(ss,:) = squeeze(sum(Dsub.pred,1)); % summed prediction per condition

            % Data are the same for all models
            ieeg_data(ss,:) = sum(Dsub.data,1);
        else
            error('File does not exist: %s', dn_path)
        end
    end
end

% Average across subjects
mean_pred_dn = mean(sum_preds{1},1);
mean_pred_lin = mean(sum_preds{2},1);

% Scale up by the mean of the fMRI data
scaler = mean(fmri_summed_data(:))/mean(ieeg_data(:));
% model_pred2: model x timepoint x trial/condition
% Sum across timeponts
fmri_one_pulse_pred =  sum(one_pulse_pred2,2) * scaler;
fmri_two_pulse_pred =  sum(two_pulse_pred2,2) * scaler;

% Duration and ISI manipulation is the same
x = xISI;

%% Plot: fMRI + scaled summed broadband model predictions

fmriSummaryFig = figure('Color', [1 1 1], 'Position', [0 0 1000 250]);
set(fmriSummaryFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[1000 250])
T = tiledlayout(1, 3,'TileIndexing','rowmajor');
ybounds_fmri = [-1 5];
breakpoint = 0.02;

%-- one pulse
t1 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
t1.Layout.Tile = 1;
ax1 = nexttile(t1,[1 1]);
hold on
plot([-0.1 breakpoint], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(x(1) .* [1; 1],fmri_ci_95(1,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x(1) .* [1; 1], fmri_ci_68(1,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x(1),  fmri_summed_data(1), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

ax2 = nexttile(t1,[1 2]);
hold on
plot([breakpoint 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(x(2:end) .* [1; 1], fmri_ci_95(2:numel(f_onepulse),:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x(2:end) .* [1; 1], fmri_ci_68(2:numel(f_onepulse),:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x(2:end), fmri_summed_data(2:numel(f_onepulse)), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

%-- paired pulse
t2 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
t2.Layout.Tile = 2;
ax3 = nexttile(t2,[1 1]);
hold on
plot([-0.1 breakpoint], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(x(1) .* [1; 1], fmri_ci_95((1)+numel(f_onepulse),:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x(1) .* [1; 1], fmri_ci_68((1)+numel(f_onepulse),:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x(1),  fmri_summed_data((1)+numel(f_onepulse)), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

ax4 = nexttile(t2,[1 2]);
hold on
plot([breakpoint 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(x(2:end) .* [1; 1], fmri_ci_95((2+numel(f_onepulse)):end,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x(2:end) .* [1; 1], fmri_ci_68((2+numel(f_onepulse)):end,:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x(2:end),  fmri_summed_data((2+numel(f_onepulse)):end), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')

%-- linear xaxis
t3 = tiledlayout(T,2,1,'TileIndexing','columnmajor');
t3.Layout.Tile = 3;
ax5 = nexttile(t3,[1 1]); hold on
plot([-0.05 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(x.* [1; 1], fmri_ci_95(1:numel(f_onepulse),:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x.* [1; 1], fmri_ci_68(1:numel(f_onepulse),:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x,  fmri_summed_data(1:numel(f_onepulse)), '.k', 'MarkerSize', 10,  'HandleVisibility', 'off')

ax6 = nexttile(t3,[1 1]); hold on
plot([-0.05 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
plot(x .* [1; 1], fmri_ci_95((1:numel(f_onepulse))+numel(f_onepulse),:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x .* [1; 1], fmri_ci_68((1:numel(f_onepulse))+numel(f_onepulse),:)', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
plot(x,  fmri_summed_data((1:numel(f_onepulse))+numel(f_onepulse)), '.k', 'MarkerSize', 10,  'HandleVisibility', 'off')

for kk = 1:2


    % one pulse
    x_ax2 = one_pulse_x2<breakpoint;
    x_ax3 = one_pulse_x2>=breakpoint;
    plot(ax1, one_pulse_x2(x_ax2), fmri_one_pulse_pred(kk,x_ax2), 'Color', modelColors(kk,:), 'LineWidth', 2)
    plot(ax2, one_pulse_x2(x_ax3), fmri_one_pulse_pred(kk,x_ax3), 'Color', modelColors(kk,:), 'LineWidth', 2)
    plot(ax5, one_pulse_x2, fmri_one_pulse_pred(kk,:), 'Color', modelColors(kk,:), 'LineWidth', 2)

    % two pulse 
        x_ax2 = two_pulse_x2<breakpoint;
    x_ax3 = two_pulse_x2>=breakpoint;
    plot(ax3, two_pulse_x2(x_ax2), fmri_two_pulse_pred(kk,x_ax2), 'Color', modelColors(kk,:), 'LineWidth',2)
    plot(ax4, two_pulse_x2(x_ax3), fmri_two_pulse_pred(kk,x_ax3), 'Color', modelColors(kk,:), 'LineWidth',2)
    plot(ax6, two_pulse_x2, fmri_two_pulse_pred(kk,:), 'Color', modelColors(kk,:), 'LineWidth', 2)
end

tickspace = 1;
ax2.Box = 'off';
t1.Title.String = 'Single-pulse conditions';
t1.XLabel.String = 'Stimulus duration (s)';
t1.YLabel.String = {'Summed BOLD time series';'(%SC)'};

set(ax1, 'TickDir', 'out', 'XTick', 0, 'XTickLabel', 0, ...
    'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(tickspace):1:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [-0.005 0.01]);
set(ax2, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
    'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):tickspace:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [breakpoint 1.3]); ax2.YAxis.Visible = 'off';
set(ax3, 'TickDir', 'out', 'XTick', 0, 'XTickLabel', 0, ...
    'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):tickspace:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [-0.005 0.01]);
set(ax4, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
    'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):tickspace:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [breakpoint 1.3]); ax4.YAxis.Visible = 'off';
set(ax5, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
    'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):tickspace:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [-0.1 1.3]);
set(ax6, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
    'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
    'YTick', ybounds_fmri(1):tickspace:ybounds_fmri(2), 'YLim', ybounds_fmri, 'XLim', [-0.1 1.3]);
ax5.Title.String = 'Single-pulse conditions'; ax5.XLabel.String = 'Duration (s)';
ax6.Title.String = 'Paired-pulse conditions'; ax6.XLabel.String = 'ISI (s)';
t3.YLabel.String = {'Summed BOLD time series';'(%SC)'};

saveas(fmriSummaryFig, fullfile(fig_dir, sprintf('predicted_fmri_summed_responses')), 'pdf');
