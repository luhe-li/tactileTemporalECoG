function s6_createFigures_ECoG(saveFig, seed)
% Summary figure illustrating temporal dynamics in ECoG data
%
% 2025 - Ilona Bloem
%
% Dependencies: ECoG_utils toolbox

if ~exist('seed', 'var') || isempty(seed)
   seed         = rng('shuffle', 'twister'); % random seed based on current time
end

if ~exist('saveFig', 'var') || isempty(saveFig)
    saveFig     = false;
end

bounds          = @(x) [floor(min(x(:))*10) ceil(max(x(:))*10)]/10;

%-- define seed for reproducibility
rng(seed);

%-- Set up paths
dataStr         = 'ECoG';
fileName        = 'individualelecs';
[~, dataRootDir]= rootPath(false);
figRoot         = fullfile(dataRootDir, '..', 'Figures');
figDir          = fullfile(figRoot, dataStr);
if ~exist(fullfile(figDir), 'dir'), mkdir(fullfile(figDir)); end

%-- Find data
fnameList       = dir(fullfile(dataRootDir, dataStr, sprintf('DN*%s.mat', fileName)));
assert(length(fnameList) == 1, 'ECoG data files not found, verify paths')

%-- Load results
%modelColor      = [0.1540    0.5902    0.9218];
results         = load(fullfile(dataRootDir, dataStr, fnameList.name), ...
                    'data', 'pred', 'stim_info', 'stim', 't', 'channels');

x_data          = results.t;
stim_ts         = results.stim;
stim_info       = results.stim_info;
xaxisScale      = 'lin'; % 'log' or 'lin'

%-- bootstrap across electrodes
numBoot         = 500;
CIrange         = 32; % to compute 68% CI
[data, data_se] = bootstrWithinArea(results.data, @median, numBoot, CIrange);

% setup indexes
one_idx         = find(contains(stim_info.name, 'ONEPULSE'));
xDur            = stim_info.duration(one_idx); % in s
pair_idx        = find(contains(stim_info.name, {'ONEPULSE-4', 'TWOPULSE'}));
xISI            = stim_info.ISI(pair_idx); % 

%% - ECoG figure

% Panel A: ECoG time courses
ecogFig      = figure('Color', [1 1 1], 'Position', [30 300 700 400]);
set(ecogFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[700 400])
T1 = tiledlayout(2, 1,'TileIndexing','rowmajor');

% combine one-pulse conditions
onePulseIndx    = find(contains(results(1).stim_info.name, 'ONE'));

% combine two-pulse conditions  
twoPulseIndx    = cat(1, find(contains(results(1).stim_info.name, 'ONEPULSE-4')), ...
find(contains(results(1).stim_info.name, 'TWO')));

allPulsesIndx   = cat(1, onePulseIndx, twoPulseIndx);

ybounds         = bounds(prctile([data, data], [1 99], [1 3]));
if mod(ybounds(1), 0.2) > 0; ybounds(1) = ybounds(1)-0.2; end
if mod(ybounds(2), 0.2) > 0; ybounds(2) = ybounds(2)+0.5; end
Responses       = NaN(numel(allPulsesIndx), numel(results(1).t));
ciResponses     = NaN(numel(allPulsesIndx), numel(results(1).t), 2);

t = tiledlayout(T1,1,7,'TileIndexing','columnmajor');
t.Layout.Tile = 1;
t.XLabel.String = 'Time (s)';
t.Title.String = 'Single Pulse Conditions';
nexttile(t,[1 1])
ylim(ybounds)

%-- Extract data for all pulse conditions
for ii = 1:numel(allPulsesIndx)

    % time series data & predictions
    Responses(ii,:)         = data(:,allPulsesIndx(ii)); %results(1).data(:,allPulsesIndx(ii));

    % confidence interval data based on bootstrapped data
    ciResponses(ii,:,:)     = data_se(:,allPulsesIndx(ii),:);

    %-- Visualize time courses for all pulse conditions
    figure(ecogFig)

    if ii == 7
        t = tiledlayout(T1,1,7,'TileIndexing','columnmajor');
        t.Layout.Tile = 2;
        t.XLabel.String = 'Time (s)';
        t.Title.String = 'Paired Pulse Conditions';

    end

    nexttile(t,[1 1])
    hold on,
    plot(x_data, stim_ts(:, allPulsesIndx(ii)) * ybounds(end), 'Color', [.5 .5 .5], 'HandleVisibility', 'off')

    plot(x_data([1 end]), [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
    hCI = fill([x_data', fliplr(x_data')], [squeeze(ciResponses(ii,:,1)), fliplr(squeeze(ciResponses(ii,:,2)))], 'k', 'HandleVisibility', 'off');
    hCI.FaceAlpha = 0.2; hCI.EdgeColor = [1 1 1];
    plot(x_data, Responses(ii,:)', 'k-', 'LineWidth', 2, 'HandleVisibility', 'off')
    if ii <= numel(onePulseIndx)
        title(sprintf('Dur %.2fs', xDur(ii)), 'FontSize', 10)
    else
        title(sprintf('ISI %.2fs', xISI(ii-numel(onePulseIndx))), 'FontSize', 10)
    end
    set(gca,'TickDir', 'out', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1)
    ylim(ybounds)
    % set(gca, 'ytick', ybounds(1):0.4:ybounds(2), 'tickdir', 'out')
    box off

end

if saveFig > 0
    print(ecogFig, fullfile(figDir, sprintf('panelA_ECoGtimecourses')), '-dpdf')
end

%% -- panel B: temporal sub-additivity
subAddFig      = figure('Color', [1 1 1], 'Position', [30 300 350 400]);
set(subAddFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[350 400])
T1 = tiledlayout(4, 1,'TileIndexing','rowmajor');

% create linear system prediction
allResults              = struct;
allResults.condNames    = stim_info.name(one_idx);
allResults.stimDur      = xDur;
allResults.resp         = data(:, one_idx)'; 
allResults.x_data       = x_data;

createLinearSystemPredictions(allResults, subAddFig);

if saveFig > 0
    print(subAddFig, fullfile(figDir, sprintf('panelB_ECoGsubAdditivity')), '-dpdf')
end

%{
%-- panel B: One pulse stimuli to show compressive temporal summation
t = tiledlayout(T1,1,2,'TileIndexing','columnmajor');
t.Layout.Tile       = 1;

% Select two conditions to plot
conditionsOfInterest = {'ONEPULSE-1', 'ONEPULSE-2'};

% Select time window to plot
timepointsOfInterest = [-0.05 0.3];
x_idx                   = x_data>timepointsOfInterest(1) & x_data<=timepointsOfInterest(2);

% Look up corresponding indices in the data
tmp_idx             = find(contains(stim_info.name, conditionsOfInterest));
maxResp             = max(data(:, tmp_idx), [], 'all');

% create time shifted copy of the response to create a linear prediction
stimDur             = stim_info.duration(tmp_idx) * 1000;
padZeros            = 1000; % add a second of zeros as padding
tmpPred             = cat(1, data(:,tmp_idx(1)), zeros(padZeros,1));

% how many shifted copies of the shorter pulse to create
for numShifts = 1:(stimDur(2) / stimDur(1) -1)   
    tmpPred         = tmpPred + cat(1, zeros(numShifts * stimDur(1), 1), ...
                    data(:,tmp_idx(1)), zeros(padZeros - (numShifts * stimDur(1)), 1));
end

% remove padding 
linPred             = tmpPred(1:size(data,1)); 
ybounds             = bounds(linPred ./ maxResp);

% plot time courses + stimulus duration
for ii = 1:numel(conditionsOfInterest)
    nexttile(t,[1,1]);
    hold on
    plot(x_data(x_idx), stim_ts(x_idx, tmp_idx(ii)), 'Color', [.5 .5 .5], 'HandleVisibility', 'off')
    
    % se_conc  = [data_se(x_idx,stim_idx(ii),1); flipud(data_se(x_idx,stim_idx(ii),2))]';

    % hCI = fill([x_data(x_idx); flipud(x_data(x_idx))]', se_conc./maxResp, 'k', 'HandleVisibility', 'off');
    % hCI.FaceAlpha = 0.2; hCI.EdgeColor = [1 1 1];
    plot(x_data(x_idx), data(x_idx, tmp_idx(ii))./maxResp, 'k', 'LineWidth', 2)
    ylim(ybounds)

    if ii == numel(conditionsOfInterest)

        plot(x_data(x_idx), linPred(x_idx) ./ maxResp, 'k:', 'LineWidth', 2)
        legend({'Response', 'Linear prediction'}, 'Location', 'northwest')
    end
end
t.Title.String = 'Single pulse examples';
t.XLabel.String = 'Time (s)';
t.YLabel.String = 'Response magnitude';

%}
%% -- panel C: recovery from adaptation
recFig      = figure('Color', [1 1 1], 'Position', [30 300 350 400]);
set(recFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[350 400])
T1 = tiledlayout(2, 1,'TileIndexing','rowmajor');
ta = tiledlayout(T1,2,1,'TileIndexing','columnmajor');

% Select two conditions to plot
conditionsOfInterest = {'TWOPULSE-2', 'TWOPULSE-5'};
timepointsOfInterest = [-0.10 1.4];

tmp_idx        = find(contains(stim_info.name, conditionsOfInterest));
x_idx           = x_data>timepointsOfInterest(1) & x_data<=timepointsOfInterest(2);

maxResp         = max(data(x_idx,tmp_idx(1))); % scale stimulus to max of lowest duration
ybounds         = bounds(data(x_idx,tmp_idx) ./ maxResp);

% plot time courses + stimulus duration
for ii = 1:numel(conditionsOfInterest)
    nexttile(ta,[1,1])
    hold on
    plot(x_data(x_idx), stim_ts(x_idx, tmp_idx(ii)), 'Color', [.5 .5 .5], 'HandleVisibility', 'off')
    
    % se_conc  = [data_se(x_idx,stim_idx(ii),1); flipud(data_se(x_idx,stim_idx(ii),2))]';

    % hCI = fill([x_data(x_idx); flipud(x_data(x_idx))]', se_conc./maxResp, 'k', 'HandleVisibility', 'off');
    % hCI.FaceAlpha = 0.2; hCI.EdgeColor = [1 1 1];
    plot(x_data(x_idx), data(x_idx, tmp_idx(ii))./maxResp, 'k', 'LineWidth', 2)
    ylim(ybounds)

end
ta.Title.String = 'Paired pulse examples';
ta.XLabel.String = 'Time (s)';
ta.YLabel.String = 'Response magnitude';


% Compute recovery per electrode
srate           = results.channels.sampling_frequency(1);
[m, ~]         = sx_computeISIrecovery(results.data,results.t,results.stim_info,srate, 0.4, [], 'max');

% Bootstrap average parameter values
[m, se]         = bootstrWithinArea(m, @median, numBoot, CIrange);
ybounds         = bounds(se);

switch xaxisScale
    case 'log'

        % panel B - single pulses
        tb = tiledlayout(T1,1,3,'TileIndexing','columnmajor');
        tb.Layout.Tile = 2;

        nexttile(tb,[1 1])
        hold on
        plot(xISI(1), [1 1], 'k:', 'LineWidth', 2);
        plot(xISI(1)' .* [1; 1], se(1,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
        plot(xISI(1),  m(1), '.k', 'MarkerSize', 25)
        set(gca, 'TickDir', 'out', 'XTick', 0, 'XTickLabel', 0, ...
                'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
                'YTick', ybounds(1):1:ybounds(2), 'YLim', ybounds, 'XLim', [-0.005 0.01]);

        box off
        
        nexttile(tb,[1 2])
        hold on
        plot(xISI([2 end]), [1 1], 'k:', 'LineWidth', 2);
        plot(xISI(2:end)' .* [1; 1], se(2:end,:)', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
        plot(xISI(2:end),  m(2:end), '.k', 'MarkerSize', 25)
        set(gca, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
                'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
                'YTick', ybounds(1):1:ybounds(2), 'YLim', ybounds, 'XLim', [0.04 1.3]);


%         plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
%         plot(xISI(2:end)' .* [1; 1], sumResp_se(pair_idx(2:end),:)' ./ maxSum, 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
%         plot(xISI(2:end)',  sumResp(pair_idx(2:end)) ./ maxSum, '.k', 'MarkerSize', 25)
%         set(gca, 'xscale', 'log', 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
%             'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
%             'YTick', ybounds(1):1000:ybounds(2));  xlim([0.1 1])
%         legend({'summed Resp', 'linear pred'})
        box off
        tb.YLabel.String = {'Ratio second stimulus /', 'first stimulus'};
        tb.Title.String = 'Recovery from adaptation';
        tb.XLabel.String = 'Inter stimulus interval (s)';

    case 'lin'

        % recovery from adaptation
        tc = tiledlayout(T1,1,1,'TileIndexing','columnmajor');
        tc.Layout.Tile       = 2;
        tc.YLabel.String = {'Ratio second stimulus /', 'first stimulus'};
        tc.XLabel.String = 'Inter stimulus interval (s)';
        tc.Title.String = 'Recovery from adaptation';

        nexttile(tc,[1 1])
        hold on
        plot(xISI([1 end]), [1 1], 'k:', 'LineWidth', 2);
        plot(xISI' .* [1; 1], se', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
        plot(xISI,  m, '.k', 'MarkerSize', 25)
        set(gca, 'TickDir', 'out', 'XTick', 0:0.2:1.2, 'XTickLabel', 0:0.2:1.2, ...
            'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1); 
        xlim([-0.1 1.2])

end

if saveFig > 0
    print(subAddFig, fullfile(figDir, sprintf('panelC_ECoGadaptationRecovery')), '-dpdf')
end


end