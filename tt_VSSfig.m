% VSS figure: compare DN fits and linear fits to the data of one electrode

clear; close all;

modelfuns = {@LINEAR, @DN};
xvalmode = 0;
datatype = 'individualelecs';
brclt = [30, 120, 180; 227, 27, 27]./255;

%% load data and model fits
for zz = 1:length(modelfuns)
    [allD(zz)] = tt_loadDataForFigure(modelfuns{zz}, xvalmode, datatype);
end
[results] = tt_evaluateModelFit(allD,0);

D = allD(1);

%% extract info

data = D.data;
channels = D.channels;
stim_ts = D.stim;
stim_info = D.stim_info;
t = D.t;
conditionsOfInterest = {'ONEPULSE', 'TWOPULSE'};
timepointsOfInterest = [t(1) t(end)];

nModels = length(modelfuns);
nDatasets   = size(data,3);
nCond       = length(conditionsOfInterest);
t_ind       = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);
ii          = nDatasets; % select the last channel

figure;
set(gcf, 'Position', [0 0 2400 700]);
set(0, 'DefaultAxesFontName', 'Helvetica Neue');
set(0, 'DefaultTextFontName', 'Helvetica Neue');

d = data(t_ind,:,ii);
maxresp = max(max(d(:))); % scale stimulus to max across conditions and dataset

nRows = 2; nCols = 6;

% Precompute R2 for legend if needed
R2val = nan(1, nModels);
for kk = 1:nModels
    R2val(kk) = results(kk).R2.concat_all(ii);
end

for idx = 1:length(stim_info.name)
    name = stim_info.name{idx};
    % Determine row and column
    if contains(name, 'ONEPULSE')
        row = 1;
        % Extract the number at the end of the name for column
        colnum = regexp(name, '(\d+)$', 'tokens');
        if isempty(colnum)
            continue; % skip if no number found
        end
        col = str2double(colnum{1}{1});
    elseif contains(name, 'TWOPULSE')
        row = 2;
        colnum = regexp(name, '(\d+)$', 'tokens');
        if isempty(colnum)
            continue;
        end
        col = str2double(colnum{1}{1});
    else
        continue; % skip if not one of the two conditions
    end
    if col > nCols
        continue; % skip if column out of range
    end
    subplot(nRows, nCols, (row-1)*nCols + col); hold on; box off

    % Plot stimulus
    plot(t(2:end), stim_ts(t_ind,idx)*maxresp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

    % Plot data
    plot(t(2:end), smooth(d(:,idx),15), 'Color', 'k', 'LineWidth', 2);

    % Plot model predictions
    for kk = 1:nModels
        pred = results(kk).pred(t_ind,idx,ii);
        plot(t(2:end),smooth(pred,15), 'Color', brclt(kk,:), 'LineWidth', 3);
    end

    set(gca, 'LineWidth', 2, 'FontSize', 16, 'TickDir', 'out');
    % if row == nRows
    %     xlabel('Time (samples)', 'FontSize', 18);
    % end
    if col == 1
        ylabel('Broadband power change (X-fold)', 'FontSize', 18);
    else
        set(gca, 'YTick', [], 'YTickLabel', []);
    end

    % Title with condition info
    if contains(name, 'ONEPULSE')
        ttl = sprintf('ONEPULSE %ds', stim_info.duration(idx));
    else
        ttl = sprintf('TWOPULSE ISI %ds', stim_info.ISI(idx));
    end

    % Optionally, add legend to first subplot
    % if (row == 1 && col == 1)
    %     l = {'Stimulus', 'Neural response', ...
    %         sprintf('Linear pred, R^2=%.2f', R2val(1)), ...
    %         sprintf('DN pred, R^2=%.2f', R2val(2))};
    %     legend(l, 'Location', 'best');
    % end
    % xlim(timepointsOfInterest)
    
end

sprintf('R^2 = %.2f, %.2f', R2val)
%% save

figureName = sprintf('fits_%s', channels.name{ii});
% saveDir = fullfile(tt_bidsRootPath, 'derivatives', 'modelFit', 'VSSfigure');
saveDir = fullfile(pwd, 'figures');
if ~exist(saveDir, 'dir'), mkdir(saveDir), end
    print(gcf, fullfile(saveDir, figureName), '-dsvg');
close;

timepointsOfInterest = [0, 2];
tt_plotSumDataAndFits(results, D.data, D.channels, D.stim, D.stim_info, D.t, D.options, saveDir, {'ONEPULSE', 'TWOPULSE'},timepointsOfInterest)
