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

figure; hold on
set(gcf, 'Position', [0 0 1600 500]);
% Subplot positions: % [left bottom width height]
pos(1,:) = [0.05 0.62 0.9 0.30];
pos(2,:) = [0.05 0.12 0.9 0.30];

ylabel('Broadband power change (X-fold)','FontSize', 30);

d = data(t_ind,:,ii);
maxresp = max(max(d(:))); % scale stimulus to max across conditions and dataset
set(0, 'DefaultAxesFontName', 'Helvetica Neue');
set(0, 'DefaultTextFontName', 'Helvetica Neue');

% Create a 2x6 grid of subplots
nRows = 2;
nCols = 6;

% Find all unique ONEPULSE and TWOPULSE conditions and their numbers
onepulse_idx = find(contains(conditionsOfInterest, 'ONEPULSE'));
twopulse_idx = find(contains(conditionsOfInterest, 'TWOPULSE'));

% Extract the numbers at the end of each condition string
get_col = @(str) sscanf(regexp(str, '\d+$', 'match', 'once'), '%d');

% Precompute R2 values for legend
R2val = nan(nModels, nRows, nCols);

for jj = 1:length(conditionsOfInterest)
    condStr = conditionsOfInterest{jj};
    % Determine row and column
    if contains(condStr, 'ONEPULSE')
        row = 1;
    elseif contains(condStr, 'TWOPULSE')
        row = 2;
    else
        continue; % skip if not recognized
    end
    col = get_col(condStr);
    if isempty(col) || col < 1 || col > nCols
        continue; % skip if column not valid
    end
    subplot_idx = (row-1)*nCols + col;
    subplot(nRows, nCols, subplot_idx); cla; hold on

    % Find trials for this condition
    inx = contains(stim_info.name, condStr);

    % plot stimulus
    h = plot(flatten(stim_ts(t_ind,inx))*maxresp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    if isfield(h, 'Annotation')
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end

    % plot data
    dat = flatten(d(:,inx)); 
    plot(smooth(dat,15), 'Color', 'k', 'LineWidth', 2);

    % plot models
    for kk = 1:nModels
        pred = flatten(results(kk).pred(t_ind,inx,ii));
        plot(smooth(pred,15),'Color', brclt(kk,:), 'LineWidth', 3);
        R2val(kk, row, col) = results(kk).R2.concat_all(ii);
    end

    set(gca, 'LineWidth', 2, 'FontSize', 18, 'TickDir', 'out');

    % Set x-ticks and labels
    set(gca, 'XTick',1:size(d,1):length(find(inx))*size(d,1));
    if row == 1
        set(gca, 'XTickLabel', stim_info.duration(inx))
        xlabel('Duration (s)','FontSize',18)
        if col == 1
            ylabel('ONEPULSE','FontSize',18,'FontWeight','bold');
        end
    else
        set(gca, 'XTickLabel', stim_info.ISI(inx))
        xlabel('Stimulus ISI (s)', 'FontSize', 18)
        if col == 1
            ylabel('TWOPULSE','FontSize',18,'FontWeight','bold');
        end
    end

%     l = {'Neural response',...
%         sprintf('Linear prediction, R^2 = %.2f', R2val(1)),...
%         sprintf('Delayed normalization prediction, R^2 = %.2f', R2val(2))};
%     if jj == 2, legend(l); end
    
end

sprintf('R^2 = %.2f, %.2f', R2val)
%% save

figureName = sprintf('fits_%s', channels.name{ii});
saveDir = fullfile(tt_bidsRootPath, 'derivatives', 'modelFit', 'VSSfigure');
if ~exist(saveDir, 'dir'), mkdir(saveDir), end
    print(gcf, fullfile(saveDir, figureName), '-dsvg');
close;