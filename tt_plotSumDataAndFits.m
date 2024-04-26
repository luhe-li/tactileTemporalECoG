function tt_plotSumDataAndFits(results, data, channels, stim_ts, stim_info, t, options, saveDir, conditionsOfInterest, timepointsOfInterest)

% Will generate plots of summed response for data (in black) and model predictions (in colours).

if ~exist('saveDir', 'var') || isempty(saveDir), saveDir = []; end
if ~exist('conditionsOfInterest', 'var') || isempty(conditionsOfInterest)
    conditionsOfInterest = {'ONEPULSE', 'TWOPULSE'};
end
if ~exist('timepointsOfInterest', 'var') || isempty(timepointsOfInterest)
    timepointsOfInterest = [t(1) t(end)];
end

nModels     = size(results,2);
nDatasets   = size(data,3);
nCond       = length(conditionsOfInterest);
%stim_info   = stim_info(contains(stim_info.name, conditionsOfInterest),:);
t_ind       = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

% Determine if data was averaged across elecs prior to fit
if options.average_elecs
    dataWasAveraged = true;
else
    dataWasAveraged = false;
end

%% Plot summed data and predictions
colors = {'r', 'b', 'c', 'm', 'g', 'y'}; % assuming we'll never plot >6 model fits at a time

% Prepare legend
l = cell(1,nModels+1);
l{1} = 'data';
for kk = 1:nModels, l{kk+1} = func2str(results(kk).model); end

% Average and get s.d. of channels or channel averages
figure;
set(gcf, 'Position', [0 0 1100 400]);
t = tiledlayout(1, nCond);
ylabel(t, 'Summed broadband power (0-1s)','FontSize', 20);

d = data(t_ind,:,:);

% Loop over conditions
for jj = 1:length(conditionsOfInterest)

    nexttile
    hold on
    axis tight
    inx = contains(stim_info.name, conditionsOfInterest{jj});
    cond = unique(stim_info.condition(inx));

    if cond == 1
        xlabel('Stimulus duration (s)');
        lin_x = stim_info.duration(inx)';
        ylim([0, 1.5])
    elseif cond == 2
        xlabel('Stimulus ISI (s)');
        lin_x = stim_info.ISI(inx)';
        ylim([0, 4]);
    end

    % get log x
    log_x = log10(cat(2, lin_x(2)/2, lin_x(2:end)));
    x = log_x;
    xbounds = [x(1) x(end)] + diff(x(1:2)) * [-.5 .5];
    xlim(xbounds)

    % plot data
    slc_d = squeeze(sum(d(:,inx,:),1));
    m = mean(slc_d, 2, 'omitnan');
    se = squeeze(prctile(slc_d,[15.87 84.13],2));
    tt_plotPoints(m, se, x, 'errbar', 1,[],50);

    % plot DN model prediction
    kk=1;
    pred = squeeze(sum(results(kk).pred(t_ind,inx,:),1));
    m = mean(pred, 2, 'omitnan');
    se = squeeze(prctile(pred,[15.87 84.13],2));
    tt_plotPoints(m, se, x, 'ci', 1,[], 50, colors{kk});

    set(gca, 'XTick',x, 'XTickLabel', lin_x,'xticklabelrotation', 45);
    set(gca, 'FontSize', 20);

    % add legend
    %if jj == 1, legend(l); end
end

% Determine how to name the plot
if ~dataWasAveraged
    figureName = sprintf('sumFits_indivelecs_%s_%s', [channels.name{:}], [l{2:end}]);
else
    figureName = sprintf('sumFits_avgelecs_%s_%s', [channels.name{:}], [l{2:end}]);
end

set(gcf, 'Name', figureName);

% Determine whether to save it and if so where
if ~isempty(saveDir)
    if ~dataWasAveraged
        figDir = fullfile(saveDir, 'individualelectrodes');
    else
        figDir = fullfile(saveDir, 'electrodeaverages');
    end
    if ~exist(figDir, 'dir'), mkdir(figDir), end
    saveas(gcf, fullfile(figDir, figureName), 'png');
    close;
end

end