function tt_plotDataAndFits(results, data, channels, stim_ts, stim_info, t, options, saveDir, conditionsOfInterest, timepointsOfInterest)

% Will generate plots of concatenated stimulus time courses (in gray),
% concatenated data time courses (in black) and model predictions (in colours).
%
% 2020 Iris Groen 

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

%% Plot data and predictions of timecourse
colors = {'r', 'b', 'c', 'm', 'g', 'y'}; % assuming we'll never plot >6 model fits at a time

% Prepare legend
l = cell(1,nModels+1);
l{1} = 'data';
for kk = 1:nModels, l{kk+1} = func2str(results(kk).model); end

% Loop over channels or channel averages
for ii = 1:nDatasets
    
    figure;
    set(gcf, 'Position', [400 200 1800 1200]);
    tt = tiledlayout(nCond, 1);
    ylabel(tt, 'X-fold increase in broadband power','FontSize', 20);

    d = data(t_ind,:,ii);
	maxresp = max(max(d(:))); % scale stimulus to max across conditions and dataset

    % Loop over conditions 
    for jj = 1:length(conditionsOfInterest)

        nexttile
        hold on
        inx = contains(stim_info.name, conditionsOfInterest{jj});
        cond = unique(stim_info.condition(inx));

        % plot stimulus
        h = plot(flatten(stim_ts(t_ind,inx))*maxresp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        % plot data
        plot(flatten(d(:,inx)), 'Color', 'k', 'LineWidth', 2); 
        titlestr = cell(1,nModels);
        
        % plot models
        for kk = 1:nModels           
            pred = results(kk).pred(t_ind,inx,ii);
            %plot(flatten(pred), 'Color', colors{kk}, 'LineStyle', '-.', 'LineWidth', 2);
            plot(flatten(pred), 'Color', colors{kk},  'LineWidth', 2);
            if isfield(results(kk).R2, 'concat_cond')
                R2val = mean(results(kk).R2.concat_cond(cond,ii));
            else
                R2val = mean(results(kk).R2.stim(inx,ii));
            end
            titlestr{kk} = sprintf('%s r2 = %0.2f   ', func2str(results(kk).model), R2val);
        end
        
        % add title
        title(sprintf('%s: %s', conditionsOfInterest{jj}, [titlestr{:}]));
        set(gca, 'XTick',1:size(d,1):length(find(inx))*size(d,1));
        if contains(conditionsOfInterest{jj}, 'ONEPULSE')
            set(gca, 'XTickLabel', stim_info.duration(inx))
            xlabel('Stimulus duration (s)')
        else
            set(gca, 'XTickLabel', stim_info.ISI(inx))
            xlabel('Stimulus ISI (s)')
        end

        set(gca, 'FontSize', 20);

        % add legend
        %if jj == 1, legend(l); end
    end

    % Determine how to name the plot 
    if ~dataWasAveraged
        figureName = sprintf('fits_%s_%s', channels.name{ii}, [l{2:end}]);
    else
        figureName = sprintf('fits_avg%s_%s', [channels.name{:}], [l{2:end}]);
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

end