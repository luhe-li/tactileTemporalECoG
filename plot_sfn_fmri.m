% Run from general Code folder

clearvars; close all

%% Set up parameters
plotIndvFit     = false;
saveFig         = false;

%-- Experimental details
projectName     = 'TemporalTactileCounting';
savestr         = 'bads';
modelStr        = 'modelOutput';
models          = {'NORM', 'HRF'}; %{'NORM', 'HRF'};

ROInames        = {'localizerROI-S1'};
numROIs         = length(ROInames);

CIrange         = 5; % to compute 95% CI

%-- Setup 
bounds          = @(x) [floor(min(x(:))*10) ceil(max(x(:))*10)]/10;

%-- Set up paths
projectDir      = pwd;
dataRootDir     = fullfile('/Volumes', 'server', 'Projects', projectName,'Data');

figRoot         = fullfile(dataRootDir, '..', 'Figures');
figDir          = fullfile(figRoot, 'SFN',modelStr);
if ~exist(fullfile(figDir), 'dir'), mkdir(fullfile(figDir)); end

%-- Find participants
fnameList       = fullfile(dataRootDir, sprintf('participants.tsv'));
assert(isempty(fnameList) == 0, 'participants.tsv not found, verify paths')
T               = readtable(fnameList, 'FileType', 'text');
numSubjects     = numel(T.participant_id);
subjNames       = T.participant_id;

%-- Preallocate variables
allResults      = struct('ROI', cell(1,numROIs), ...
    'models', cell(1,numROIs), ...
    'resp', cell(1,numROIs), ...
    'c_resp', cell(1,numROIs), ...
    'sumResp', cell(1,numROIs));


%% -- Loop over all ROIs
for roi = 1:numROIs

    fprintf('[%s] Extracting model predictions for ROI %s ... \n', mfilename, ROInames{roi})

    allResults(roi).ROI     = ROInames{roi};
    allResults(roi).models  = models;

    subject         = 'group';
    fprintf('[%s] sub-%s ... \n', mfilename, subject)

    %-- create figure handles
    stimuliFig      = figure('Color', [1 1 1], 'Position', [30 300 550 300]);
    set(stimuliFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[550 300])
    summationFig    = figure('Color', [1 1 1], 'Position', [30 300 500 500]);
    set(summationFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[500 500])
    tcourseFig      = figure('Color', [1 1 1], 'Position', [30 300 700 390]);
    set(tcourseFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[700 390])
    paramFig        = figure('Color', [1 1 1], 'Position', [30 300 700 250]);
    set(paramFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize', [700 250])
    summaryFig      = figure('Color', [1 1 1], 'Position', [30 300 550 260]);
    set(summaryFig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[550 260])
    T = tiledlayout(1, 3,'TileIndexing','rowmajor');

    for wModel = 1:numel(models)

        model           = models{wModel};
        modelSett       = visualizationSettings(models, model);

        %-- Necessary paths
        resultsDir      = fullfile(dataRootDir, modelStr, model);

        %-- Results file names (bootstrapped, avg fits, and cross val R2)
        fileNames       = dir(fullfile(resultsDir, sprintf('sub-%s_%s_model-%s_crossval-*_optimizer-%s_*%s.mat', subject, ROInames{roi}, model, savestr, modelStr)));

        %-- Skip if not all model output files exist
        if length(fileNames) < 3

            fprintf('[%s] No or not all model results found: skipping %s - %s \n', mfilename, subject, ROInames{roi})
            continue

        end

        %-- Load results
        clearvars results
        for n = 1:numel(fileNames)
            if contains(fileNames(n).name, 'bts')
                % bootstrap, nocross results
                results(2) = load(fullfile(resultsDir, fileNames(n).name), 'subject', 'model', 'currModel', 'ROIname', 'condOrder', 'x_data', 'y_data', 'y_est', 'params', 'R2');
            elseif contains(fileNames(n).name, 'withCross')
                % fit average, cross results
                results(3) = load(fullfile(resultsDir, fileNames(n).name), 'subject', 'model', 'currModel', 'ROIname', 'condOrder', 'x_data', 'y_data', 'y_est', 'params', 'R2');
            else
                % fit average results
                results(1) = load(fullfile(resultsDir, fileNames(n).name), 'subject', 'model', 'currModel', 'ROIname', 'condOrder', 'x_data', 'y_data', 'y_est', 'params', 'R2');
            end
        end

        % combine one-pulse conditions
        onePulseIndx    = cat(1, find(contains(results(1).condOrder, 'BLANK')), ...
            find(contains(results(1).condOrder, 'ONE')));

        % combine two-pulse conditions
        twoPulseIndx    = cat(1, find(contains(results(1).condOrder, 'ONE-PULSE-4')), ...
            find(contains(results(1).condOrder, 'TWO')));

        allPulsesIndx   = cat(1, onePulseIndx, twoPulseIndx);

        %-- Preallocate variables
        ybounds         = bounds(prctile([results(2).y_data, results(2).y_est], [1 99], [1 3]));
        if mod(ybounds(1), 0.2) > 0; ybounds(1) = ybounds(1)-0.1; end
        if mod(ybounds(2), 0.2) > 0; ybounds(2) = ybounds(2)+0.1; end
        Responses       = NaN(modelSett.numCond, length(results(1).x_data));
        cResponses68    = NaN(modelSett.numCond, length(results(1).x_data), 2);
        cResponses95    = NaN(modelSett.numCond, length(results(1).x_data), 2);
        %modelPrediction = NaN(modelSett.numCond, length(results(1).x_data));
        btstrSumResp    = NaN(modelSett.numCond, size(results(2).y_data,3));

        % create finer sampling for model prediction
        out             = createSmoothPrediction(results(1).params, results(1).model);
        modelPrediction = out.pred; 

        % %-- Extract data for all pulse conditions
        % for ii = 1:modelSett.numCond

        %     % time series data & predictions
        %     Responses(ii,:)         = results(1).y_data(:,allPulsesIndx(ii));
        %     btstrSumResp(ii,:)      = sum(squeeze(results(2).y_data(:,allPulsesIndx(ii),:)),1);
        %     %modelPrediction(ii,:)   = results(1).y_est(:,allPulsesIndx(ii));

        %     % confidence interval data based on bootstrapped data
        %     cResponses95(ii,:,:)    = prctile(squeeze(results(2).y_data(:,allPulsesIndx(ii),:)), [0 100] + (5/2 * [1 -1]), 2);
        %     cResponses68(ii,:,:)    = prctile(squeeze(results(2).y_data(:,allPulsesIndx(ii),:)), [0 100] + (32/2 * [1 -1]), 2);

        %     %-- Visualize time courses for all pulse conditions
        %     figure(tcourseFig)
        %     subplot(2,numel(onePulseIndx),ii)
        %     hold on,
        %     if wModel == 1
        %         plot([0 results(1).x_data(end)], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
        %         hCI = fill([results(1).x_data, fliplr(results(1).x_data)], [squeeze(cResponses95(ii,:,1)), fliplr(squeeze(cResponses95(ii,:,2)))], 'k', 'HandleVisibility', 'off');
        %         hCI.FaceAlpha = 0.2; hCI.EdgeColor = [1 1 1];
        %         plot(results(1).x_data, Responses(ii,:)', '.k-', 'LineWidth', 2, 'MarkerSize', 15, 'HandleVisibility', 'off')
        %         if ii <= numel(onePulseIndx)
        %             title(sprintf('Dur %.2fs', modelSett.stimDur(ii)), 'FontSize', 10)
        %         else
        %             title(sprintf('ISI %.2fs', modelSett.stimDur(ii-numel(onePulseIndx))), 'FontSize', 10)
        %         end
        %         set(gca,'TickDir', 'out', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1)

        %     end
        %     if ii <= numel(onePulseIndx)
        %         stimIndx  = find(out.stimdur == modelSett.stimDur(ii));
                
        %     else
        %         stimIndx  = find(out.stimdur == modelSett.stimDur(ii-numel(onePulseIndx))) + numel(out.stimdur);
        %     end
        %     plot(results(1).x_data, modelPrediction(stimIndx,:), 'Color', modelSett.color, 'LineWidth', 2)
        %     ylim(ybounds)
        %     set(gca, 'ytick', ybounds(1):0.4:ybounds(2), 'tickdir', 'out')
        %     box off

        % end

        % % Create figure title
        % if wModel == 1
        %     titleStr    = sprintf('%s %s \n %s{%f %f %f} %smodel crossval R2: %.2f \n', subject, ROInames{roi}, '\color[rgb]', modelSett.color, model, results(3).R2);
        % else
        %     titleStr    = cat(2, titleStr, sprintf('%s{%f %f %f} %smodel crossval R2: %.2f \n', '\color[rgb]', modelSett.color, model, results(3).R2));
        % end

        % %-- Visualize reliability of param estimates
        % figure(paramFig)
        % c = 1;
        % axOrder     = [5 6 1 2 3 4];

        % for ii = 1:modelSett.totParams

        %     subplot(1,6,axOrder(ii))

        %     % skip plotting if current model does not have this param
        %     if strcmp(modelSett.allLabels{ii}, modelSett.mLabels{c})

        %         hold on,
        %         % median with 68 and 95 CI intervals
        %         plot(wModel * ones(1,2), prctile(results(2).params(c,:), [0 100] + (5/2 * [1 -1])), 'Color', [0.8 0.8 0.8], 'LineWidth', 3)
        %         plot(wModel * ones(1,2), prctile(results(2).params(c,:), [0 100] + (32/2 * [1 -1])), 'Color', modelSett.color, 'LineWidth', 3)
        %         scatter(wModel, median(results(2).params(c,:)), 100, modelSett.color, 'filled')

        %         c = c+1;

        %     end

        %     if wModel == numel(models)
        %         box off; axis square
        %         title(modelSett.labels{ii}, 'Interpreter', 'none')
        %         % ensure 2 gammas have same axes
        %         if strcmp(modelSett.allLabels{ii}, 'gamma2')
                    
        %             oldLim = currLim;
        %             currLim = get(gca, 'YLim'); currLim = [floor(currLim(1)) ceil(currLim(2))];

        %             currLim(1) = floor(min(oldLim(1), currLim(1)));
        %             currLim(2) = ceil(max(oldLim(2), currLim(2)));
        %             subplot(1,6,axOrder(ii-1))
        %             ylim(currLim)
        %             set(gca, 'xtick', 1:numel(models), 'xticklabel', models, ...
        %                     'tickdir', 'out', 'ytick', linspace(currLim(1),currLim(2), 5));
        %             subplot(1,6,axOrder(ii))
        %             ylim(currLim)
        %             set(gca, 'xtick', 1:numel(models), 'xticklabel', models, ...
        %                     'tickdir', 'out', 'ytick', linspace(currLim(1),currLim(2), 5));
        %         else
        %             currLim = get(gca, 'YLim'); currLim =  [floor(currLim(1)*10)/10 ceil(currLim(2)*10)/10];
        %             ylim(currLim)
        %             set(gca, 'xtick', 1:numel(models), 'xticklabel', models, ...
        %                     'tickdir', 'out', 'ytick', linspace(currLim(1),currLim(2), 5));
        %         end
        %         %                 ylim([modelSett.mBounds(1,ii) currLim(2)])
        %         xlim([0 numel(models)+1])
               
        %     end
        % end

        %-- summed response figure
        figure(summaryFig)

        if wModel == 1

            xvalues         = modelSett.stimDur;
            ybounds = [-1 5];
            
            %-- one pulse
            t1 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
            t1.Layout.Tile = 1;
            ax1 = nexttile(t1,[1 1]);
            hold on,
            plot([-0.1 0.1], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
            plot(0 .* [1; 1], squeeze(sum(cResponses95(1,:,:),2))', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(0 .* [1; 1], squeeze(sum(cResponses68(1,:,:),2))', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(0,  sum(Responses(1, :),2), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')
          
            ax2 = nexttile(t1,[1 2]);
            hold on
            plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
            plot(xvalues(2:end) .* [1; 1], squeeze(sum(cResponses95(2:numel(modelSett.stimDur),:,:),2))', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(xvalues(2:end) .* [1; 1], squeeze(sum(cResponses68(2:numel(modelSett.stimDur),:,:),2))', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(xvalues(2:end),  sum(Responses(2:numel(modelSett.stimDur), :),2), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')
            
            %-- paired pulse             
            t2 = tiledlayout(T,1,3,'TileIndexing','columnmajor');
            t2.Layout.Tile = 2;
            ax3 = nexttile(t2,[1 1]);
            hold on
            plot([-0.1 0.1], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
            plot(0 .* [1; 1], squeeze(sum(cResponses95((1)+numel(onePulseIndx),:,:),2))', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(0 .* [1; 1], squeeze(sum(cResponses68((1)+numel(onePulseIndx),:,:),2))', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(0,  sum(Responses((1)+numel(onePulseIndx), :),2), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')
            
            ax4 = nexttile(t2,[1 2]);
            hold on
            plot([0.04 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
            plot(xvalues(2:end) .* [1; 1], squeeze(sum(cResponses95((2:numel(modelSett.stimDur))+numel(modelSett.stimDur),:,:),2))', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(xvalues(2:end) .* [1; 1], squeeze(sum(cResponses68((2:numel(modelSett.stimDur))+numel(modelSett.stimDur),:,:),2))', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(xvalues(2:end),  sum(Responses((2:numel(modelSett.stimDur))+numel(modelSett.stimDur), :),2), '.k', 'MarkerSize', 25,  'HandleVisibility', 'off')


            %-- linear xaxis
            t3 = tiledlayout(T,2,1,'TileIndexing','columnmajor');
            t3.Layout.Tile = 3;
            ax5 = nexttile(t3, [1 1]);
            hold on
            plot([-0.05 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
            plot(xvalues .* [1; 1], squeeze(sum(cResponses95(1:numel(modelSett.stimDur),:,:),2))', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(xvalues .* [1; 1], squeeze(sum(cResponses68(1:numel(modelSett.stimDur),:,:),2))', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(xvalues,  sum(Responses(1:numel(modelSett.stimDur), :),2), '.k', 'MarkerSize', 10,  'HandleVisibility', 'off')
            
            ax6 = nexttile(t3,[1 1]);
            hold on
            plot([-0.05 1.3], [0 0], 'k', 'LineWidth', 1, 'HandleVisibility','off')
            plot(xvalues .* [1; 1], squeeze(sum(cResponses95((1:numel(modelSett.stimDur))+numel(modelSett.stimDur),:,:),2))', 'Color', [0.8 0.8 0.8], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(xvalues .* [1; 1], squeeze(sum(cResponses68((1:numel(modelSett.stimDur))+numel(modelSett.stimDur),:,:),2))', 'Color', [0 0 0], 'LineWidth', 2, 'HandleVisibility', 'off')
            plot(xvalues,  sum(Responses((1:numel(modelSett.stimDur))+numel(modelSett.stimDur), :),2), '.k', 'MarkerSize', 10,  'HandleVisibility', 'off')


        end

        onePulsePred        = modelPrediction(1:numel(out.stimdur), :);
        twoPulsePred        = modelPrediction(numel(out.stimdur)+1:end, :);
        plot(ax1, out.stimdur(out.stimdur < 0.01), sum(onePulsePred((out.stimdur < 0.01), :),2), 'Color', modelSett.color, 'LineWidth', 2)
        plot(ax2, out.stimdur(out.stimdur >= xvalues(2)-0.01), sum(onePulsePred((out.stimdur >= xvalues(2)-0.01), :),2), 'Color', modelSett.color, 'LineWidth', 2)
        
        plot(ax3, out.stimdur(out.stimdur < 0.01), sum(twoPulsePred((out.stimdur < 0.01), :),2), 'Color', modelSett.color, 'LineWidth', 2)
        plot(ax4, out.stimdur(out.stimdur >= xvalues(2)-0.01), sum(twoPulsePred((out.stimdur >= xvalues(2)-0.01), :),2), 'Color', modelSett.color, 'LineWidth', 2)

        plot(ax5, out.stimdur, sum(onePulsePred,2), 'Color', modelSett.color, 'LineWidth', 1)
        plot(ax6, out.stimdur, sum(twoPulsePred,2), 'Color', modelSett.color, 'LineWidth', 1)
        
        if wModel == numel(models)
            
            ybounds = get(ax2, 'YLim');

            set(ax1, 'TickDir', 'out', 'XTick', 0, 'XTickLabel', 0, ...
                'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
                'YTick', ybounds(1):1:ybounds(2), 'YLim', ybounds, 'XLim', [-0.005 0.01]);
            ax2.Box = 'off';
            t1.Title.String = 'Single pulse conditions';
            t1.XLabel.String = 'Stimulus duration (s)';
            t1.YLabel.String = {'Summed BOLD time series';'(%SC)'};
            

            set(ax2, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
                'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
                'YTick', ybounds(1):1:ybounds(2), 'YLim', ybounds, 'XLim', [0.04 1.3]);
            ax2.Box = 'off';
            ax2.YAxis.Visible = 'off';
          
            set(ax3, 'TickDir', 'out', 'XTick', 0, 'XTickLabel', 0, ...
                'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
                'YTick', ybounds(1):1:ybounds(2), 'YLim', ybounds, 'XLim', [-0.005 0.01]);
            ax3.Box = 'off';
            t2.Title.String = 'Paired pulse conditions';
            t2.XLabel.String = 'Interstimulus interval (s)';
            t2.YLabel.String = {'Summed BOLD time series';'(%SC)'};
            
            set(ax4, 'TickDir', 'out', 'XTick', [0.1 1], 'XTickLabel', [0.1 1], ...
                'Xscale', 'log', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
                'YTick', ybounds(1):1:ybounds(2), 'YLim', ybounds, 'XLim', [0.04 1.3]);
            ax4.Box = 'off';
            ax4.YAxis.Visible = 'off';
             

            set(ax5, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
                'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
                'YTick', ybounds(1):1:ybounds(2), 'YLim', ybounds, 'XLim', [-0.1 1.3]);
            ax5.Box = 'off';
            ax5.Title.String = 'Single pulse conditions';
            ax5.XLabel.String = 'Duration (s)';
            ax5.DataAspectRatio = [1 4 1];
            t3.YLabel.String = {'Summed BOLD time series';'(%SC)'};

            set(ax6, 'TickDir', 'out', 'XTick', 0:0.4:1.2, 'XTickLabel', 0:0.4:1.2, ...
                'FontSize', 8, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, ...
                'YTick', ybounds(1):1:ybounds(2), 'YLim', ybounds, 'XLim', [-0.1 1.3]);
            ax6.Box = 'off';
            ax6.Title.String = 'Paired pulse conditions';
            ax6.DataAspectRatio = [1 4 1];
            ax6.XLabel.String = 'ISI (s)';

            title(T, sprintf('Summed response: %s %s', subject, ROInames{roi}), 'Fontsize', 12)

            if saveFig > 0
                print(summaryFig, fullfile(figDir, sprintf('sub-%s_%s_sumResponses', subject, ROInames{roi})), '-dpdf')
            end
        end


        %-- Save responses for all participants [cond tr sub]
        if wModel == 1
            allResults(roi).resp        = Responses;
            allResults(roi).cResp95     = cResponses95;
            allResults(roi).cResp68     = cResponses68;
            allResults(roi).sumResp     = sum(Responses,2);
            allResults(roi).btstrSumResp = btstrSumResp;
        end

        %-- Extract model params and CI
        allResults(roi).(sprintf('%spred', model))      = cat(1, onePulsePred(ismember(out.stimdur, modelSett.stimDur), :), ...
                                                                 twoPulsePred(ismember(out.stimdur, modelSett.stimDur), :));
        allResults(roi).(sprintf('sum%spred', model))   = sum(allResults(roi).(sprintf('%spred', model)),2);
        allResults(roi).(sprintf('%sparams', model))    = results(1).params;
        allResults(roi).(sprintf('c%sparams', model))   = prctile(results(2).params, [5 95], 2);
        allResults(roi).(sprintf('%scrossR2', model))   = results(3).R2;

        allResults(roi).(sprintf('mse%s', model))       = mean((Responses - allResults(roi).(sprintf('%spred', model))).^2,2);

    end

    % %-- wrap up timecourses figure
    % figure(tcourseFig)
    % subplot(2,numel(onePulseIndx),modelSett.numCond)
    % legend(models)

    % sgtitle(titleStr, 'fontsize', 20)

    % if saveFig > 0
    %     print(tcourseFig, fullfile(figDir, sprintf('sub-%s_%s_timeCourses', subject, ROInames{roi})), '-dpdf')
    % end

    % %-- wrap up parameter figure
    % figure(paramFig)
    % sgtitle(sprintf('%s %s \n Bootstrapped parameter estimates \n', subject, ROInames{roi}), 'Fontsize', 20)
    % colororder(modelSett.mColors)
    % if saveFig > 0
    %     print(paramFig, fullfile(figDir, sprintf('sub-%s_%s_modelParams', subject, ROInames{roi})), '-dpdf')
    % end

    % %-- create stimuli figure
    % figure(stimuliFig)
    % stimSeqs    = createStimSeq(modelSett, stimuliFig);
    % if saveFig > 0
    %     print(stimuliFig, fullfile(figDir, sprintf('sub-%s_%s_stimuliSeq', subject, ROInames{roi})), '-dpdf')
    % end

    % %-- create linear system prediction
    % figure(summationFig)
    % allResults.condNames    = modelSett.condNames;
    % allResults.x_data       = modelSett.x_data;
    % allResults.stimDur      = modelSett.stimDur;
    % createLinearSystemPredictions(allResults, summationFig);

    if saveFig > 0
        print(summationFig, fullfile(figDir, sprintf('sub-%s_%s_temporalLinearity', subject, ROInames{roi})), '-dpdf')
    end


end


