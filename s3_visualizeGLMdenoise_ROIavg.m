
clearvars; close all

%% Setup params
projectName     = 'TemporalTactileCounting';

%-- Setup paths
[projectRootDir, dataRootDir] = rootPath(false, projectName);
dataRootDir     = fullfile('/data1/projects/dumoulinlab/Lab_members/Bloem/projects', projectName, 'Data');
expPrms         = expDetailsTemporalTactile(projectName);

taskNames       = expPrms.taskNames;
hrfmodel        = 'fir'; % 'optimize' or 'fir'
glmFolder       = 'corticalRibbon';

stimdur         = [0 0.05 0.1 0.2 0.4 0.8 1.2]; % seconds
twoPulseDur     = 0.2; % pulse duration in two-pulse condition, in seconds
numPulses       = numel(stimdur);
tr              = 1;

bounds          = @(x) [floor(min(x(:))*10) ceil(max(x(:))*10)]/10;

Glasserlabels   = { 'S1', [9 51 52]
                    'BA3b', 9; ...
                    'BA1', 51; ...
                    'BA2', 52; ...
                    };

allROInames     = { 'localizerROI', 'S1';
                    'localizerROI', 'BA3b';
                    'localizerROI', 'BA1';
                    'localizerROI', 'BA2';
                    };

ROInames        = allROInames(1,:);
numROIs         = size(ROInames,1);
roiColors       = parula(numROIs+1);

ConditionNames  = expPrms.ConditionNames;
condOrder       = cat(1, find(contains(ConditionNames, 'BLANK')), ...
    find(contains(ConditionNames, 'ONE')), ...
    find(contains(ConditionNames, 'ONE-PULSE-4')), ...
    find(contains(ConditionNames, 'TWO')));

numConditions   = numel(ConditionNames); %sum(contains(ConditionNames, 'PULSE'));

%% Setup directories
glmDir          = fullfile(dataRootDir, 'GLMdenoise', glmFolder);

%% find all participants
fnameList       = fullfile(dataRootDir, '..', sprintf('participants.tsv'));
T               = readtable(fnameList, 'FileType', 'text');
subjNames       = T.participant_id;
numSubjects     = numel(subjNames);

%% Preallocate variables
allFirBetas     = cell(1,size(ROInames,1));
corrFirBetas    = cell(1,size(ROInames,1));
sumBetas        = cell(1,size(ROInames,1));
corrsumBetas    = cell(1,size(ROInames,1));

%-- Where to save the figures
figDir          = fullfile(fileparts(dataRootDir), 'Figures', 'GLMoutput', glmFolder);
if ~exist(figDir, 'dir'), mkdir(figDir); end

%%
for sub = 1:numSubjects

    %-- Current subject
    subject         = subjNames{sub};

    %-- Load deconvolution GLMdenoise results
    dataFile        = dir(fullfile(glmDir, subject, hrfmodel, sprintf('results_%s_task-%s.mat', subject, taskNames{1})));
    assert(~isempty(dataFile), 'GLM results file not found')

    load(fullfile(dataFile(1).folder, dataFile(1).name));

    %-- Setup params
    lastTR          = size(results.modelmd,3);
    x_data          = 0:tr:lastTR-1;

    %-- Setup figure legend
    globalMin = Inf;
    globalMax = -Inf;

    figHandle       = figure('Color', [1 1 1], 'Position', [30+sub*10 300 1400 800]);
    set(figHandle,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[1400 800])

    for roi = 1:size(ROInames,1)

        if sub == 1
            allFirBetas{roi}    = NaN(numConditions, numel(x_data), numSubjects);
            sumBetas{roi}       = NaN(numConditions, numSubjects);
            corrFirBetas{roi}   = NaN(numConditions, numel(x_data), numSubjects);
            corrsumBetas{roi}   = NaN(numConditions, numSubjects);
        end

        %-- Load ROI to mask data
        maskFile        = dir(fullfile(dataRootDir, 'roiVols', subject, sprintf('*h.%s.VOL.nii.gz', ROInames{roi,1})));
        if ~isempty(maskFile)
            mask            = MRIread(fullfile(maskFile(1).folder, maskFile(1).name));
            maskVec         = reshape(mask.vol, [prod(results.volDims), 1]);
            roiMask         = maskVec(results.brainMask);
        else
            roiMask         = true(sum(results.brainMask,1),1);
        end

        %-- Use voxels within ROI to average
        whichLabels     = strcmp(Glasserlabels(:,1), ROInames(roi,2));
        glasserMask     = ismember(results.ROIlabels, Glasserlabels{whichLabels, 2}) & roiMask;

        ROIindx         = roiMask & glasserMask;
        coThreshold     = prctile(results.locCoherence(ROIindx), 75);
        selectedVox     = results.locCoherence >= coThreshold & ROIindx;
        firBetas        = squeeze(mean(results.modelmd(selectedVox,:,:),1)); % voxels * cond * timepoints

        %-- compute bootstrapped CI
        firBetasCI      = NaN(numConditions, numel(x_data), 2);
        for cond = 1:numConditions
            
            voxFIR                  = squeeze(results.modelmd(selectedVox,cond,:));
            firBetasCI(cond,:,:)    = bootci(100, @mean, voxFIR)'; %bootci samples from rows

        end

        %-- combine across observers
        allFirBetas{roi}(:,:,sub)   = firBetas;

        %-- sum responses
        sumBetas{roi}(:,sub)        = sum(firBetas,2);

        %-- save out some stats about voxel selection
        fname           = sprintf('%s_%s_GLM%s_%s_voxelStats.json', projectName, subject, ROInames{roi,2}, glmFolder);
        inputVar        = struct('scriptName', mfilename, 'subId', subject, ...
                            'roiName', ROInames(roi,2), 'glasserLabels', Glasserlabels(whichLabels, 2), ...
                            'nVoxlocROI', sum(roiMask), 'glasserLabelslocROI', unique(results.ROIlabels(roiMask>0)), ...
                            'nVoxlocROIbyGlasserLabel', histc(results.ROIlabels(roiMask > 0), unique(results.ROIlabels(roiMask>0))), ...
                            'nVoxGlasser', sum(ROIindx), 'nVoxthresCo', sum(selectedVox), ...
                            'coThreshold', coThreshold, ...
                            'nVoxthresCobyGlasserLabel', histc(results.ROIlabels(selectedVox > 0), unique(results.ROIlabels(selectedVox>0))));        
        savejson('',inputVar,fullfile(dataFile.folder,fname));

        %-- Plot avg FIR      
        ybounds = bounds(firBetasCI);
        globalMin = min(globalMin, ybounds(1));
        globalMax = max(globalMax, ybounds(2));

        roiHandle       = figure('Color', [1 1 1], 'Position', [30+sub*10 300 1400 800]);
        set(roiHandle,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[1400 800])

        for cond = 1:numel(condOrder)
    
            cBoundsFIR      = squeeze(firBetasCI(condOrder(cond),:,:))';

            figure(figHandle)
            subplot(2, 7, cond);
            hold on;
            h = fill([x_data, fliplr(x_data)], [cBoundsFIR(1,:), fliplr(cBoundsFIR(2,:))], roiColors(roi,:), 'HandleVisibility', 'off');
            h.FaceAlpha = 0.2;   
            plot(x_data, firBetas(condOrder(cond), :), '.-', 'Color', roiColors(roi,:), 'LineWidth', 2, 'MarkerSize', 20);
            box off;
            ylim([globalMin, globalMax]);
            title(ConditionNames{condOrder(cond)}, 'Interpreter','none');

            figure(roiHandle)
            subplot(2, 7, cond);
            hold on;
            h = fill([x_data, fliplr(x_data)], [cBoundsFIR(1,:), fliplr(cBoundsFIR(2,:))], roiColors(roi,:), 'HandleVisibility', 'off');
            h.FaceAlpha = 0.2;   
            plot(x_data, firBetas(condOrder(cond), :), '.-', 'Color', roiColors(roi,:), 'LineWidth', 2, 'MarkerSize', 20);
            box off;
            ylim(ybounds);
            title(ConditionNames{condOrder(cond)}, 'Interpreter','none');
        end
        figure(roiHandle)
        sgtitle(sprintf('%s: %s avgFIR results, nVoxel: %i', subject, ROInames{roi, 2}, sum(selectedVox)));
        
        subjFigDir = fullfile(figDir, subject);
        if ~exist(subjFigDir, 'dir'), mkdir(subjFigDir); end
        print(roiHandle, fullfile(subjFigDir, sprintf('compareResponses_%s_%s', subject, ROInames{roi,2})), '-dpdf')
        close(roiHandle);

    end

    %-- figure across all ROIs
    figure(figHandle)
    sgtitle(sprintf('%s: avgFIR results', subject));
    subplot(2,7,1)
    legend(ROInames(:, 2), 'FontSize', 16, 'Location', 'northwest');
    legend boxoff
    print(figHandle, fullfile(subjFigDir, sprintf('compareResponses_%s_all', subject)), '-dpdf')

end


%% Group figure time course
groupFigDir = fullfile(figDir, 'groupAvg');
if ~exist(groupFigDir, 'dir'), mkdir(groupFigDir); end

figHandle       = figure('Color', [1 1 1], 'Position', [30 300 1400 800]);
set(figHandle,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[1400 800])
colororder(roiColors);

globalMin       = Inf;
globalMax       = -Inf;

for roi = 1:numROIs
    ybounds     = bounds(allFirBetas{roi});
    globalMin   = min(globalMin, ybounds(1));
    globalMax   = max(globalMax, ybounds(2));
end
globalBounds    = [globalMin, globalMax];

for roi = 1:numROIs
    
    %-- roi figure
    roiHandle       = figure('Color', [1 1 1], 'Position', [30 300 1400 800]);
    set(roiHandle,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[1400 800])

    c       = 0;
    for cond = 1:numel(condOrder)
        
        c               = c + 1;
        %-- 95% confidence interval
        tmpFIR          = squeeze(allFirBetas{roi}(condOrder(cond),:,:))';
        cBoundsFIR      = bootci(100, @mean, tmpFIR); %bootci samples from rows

        figure(figHandle)
        subplot(2, 7, c);
        hold on;
        h = fill([x_data, fliplr(x_data)], [cBoundsFIR(1, :), fliplr(cBoundsFIR(2, :))], roiColors(roi, :), 'HandleVisibility', 'off');
        h.FaceAlpha = 0.2;
        plot(x_data, mean(tmpFIR, 1), '.-', 'Color', roiColors(roi, :), 'LineWidth', 2, 'MarkerSize', 20);

        ylim(globalBounds);
        box off;
        title(ConditionNames{condOrder(cond)}, 'Interpreter','none');

        figure(roiHandle)
        subplot(2, 7, c);
        hold on;
        h = fill([x_data, fliplr(x_data)], [cBoundsFIR(1, :), fliplr(cBoundsFIR(2, :))], roiColors(roi, :), 'HandleVisibility', 'off');
        h.FaceAlpha = 0.2;
        plot(x_data, mean(tmpFIR, 1), '.-', 'Color', roiColors(roi, :), 'LineWidth', 2, 'MarkerSize', 20);

        ylim(bounds(allFirBetas{roi}));
        box off;
        title(ConditionNames{condOrder(cond)}, 'Interpreter','none');

    end

    figure(roiHandle)
    sgtitle(sprintf('%s groupAvg', ROInames{roi, 2}))
    print(roiHandle, fullfile(groupFigDir, sprintf('%s_ROIResponses_%s', ROInames{roi,2}, 'groupAvg')), '-dpdf')
    close(roiHandle);

end

figure(figHandle)
legend(ROInames(:, 2), 'box', 'off', 'FontSize', 16);
sgtitle('Compare ROIs avgFIR results', 'FontSize', 20);
print(figHandle, fullfile(groupFigDir, sprintf('allROIResponses_%s', 'groupAvg')), '-dpdf')

%% Group figure summary metric
figHandle       = figure('Color', [1 1 1], 'Position', [30 300 700 1000]);
set(figHandle,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[700 1000])

onePulseIndx    = cat(1, find(contains(ConditionNames, 'BLANK')), find(contains(ConditionNames, 'ONE')));
twoPulseIndx    = cat(1, find(contains(ConditionNames, 'ONE-PULSE-4')), find(contains(ConditionNames, 'TWO')));

xvalues         = log10(cat(2, stimdur(2)/2, stimdur(2:end)));
xbounds         = [xvalues(1) xvalues(end)] + diff(xvalues(1:2)) * [-.5 .5]; % [-0.1 1.3] or [0 8]
oldbounds       = [0 0];
globalMin       = Inf;
globalMax       = -Inf;

% Visualize summed responses
for roi = 1:numROIs

    %-- reorder data so that one pulse and two pulse conditions are grouped 
    meanSum     = sumBetas{roi}(condOrder, :);
    ci          = bootci(100, @mean, meanSum');

    figure(figHandle)
    subplot(1,2,1)
    hold on,

    if size(sumBetas{roi},2) > 1
        errorbar(xvalues, mean(meanSum(1:numel(stimdur),:),2),  std(meanSum(1:numel(stimdur),:), [], 2)/sqrt(numSubjects), ...
            '.', 'Color', roiColors(roi,:), 'CapSize', 0, 'LineWidth', 3, 'MarkerSize', 50)

        % plot 95% bootstrapped CI
        plot(xvalues .* [1; 1],  ci(:,1:numel(stimdur)), 'Color', roiColors(roi,:), 'LineWidth', 2, 'HandleVisibility', 'off')
        plot(xvalues,  mean(meanSum(1:numel(stimdur), :),2), '.', 'Color', roiColors(roi,:), 'MarkerSize', 50,'HandleVisibility', 'off')
    else
        plot(xvalues, meanSum(1:numel(stimdur), :), ...
            '.', 'Color', roiColors(roi,:), 'LineWidth', 3, 'MarkerSize', 50)
    end
    set(gca, 'XTick', xvalues, 'XTickLabel', num2str(stimdur'))
    ybounds = get(gca, 'YLim');

    globalMin = min(globalMin, ybounds(1));
    globalMax = max(globalMax, ybounds(2));

    subplot(1,2,2)
    hold on,
    if size(sumBetas{roi},2) > 1
        errorbar(xvalues, mean(meanSum((1:numel(stimdur))+numel(stimdur),:),2),  std(meanSum((1:numel(stimdur))+numel(stimdur),:), [], 2)/sqrt(numSubjects), ...
            '.', 'Color', roiColors(roi,:), 'CapSize', 0, 'LineWidth', 3, 'MarkerSize', 50)
        % plot 95% bootstrapped CI
        plot(xvalues .* [1; 1],  ci(:,(1:numel(stimdur))+numel(stimdur)), 'Color', roiColors(roi,:), 'LineWidth', 2, 'HandleVisibility', 'off')
        plot(xvalues,  mean(meanSum ((1:numel(stimdur))+numel(stimdur), :),2), '.', 'Color', roiColors(roi,:), 'MarkerSize', 50,'HandleVisibility', 'off')
    else
        plot(xvalues, sumBetas{roi}(twoPulseIndx,:), ...
            '.', 'Color', roiColors(roi,:), 'LineWidth', 3, 'MarkerSize', 50)

    end


end
figure(figHandle)
subplot(1,2,1)
ylim(ybounds); box off; xlabel('Duration'); xlim(xbounds)
title(sprintf('one pulse conditions'))

subplot(1,2,2)
set(gca, 'XTick', xvalues, 'XTickLabel', num2str(stimdur'))
ylim(ybounds); box off; xlabel('Interstimulus interval'); xlim(xbounds)
title(sprintf('two pulse conditions'))
legend(ROInames(:,2), 'FontSize', 16,'Location','best')

sgtitle(sprintf('Compare ROIs summed response'), 'fontsize', 20)
print(figHandle, fullfile(groupFigDir, sprintf('compareROISumResp_%s', 'groupAvg')), '-dpdf')

