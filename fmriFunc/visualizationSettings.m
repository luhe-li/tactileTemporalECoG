
function out = visualizationSettings(models, currModel)

out.condNames   = {'ONE-PULSE-0'; % blank
                   'ONE-PULSE-1';
                   'ONE-PULSE-2';
                   'ONE-PULSE-3';
                   'ONE-PULSE-4';
                   'ONE-PULSE-5';
                   'ONE-PULSE-6';
                   'TWO-PULSE-0'; % ONE_PULSE-4
                   'TWO-PULSE-1';
                   'TWO-PULSE-2';
                   'TWO-PULSE-3';
                   'TWO-PULSE-4';
                   'TWO-PULSE-5';
                   'TWO-PULSE-6';
                   };
out.numCond     = numel(out.condNames); 
    
out.stimDur     = [0 0.05 0.1 0.2 0.4 0.8 1.2]; % seconds
out.twoPulseDur = 0.2; % pulse duration in two-pulse condition, in seconds

out.x_data      = 0:1:10; 
out.xvalues     = log10(cat(2, out.stimDur(2)/2, out.stimDur(2:end)));
out.xbounds     = [out.xvalues(1) out.xvalues(end)] + diff(out.xvalues(1:2)) * [-.5 .5]; % [-0.1 1.3] or [0 8]

if nargin > 0

    modelNames      = {'DN', 'NORM', 'HRF'};
    wModel          = ismember( modelNames, currModel);
    assert(sum(wModel) > 0, sprintf('[%s] %s Model name does not match visualization settings', mfilename, currModel))

    % Colors used for the models
    modelColors     = parula(numel(modelNames)+1);
    [~,colorOrder]  = ismember(models, modelNames);

    % Parameter bounds for visualization
    modelStr        = str2func(sprintf('fit%smodel', currModel));
    modelSpecs      = modelStr([], {'initialize'});
    NORMmodel       = fitNORMmodel([], {'initialize'}); % model with most params

    modelBounds     = cat(1, NORMmodel.lb, NORMmodel.ub); % use model constrainst as axis bounds
    modelBounds(isinf(modelBounds)) = [8 15]; % upper axis bounds overall gain and negative gamma gain
    modelBounds(2,strcmp(NORMmodel.labels, 'sigma')) = 0.3; % upper axis bound

    % output
    out.color       = modelColors(wModel,:);
    out.mColors     = modelColors(colorOrder,:);
    out.mBounds     = modelBounds;
    out.mLabels     = modelSpecs.labels;
    out.allLabels   = NORMmodel.labels;
    out.labels      = {'tau', 'sigma', 'gain', 'p1', 'p2', 'w'}; % match names in manuscript
    out.totParams   = numel(NORMmodel.labels);
    out.model       = modelStr;

end

