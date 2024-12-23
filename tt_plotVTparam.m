clear; close all;

% Compare parameter estimates from the same model and fitting procedure,
% without cross-validation
modelfun = @DN;
xvalmode = 0;

% Add util folders copied from temporalECoG repository to load visual data
% results
addpath(pwd,'utils');

%% 1. Load tactile averaged electrode fitting results

% Load data and fits
dataPath = fullfile(pwd, 'modelFitsResults','tactile');
dataType = 'electrodeaverages';
D{1} = tt_loadDataForFigure(modelfun, xvalmode, dataType, dataPath);

%% 2. Load individual electrodes fitting results

% Load data and fits
dataPath = fullfile(pwd, 'modelFitsResults','tactile');
dataType = 'individualelecs';
D{2} = tt_loadDataForFigure(modelfun, xvalmode, dataType, dataPath);

%% 3. Load visual electrode fitting results: fits of individual electrodes averaged within each visual area

% Load data and fits
dataPath = fullfile(pwd, 'modelFitsResults','visual');
dataType = 'individualelecs';
D{3} = tt_loadDataForFigure(modelfun, xvalmode, dataType, dataPath);

%% 4. Plot
saveDir = fullfile(pwd, mfilename);
if ~exist("saveDir", 'dir'); mkdir(saveDir); end
tt_plotTactileVisualParams(D, modelfun, saveDir);