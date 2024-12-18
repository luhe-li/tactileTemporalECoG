
clear; close all;

% Compare parameter estimates from the same model and fitting procedure
modelfun = @LINEAR_RECTF_EXP_NORM;
xvalmode = 0;

%% 1. Load averaged electrode fitting results

% Load data and fits
datatype = 'electrodeaverages';
[D{1}] = tt_loadDataForFigure(modelfun, xvalmode, datatype);

%% 2. Load individual electrodes fitting results

% Load data and fits
datatype = 'individualelecs';
[D{2}] = tt_loadDataForFigure(modelfun, xvalmode, datatype);

%% 3. Load visual parameters

% Load data and fits
datatype = 'individualelecs';
[D{3}] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

%% 4. Plot

saveDir = fullfile(tt_bidsRootPath, 'derivatives', 'modelFit', 'figure', 'compareTactVis');

tt_plotTactileVisualParams(D, modelfun, saveDir)