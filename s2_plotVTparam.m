clear; close all;

% Compare parameter estimates from the same model and fitting procedure,
% without cross-validation

% Dependency:
% - make sure fitting results files are in the folder: 
%    - tactile: /Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/modelFit/results/
%    - visual: temporalECoG/analysis/results
% OR, to rerun model fitting
% - run s1_fitECoG.m to obtain fits of DN to both 
%   individualelecs and electrodeaverages
% - run tde_run with DN

modelfun = @DN;
xvalmode = 0;

%% 1. Load tactile averaged electrode fitting results

% Load data and fits
datatype = 'electrodeaverages';
[D{1}] = tt_loadDataForFigure(modelfun, xvalmode, datatype);

%% 2. Load individual electrodes fitting results

% Load data and fits
datatype = 'individualelecs';
[D{2}] = tt_loadDataForFigure(modelfun, xvalmode, datatype);

%% 3. Load visual electrode fits (within each visual area) fitting results

% Load data and fits
datatype = 'individualelecs';
[D{3}] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

%% 4. Plot

saveDir = fullfile(tt_bidsRootPath, 'derivatives', 'modelFit', 'figure', 'compareTactVis');
tt_plotTactileVisualParams(D, modelfun, saveDir)