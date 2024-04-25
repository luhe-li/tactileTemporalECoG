
clear; close all

% manage path
tt_RootPath;

%% 1. Load ECoG best-fitting parameters% Load data and fits

modelfun = @DN;
xvalmode = 0;
datatype = 'electrodeaverages';
% datatype = 'individualelecs';
[D] = tt_loadDataForFigure(modelfun, xvalmode, datatype);

