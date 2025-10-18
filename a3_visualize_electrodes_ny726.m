
clear; close all;
tbUse tactileTemporalECoG;
addpath(fullfile(pwd, 'codes'));

atlasName           =  {'wang15_fplbl_norm'};
projectDir          = tt_bidsRootPath; 

specs = [];
specs.plotelecs     = 'yes';
specs.plotlabel     = 'no';

subject             = 'ny726'; 
specs.plotmesh      = 'right';
specs.plotelecrad   = 2.3;
specs.plotlabel = 'yes';
% lateral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(-60,-10);

% medial view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(50,0);

% ventral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(180,-90);