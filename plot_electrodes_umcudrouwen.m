
tbUse ECoG_utils;

bidsRootPath = tt_bidsRootPath;
subject = 'umcudrouwen';
session = 'umcuiemu01';

atlasName = {'glasser16_atlas'};
specs = [];
specs.plotelecs     = 'yes';
specs.plotlabel     = 'no';

%% 

% lateral view
bidsEcogPlotElectrodesOnMesh(bidsRootPath, subject, [], atlasName, [], specs);
view(-60,-10);
% 
% % medial view
% bidsEcogPlotElectrodesOnMesh(bidsRootPath, subject, [], atlasName, [], specs);
% view(50,0);
% 
% % ventral view
% bidsEcogPlotElectrodesOnMesh(bidsRootPath, subject, [], atlasName, [], specs);
% view(180,-90);