clear

% tbUse ECoG_utils;

bidsRootPath = tt_bidsRootPath;
subject = 'umcudrouwen';
session = 'umcuiemu01';

atlasName = 'glasser16_atlas';
specs = [];
specs.plotelecs     = 'yes';
specs.plotlabel     = 'no';

%% Visualize electrodes: not working because electrode tsv file is empty:
% check sub-umcudrouwen_ses-umcuiemu01_task-vtstemporalpattern_acq-micromed_run-1_electrodes

% % lateral view
% bidsEcogPlotElectrodesOnMesh(bidsRootPath, subject, [], atlasName, [], specs);
% view(-60,-10);
% % 
% % medial view
% bidsEcogPlotElectrodesOnMesh(bidsRootPath, subject, [], atlasName, [], specs);
% view(50,0);
% 
% % ventral view
% bidsEcogPlotElectrodesOnMesh(bidsRootPath, subject, [], atlasName, [], specs);
% view(180,-90);

%% Plot atlas only

fsDir = fullfile(bidsRootPath, 'derivatives', 'freesurfer',  sprintf('sub-%s', subject));
atlas_file_rh = fullfile(fsDir, 'surf', sprintf('rh.%s.mgz', atlasName));

atlas_r = load_mgh(atlas_file_rh);

% ROIs
% BA3b: 9
% BA1: 51
% BA2: 52
atlas_labels = unique(atlas_r(:));
roi_inds = [9 51 52];

clut = repmat(0.7, [numel(atlas_labels),3]);
alpha = repmat(0, [numel(atlas_labels),1]);
roi_color = parula(numel(roi_inds));
for ii = 1:numel(roi_inds)
    i_roi = roi_inds(ii);
    clut(i_roi,:) = roi_color(ii,:);
    alpha(i_roi,:) = 1;
end

% v = volshow(atlas_r,Colormap=clut,Alphamap=alpha);
v = volshow(atlas_r,Colormap=clut);