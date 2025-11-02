
clear; close all; clc

addpath('/Applications/freesurfer/8.1.0/matlab');

% Load electrodes in T1
all_electrodes_file = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/NY726_2_elec/NY726_2_726_2_elec_bin_T1_2025-04-29.nii.gz';
elecs = MRIread(all_electrodes_file);
elecs_vol = elecs.vol;

% Unique contacts index
contacts = unique(elecs_vol(:));
contacts(contacts==0) = []; 

% Load Glasser ROI
roi_dir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/roiVols/som726_warped';
glasser_file = fullfile(roi_dir, 'lh.Glasser2016.VOL.nii.gz');
glasser = MRIread(glasser_file);
glasser_vol = glasser.vol;

% ROIs
% BA3b: 9
% BA1: 51
% BA2: 52
glasser_labels = unique(glasser_vol(:));
roi_inds = [9 51 52];

roi_vol = zeros(size(glasser_vol));
for rr = 1:numel(roi_inds)
    roi_vol = roi_vol + (glasser_vol == roi_inds(rr));
end

% % Find overlap
% overlap = elecs_vol & glasser_vol;
% overlap_add = elecs_vol + glasser_vol;

add_mask = elecs_vol + roi_vol;
volshow(addmask)