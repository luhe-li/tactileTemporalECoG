% This script visualizes the overlay between the electrodes and postcentral
% gyrus defined in Desikan-Killiany atlas 

% Technically no conversion needs to be done before running this script,
% but I converted nifti to mgz by all_electrodes_nifti_to_mgz.sh just to
% keep consistency between ROI and electrode formats.

clear; close all; clc

addpath('/Applications/freesurfer/8.1.0/matlab');

% Load ROI
subj_dir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/freesurfer/som726_warped/mri';
atlas_file = fullfile(subj_dir, 'aparc+aseg.mgz');
atlas = MRIread(atlas_file);
atlas_vol = atlas.vol;

% Load electrodes
all_electrodes_file = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/freesurfer/som726_warped/mri/NY726_2_elec.mgz';
mri = MRIread(all_electrodes_file);
elec_vol = mri.vol;
elec_mask = elec_vol ~= 0;

% Check the full atlas brain
volshow(atlas_vol);

%% Overlay electrodes with atlas, ROI in red

% 2022 ctx-rh-postcentral
% 2024 ctx-rh-precentral
roi_idx = 2022;

% Initiate
combined_vol = zeros(size(atlas_vol));

% Label the full atlas
combined_vol(atlas_vol~=0) = 1;

% Label ROI
combined_vol(atlas_vol == roi_idx) = 2;

% Label electrodes
combined_vol(elec_vol > 0) = 3;

cmap = [
    0 0 0;  
    1 1 1;
    0.9 0.3 0.3;
    0 1 0         
];

alpha = [0; 0.01; 0.5; 1];

v = volshow(combined_vol, ...
    'Colormap', cmap, ...
    'Alphamap', alpha);

viewer = v.Parent;
viewer.BackgroundColor = [0 0 0];
viewer.BackgroundGradient = 'off';

%% Find how many voxel overzlap with the roi for each contact

roi_mask = atlas_vol == roi_idx;
elec_mask = elec_vol > 0;

contacts = unique(elec_vol(:));
contacts(contacts==0) = [];

overlapped_contacts = [];
overlap_voxels_per_contact = [];
for c = contacts'
    mask = elec_vol==c;
    overlap_vox = nnz(mask & roi_mask);
    if overlap_vox > 0
        overlapped_contacts = [overlapped_contacts, c];
        overlap_voxels_per_contact = [overlap_voxels_per_contact, overlap_vox];
    end
end

%% Organize selected electrodes information

elec_folder = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/NY726_2_elec';
elec_txt_file = fullfile(elec_folder, 'NY726_2_elec.txt');
elec_tbl = readtable(elec_txt_file, 'Delimiter', ',', 'NumHeaderLines', 1, 'FileType', 'text');
contact_labels = elec_tbl{:,1};
contact_names = elec_tbl{:,2};

num_contacts = numel(overlapped_contacts);
contact_table_data = cell(num_contacts, 4);

for i = 1:num_contacts
    idx = find(contact_labels == overlapped_contacts(i), 1);
    if ~isempty(idx)
        contact_name = contact_names{idx};
    else
        contact_name = '';
    end
    contact_table_data{i,1} = overlapped_contacts(i);
    contact_table_data{i,2} = contact_name;
    contact_table_data{i,3} = overlap_voxels_per_contact(i);
    contact_table_data{i,4} = 'postcentral (idx: 2022)';
end

T = cell2table(contact_table_data, ...
    'VariableNames', {'ContactIndex','ContactName','NumVoxels','ROI'});

T

output_csv_file = fullfile(elec_folder, 'electrodes_in_postcentral.csv');
writetable(T, output_csv_file);
