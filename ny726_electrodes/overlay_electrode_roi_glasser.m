
clear; close all; clc

addpath('/Applications/freesurfer/8.1.0/matlab');

% Load ROI
subj_id = 'som726_warped';
bids_dir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
ROI_dir = fullfile(bids_dir, 'derivatives', 'roiVols', subj_id);

glasser_file = fullfile(ROI_dir, 'rh.Glasser2016_RAS.VOL.nii.gz');
glasser = MRIread(glasser_file);
glasser_vol = glasser.vol;

% Check the full atlas brain
% volshow(glasser_vol);

% Load electrodes in T1
all_electrodes_file = fullfile(bids_dir, 'derivatives', 'NY726_2_elec', 'NY726_2_726_2_elec_bin_T1_2025-04-29.nii.gz');
elec_file = MRIread(all_electrodes_file);
elec_vol = elec_file.vol;
elec_mask = elec_vol ~= 0;

%% Check and convert file coordinates (done, just for the record)

% mri_info /Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/roiVols/som726_warped/rh.Glasser2016.VOL.nii.gz | grep Orientation
% mri_info /Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/NY726_2_elec/NY726_2_726_2_elec_bin_T1_2025-04-29.nii.gz | grep Orientation

% Found that glasser is in LIA coordiante and electrodes are in RAS
% coordinates, so they are flipped

% Convert glasser LIA to the more common RAS
% mri_convert /Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/roiVols/som726_warped/rh.Glasser2016.VOL.nii.gz /Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/roiVols/som726_warped/rh.Glasser2016_RAS.VOL.nii.gz --out_orientation RAS


%% Overlay electrodes with Glasser ROI

% ROI = S1 (BA3a, BA3b, BA1, BA2): 53, 9, 51, 52 (Glasser indices)
glasser_labels = unique(glasser_vol(:));
roi_inds = [53 9 51 52];

% Initiate
combined_vol = zeros(size(glasser_vol));

% Label the full atlas
combined_vol(glasser_vol~=0) = 1;

% Label the ROIs
for rr = 1:numel(roi_inds)
    combined_vol(glasser_vol == roi_inds(rr)) = 1+rr;
end

% Label the electrodes
combined_vol(elec_mask) = numel(roi_inds) + 2;

lut = [
    0     0     0;        % background
    1     1     1;        % whole atlas (white)
    0     0     205/255;  % BA3a 
    0     191/255 1;      % BA3b 
    1     0     1;        % BA1 
    139/255 0 139/255;    % BA2 
    0     1     0         % electrodes
];
roi_alpha = repmat(0.5, [numel(roi_inds), 1]);
alpha = [0; 0.005; roi_alpha; 1]; 

v = volshow(combined_vol, ...
    'Colormap', lut, ...
    'Alphamap', alpha);

viewer = v.Parent;
viewer.BackgroundColor = [0 0 0];
viewer.BackgroundGradient = 'off';

%% Find how many voxel overlap with the ROIs for each contact

% create a mask for S1 by combining the masks for each ROI
roi_mask = zeros(size(glasser_vol));
for rr = 1:numel(roi_inds)
    rr_roi_mask = glasser_vol == roi_inds(rr);
    roi_mask = roi_mask | rr_roi_mask;
end

% Unique contacts index
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

elec_folder = fullfile(bids_dir, 'derivatives', 'NY726_2_elec');
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
    contact_table_data{i,4} = 'S1';
end

T = cell2table(contact_table_data, ...
    'VariableNames', {'ContactIndex','ContactName','NumVoxels','ROI'});

T


output_csv_file = fullfile(elec_folder, 'electrodes_in_S1.csv');
writetable(T, output_csv_file);
