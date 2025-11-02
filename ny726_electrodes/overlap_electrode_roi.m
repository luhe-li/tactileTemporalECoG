clear; close all

addpath('/Applications/freesurfer/8.1.0/matlab');

subj_dir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/freesurfer/som726_warped/mri';
atlas_file = fullfile(subj_dir, 'aparc+aseg.mgz');
atlas = MRIread(atlas_file);
atlas_vol = atlas.vol;

% ROI definitions Desikan–Killiany default
roi_defs = struct( ...
    'rh_precentral', 2028, ...
    'rh_postcentral', 2030);

files = dir(fullfile(subj_dir, '*_elec_bin_T1_*.mgz'));

for f = 1:numel(files)
    fname = fullfile(subj_dir, files(f).name);
    mri = MRIread(fname);
    elec_vol = mri.vol;

    contacts = unique(elec_vol(:));
    contacts(contacts==0) = [];  % remove background

    for c = contacts'
        mask = (elec_vol == c);
        total_vox = nnz(mask);

        % Loop through each ROI
        for roi_name = fieldnames(roi_defs)'
            roi_id = roi_defs.(roi_name{1});
            roi_mask = (atlas_vol == roi_id);

            overlap_vox = nnz(mask & roi_mask);
            if overlap_vox > 0
                fprintf('Electrode %s contact %d overlaps %d voxels with %s\n', ...
                    files(f).name, c, overlap_vox, roi_name{1});
            end
        end
    end
end


%% Check all electrodes file

all_electrodes_file = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile/derivatives/freesurfer/som726_warped/mri/NY726_2_elec.mgz';
mri = MRIread(all_electrodes_file);
elec_vol = mri.vol;

contacts = unique(elec_vol(:));
contacts(contacts==0) = []; 

% Loop through each ROI
for roi_name = fieldnames(roi_defs)'
    roi_id = roi_defs.(roi_name{1});
    roi_mask = (atlas_vol == roi_id);
    
    overlapped_contacts = [];
    overlap_voxels_per_contact = [];
    for c = contacts'
        mask = (elec_vol == c);
        overlap_vox = nnz(mask & roi_mask);
        if overlap_vox > 0
            overlapped_contacts = [overlapped_contacts, c];
            overlap_voxels_per_contact = [overlap_voxels_per_contact, overlap_vox];
        end
    end
    
    if ~isempty(overlapped_contacts)
        fprintf('ROI %s (id %d): %d contacts overlapped:\n', roi_name{1}, roi_id, numel(overlapped_contacts));
        for ii = 1:numel(overlapped_contacts)
            fprintf('    Contact %d: %d voxels\n', overlapped_contacts(ii), overlap_voxels_per_contact(ii));
        end
    else
        fprintf('ROI %s (id %d): No contacts overlapped\n', roi_name{1}, roi_id);
    end
end




%% Check all electrodes file