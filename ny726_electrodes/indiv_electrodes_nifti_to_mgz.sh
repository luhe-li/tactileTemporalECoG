#!/bin/bash

export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export SUBJECTS_DIR=${WORK_DIR}/derivatives/freesurfer/som726_warped
T1=${SUBJECTS_DIR}/mri/T1.mgz
LABELS=${SUBJECTS_DIR}/mri/aparc+aseg.mgz
ELEC_DIR=${WORK_DIR}/derivatives/NY726_2_elec/singleLeads

# Loop through all electrode nifti files 
for ELEC_NII in ${ELEC_DIR}/*_elec_bin_T1_*.nii.gz; do
    base=$(basename ${ELEC_NII} .nii.gz)

    # Convert electrode nifti to mgz, aligned with T1 grid
    mri_convert ${ELEC_NII} \
        ${SUBJECTS_DIR}/mri/${base}.mgz \
        --like ${T1}

    # Create a binary mask of right postcentral gyrus
    #     Desikan–Killiany label IDs:
    #       1028 = lh.precentral, 2028 = rh.precentral
    #       1030 = lh.postcentral, 2030 = rh.postcentral
    if [ ! -f ${SUBJECTS_DIR}/mri/rh_postcentral_mask.mgz ]; then
        mri_binarize --i ${LABELS} \
        --match 2030 \
        --o ${SUBJECTS_DIR}/mri/rh_postcentral_mask.mgz
    else
        echo "Right postcentral mask already exists"
    fi


    # # Intersect electrode mask with the right postcentral mask
    # mri_and \
    #     ${SUBJECTS_DIR}/mri/${base}.mgz \
    #     ${SUBJECTS_DIR}/mri/rh_postcentral_mask.mgz \
    #     ${SUBJECTS_DIR}/mri/${base}_in_rh_postcentral.mgz

    # # Get label summary (which atlas labels overlap with electrodes)
    # mri_segstats \
    #     --seg ${LABELS} \
    #     --mask ${SUBJECTS_DIR}/mri/${base}.mgz \
    #     --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \
    #     --sum ${SUBJECTS_DIR}/stats/${base}_in_labels.stats
done
