export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export SUBJECTS_DIR=${WORK_DIR}/derivatives/freesurfer/som726_warped
T1=${SUBJECTS_DIR}/mri/T1.mgz
LABELS=${SUBJECTS_DIR}/mri/aparc+aseg.mgz
ALL_ELECTRODES_FILE=${WORK_DIR}/derivatives/NY726_2_elec/NY726_2_elec.nii.gz


# # Create a binary mask of the right postcentral gyrus (2030)
# if [ ! -f ${SUBJECTS_DIR}/mri/rh_postcentral_mask.mgz ]; then
#     mri_binarize --i ${LABELS} \
#       --match 2030 \
#       --o ${SUBJECTS_DIR}/mri/rh_postcentral_mask.mgz
# else
#     echo "Right postcentral mask already exists"
# fi

# Convert ALL_ELECTRODES_FILE to .mgz format aligned with T1, if not already done
ALL_ELECTRODES_MGZ=${SUBJECTS_DIR}/mri/NY726_2_elec.mgz
if [ ! -f ${ALL_ELECTRODES_MGZ} ]; then
    mri_convert ${ALL_ELECTRODES_FILE} \
      ${ALL_ELECTRODES_MGZ} \
      --like ${T1}
fi

# Check the header of the original all-electrodes nifti and the .mgz version
mri_vol2vol \
  --mov  ${ALL_ELECTRODES_MGZ} \
  --targ ${LABELS} \
  --regheader \
  --interp nearest \
  --o ${SUBJECTS_DIR}/mri/NY726_2_elec_in_aparc_space.mgz