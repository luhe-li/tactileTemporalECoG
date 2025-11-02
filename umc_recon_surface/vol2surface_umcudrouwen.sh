export SUBJID=umcudrouwen
export SESSION=umcu3t01
export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export T1_DIR=${WORK_DIR}/sub-${SUBJID}/ses-${SESSION}/anat
export SUBJECTS_DIR=${WORK_DIR}/derivatives/freesurfer

recon-all -i ${T1_DIR}/sub-${SUBJID}_ses-${SESSION}_run-1_T1w.nii.gz -s ${SUBJID} -cw256 -all