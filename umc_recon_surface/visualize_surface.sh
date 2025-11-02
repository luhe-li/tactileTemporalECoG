export SUBJID=umcudrouwen
export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export SUBJECTS_DIR=${WORK_DIR}/derivatives/freesurfer

freeview -v ${SUBJECTS_DIR}/sub-${SUBJID}/mri/T1.mgz \
             ${SUBJECTS_DIR}/sub-${SUBJID}/mri/brainmask.mgz \
          -f ${SUBJECTS_DIR}/sub-${SUBJID}/surf/lh.white:edgecolor=blue \
             ${SUBJECTS_DIR}/sub-${SUBJID}/surf/lh.pial:edgecolor=red
