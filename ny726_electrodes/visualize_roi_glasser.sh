#!/bin/bash

export SUBJID=som726_warped
export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export SUBJECTS_DIR=${WORK_DIR}/derivatives/freesurfer
export LABEL_DIR=${SUBJECTS_DIR}/${SUBJID}/label

freeview -v ${SUBJECTS_DIR}/sub-${SUBJID}/mri/T1.mgz \
-f ${SUBJECTS_DIR}/sub-${SUBJID}/surf/lh.white:annot=${SUBJECTS_DIR}/sub-${SUBJID}/label/lh.Glasser2016.annot \
-f ${SUBJECTS_DIR}/sub-${SUBJID}/surf/rh.white:annot=${SUBJECTS_DIR}/sub-${SUBJID}/label/rh.Glasser2016.annot \
-f ${SUBJECTS_DIR}/sub-${SUBJID}/surf/lh.pial:annot=${SUBJECTS_DIR}/sub-${SUBJID}/label/lh.Glasser2016.annot \
-f ${SUBJECTS_DIR}/sub-${SUBJID}/surf/rh.pial:annot=${SUBJECTS_DIR}/sub-${SUBJID}/label/rh.Glasser2016.annot &



