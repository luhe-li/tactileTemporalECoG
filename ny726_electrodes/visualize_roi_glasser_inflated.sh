#!/bin/bash

export SUBJID=som726_warped
export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export SUBJECTS_DIR=${WORK_DIR}/derivatives/freesurfer/${SUBJID}
export LABEL_DIR=${SUBJECTS_DIR}/label

freeview -f \
  $SUBJECTS_DIR/surf/rh.inflated:annot=$LABEL_DIR/rh.Glasser2016.annot

# Visualize volume
# freeview -v ${SUBJECTS_DIR}/surf/T1.mgz \
# -f ${SUBJECTS_DIR}/surf/rh.white:annot=${LABEL_DIR}/rh.Glasser2016.annot \
# -f ${SUBJECTS_DIR}/surf/rh.pial:annot=${LABEL_DIR}/rh.Glasser2016.annot &

