#!/bin/bash

export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export SUBJ_DIR=${WORK_DIR}/derivatives/freesurfer/som726_warped
T1=${SUBJ_DIR}/mri/T1.mgz
LABELS=${SUBJ_DIR}/mri/aparc+aseg.mgz
ELEC_DIR=${WORK_DIR}/derivatives/NY726_2_elec/singleLeads

OVERLAYS=()
for f in ${SUBJ_DIR}/mri/*_in_rh_postcentral.mgz; do
    OVERLAYS+=("$f:colormap=heat:opacity=0.6")
done

freeview \
  -v ${T1} \
     ${LABELS}:colormap=lut:opacity=0.4 \
     "${OVERLAYS[@]}"
