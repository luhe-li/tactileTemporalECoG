#!/bin/bash

SUBJID=som726_warped
BIDS_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
SUBJECTS_DIR=${BIDS_DIR}/derivatives/freesurfer/${SUBJID}
ROI_VOL=${BIDS_DIR}/derivatives/roiVols/${SUBJID}/rh.Glasser2016_RAS.VOL.nii.gz
LABEL_DIR=${SUBJECTS_DIR}/label/Glasser2016
ELEC_FILE=${BIDS_DIR}/derivatives/NY726_2_elec/NY726_2_726_2_elec_bin_T1_2025-04-29.nii.gz

freeview \
  -v \
    ${SUBJECTS_DIR}/mri/T1.mgz:grayscale=0,1200 \
    ${ROI_VOL}:colormap=lut:opacity=0.4 \
    ${ELEC_FILE}:color=0,1,0:opacity=1.0 \
  -label ${LABEL_DIR}/rh.1.label:color=1,0,1:opacity=0.5 \
  -label ${LABEL_DIR}/rh.2.label:color=0.55,0,0.55:opacity=0.5 \
  -label ${LABEL_DIR}/rh.3a.label:color=0,0,0.8:opacity=0.5 \
  -label ${LABEL_DIR}/rh.3b.label:color=0,0.75,1:opacity=0.5
