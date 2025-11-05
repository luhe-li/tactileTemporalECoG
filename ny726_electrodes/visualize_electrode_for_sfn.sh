SUBJID=som726_warped
BIDS_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
SUBJECTS_DIR=${BIDS_DIR}/derivatives/freesurfer/${SUBJID}
SURF_DIR=${SUBJECTS_DIR}/surf
ELEC_DIR=${BIDS_DIR}/derivatives/NY726_2_elec 
LABEL_DIR=${SUBJECTS_DIR}/label

# Define files
ELEC_FILE=${ELEC_DIR}/NY726_2_726_2_elec_bin_T1_2025-04-29.nii.gz
NY726_ELEC_LUT=${ELEC_DIR}/NY726_elec_lut.txt

freeview -f ${SURF_DIR}/rh.white \
    -f ${LABEL_DIR}/rh.Glasser2016.annot \
    -v ${ELEC_FILE}:colormap=${NY726_ELEC_LUT}:opacity=1.0

