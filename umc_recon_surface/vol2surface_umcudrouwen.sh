
export SUBJID=umcudrouwen
export SESSION=umcu3t01
export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export SUBJ_DIR=${WORK_DIR}/derivatives/freesurfer
export LABEL_DIR=${SUBJ_DIR}/sub-${SUBJID}/label

# Transform Glasser2016 atlas from fsaverage to individual subject
mri_surf2surf --srcsubject fsaverage --trgsubject sub-${SUBJID} \
  --hemi rh \
  --sval-annot ${SUBJ_DIR}/fsaverage/label/rh.Glasser2016.annot \
  --tval ${SUBJ_DIR}/sub-${SUBJID}/surf/rh.Glasser2016.annot

# Convert annotation to volumetric atlas
mri_aparc2aseg \
  --s sub-${SUBJID} \
  --annot Glasser2016 \
  --hemi rh \
  --o ${SUBJ_DIR}/sub-${SUBJID}/surf/rh.glasser16_atlas.mgz

echo "Done"