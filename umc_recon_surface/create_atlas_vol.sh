#!/bin/bash

export SUBJID=umcudrouwen
export SESSION=umcu3t01
export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export SUBJECTS_DIR=${WORK_DIR}/derivatives/freesurfer
export LABEL_DIR=${SUBJECTS_DIR}/sub-${SUBJID}/label
export ROI_SAVE_DIR=${WORK_DIR}/derivatives/roiVols/sub-${SUBJID}

export DO_IMPORT_NATIVESPACE=0

# Map Glasser2016 from fsaverage to subject space (optional)
if [ "$DO_IMPORT_NATIVESPACE" -eq 1 ]; then
    mri_surf2surf --srcsubject fsaverage --trgsubject sub-${SUBJID} \
      --hemi rh \
      --sval-annot ${SUBJECTS_DIR}/fsaverage/label/rh.Glasser2016.annot \
      --tval ${LABEL_DIR}/rh.Glasser2016.annot

    mri_surf2surf --srcsubject fsaverage --trgsubject sub-${SUBJID} \
      --hemi lh \
      --sval-annot ${SUBJECTS_DIR}/fsaverage/label/lh.Glasser2016.annot \
      --tval ${LABEL_DIR}/lh.Glasser2016.annot
fi

# Convert annotation to volume (per hemisphere)
mkdir -p ${ROI_SAVE_DIR}
for hemi in lh rh; do
    mri_label2vol \
      --o ${ROI_SAVE_DIR}/${hemi}.Glasser2016.VOL.nii.gz \
      --annot ${LABEL_DIR}/${hemi}.Glasser2016.annot \
      --temp ${SUBJECTS_DIR}/sub-${SUBJID}/mri/T1.mgz \
      --fillthresh 0.5 \
      --proj frac 0 1 0.1 \
      --subject sub-${SUBJID} \
      --hemi ${hemi} \
      --regheader ${SUBJECTS_DIR}/sub-${SUBJID}/mri/T1.mgz
done

# Convert NIfTI to MGZ
for hemi in lh rh; do
    mri_convert \
      ${ROI_SAVE_DIR}/${hemi}.Glasser2016.VOL.nii.gz \
      ${SUBJECTS_DIR}/sub-${SUBJID}/surf/${hemi}.glasser16_atlas.mgz
done

echo "Done."
