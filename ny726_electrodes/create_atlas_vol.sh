#!/bin/bash

export SUBJID=som726_warped
export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export SUBJECTS_DIR=${WORK_DIR}/derivatives/freesurfer
export LABEL_DIR=${SUBJECTS_DIR}/${SUBJID}/label
export ROI_SAVE_DIR=${WORK_DIR}/derivatives/roiVols/${SUBJID}

export DO_IMPORT_NATIVESPACE=1

if [ "$DO_IMPORT_NATIVESPACE" -eq 1 ]; then

	# Transform annotation file with Glasser2016 atlas ROIs into labels in indv subject space
	mri_surf2surf --srcsubject fsaverage --trgsubject $SUBJID --hemi rh --sval-annot $SUBJECTS_DIR/fsaverage/label/rh.Glasser2016 --tval $SUBJECTS_DIR/${SUBJID}/label/rh.Glasser2016.annot
	mri_surf2surf --srcsubject fsaverage --trgsubject $SUBJID --hemi lh --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.Glasser2016 --tval $SUBJECTS_DIR/${SUBJID}/label/lh.Glasser2016.annot

	mkdir -p $SUBJECTS_DIR/${SUBJID}/label/Glasser2016

	mri_annotation2label --subject $SUBJID --hemi rh --annotation Glasser2016 --outdir ${SUBJECTS_DIR}/${SUBJID}/label/Glasser2016
	mri_annotation2label --subject $SUBJID --hemi lh --annotation Glasser2016 --outdir ${SUBJECTS_DIR}/${SUBJID}/label/Glasser2016
fi

# # Convert annotation to volume (per hemisphere)
# mkdir -p ${ROI_SAVE_DIR}
# for hemi in lh rh; do
#     mri_label2vol \
#       --o ${ROI_SAVE_DIR}/${hemi}.Glasser2016.VOL.nii.gz \
#       --annot ${LABEL_DIR}/${hemi}.Glasser2016.annot \
#       --temp ${SUBJECTS_DIR}/${SUBJID}/mri/T1.mgz \
#       --fillthresh 0.5 \
#       --proj frac 0 1 0.1 \
#       --subject ${SUBJID} \
#       --hemi ${hemi} \
#       --regheader ${SUBJECTS_DIR}/${SUBJID}/mri/T1.mgz
# done

# # Convert NIfTI to MGZ, required for ECOGutil electrode visualization functions
# for hemi in lh rh; do
#     mri_convert \
#       ${ROI_SAVE_DIR}/${hemi}.Glasser2016.VOL.nii.gz \
#       ${SUBJECTS_DIR}/${SUBJID}/surf/${hemi}.glasser16_atlas.mgz
# done

echo "Done."
