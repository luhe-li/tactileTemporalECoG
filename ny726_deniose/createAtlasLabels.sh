#!/bin/bash -l
# JK & IB 28/10/2022
# example to run sh createAtlasLabels.sh wlsubj121 /Volumes/server/Projects/Retinotopy_NYU_3T

# Adapt to tactile Ecog Data

export SUBJID=${1}
# export WORK_DIR=${2} # e.g /CBI/Users/jankurzawski/data/Retinotopy_NYU_3T
export WORK_DIR=/Volumes/server/Projects/BAIR/Data/BIDS/tactile
export SUBJECTS_DIR=${WORK_DIR}/derivatives/freesurfer
export LABEL_DIR=${SUBJECTS_DIR}/sub-${SUBJID}/label
export ROI_SAVE_DIR=${WORK_DIR}/derivatives/roiVols/sub-${SUBJID}

git clone https://github.com/WinawerLab/atlasmgz.git $SUBJECTS_DIR/fsaverage/atlasmgz
git pull $SUBJECTS_DIR/fsaverage/atlasmgz
export DO_IMPORT_GLASSER=0 # Only need to run once, to bring atlas into fsaverage space
export DO_IMPORT_NATIVESPACE=1
export DO_CONVERT_ATLAS=1

# Load glasser2016.mgz and a colorLUT to extract labels in fsaverage space. Convert all labels into an annotation file
if [ "$DO_IMPORT_GLASSER" == 1 ]; then

mkdir -p ${SUBJECTS_DIR}/fsaverage/label/Glasser2016

for hemi in lh rh; do 
	unset labelString

	for i in {1..180}; do  
		mri_cor2label --i ${SUBJECTS_DIR}/fsaverage/atlasmgz/${hemi}.glasser16_atlas.v1_0.mgz --id ${i} \
		--l ${SUBJECTS_DIR}/fsaverage/label/Glasser2016/${hemi}.Glasser2016.${i}.label \
		--surf fsaverage $hemi inflated;

  		fileString="--l ${SUBJECTS_DIR}/fsaverage/label/Glasser2016/${hemi}.Glasser2016.${i}.label "
		labelString="${labelString}${fileString}"
	done; 

	mris_label2annot --s fsaverage --h ${hemi} --ctab ${SUBJECTS_DIR}/fsaverage/atlasmgz/Glasser2016_ColorLUT.txt --a Glasser2016 \
	${labelString}
done

fi

if [ "$DO_IMPORT_NATIVESPACE" == 1 ]; then

	# Transform annotation file with Glasser2016 atlas ROIs into labels in indv subject space
	mri_surf2surf --srcsubject fsaverage --trgsubject sub-$SUBJID --hemi rh --sval-annot $SUBJECTS_DIR/fsaverage/label/rh.Glasser2016 --tval $SUBJECTS_DIR/sub-${SUBJID}/label/rh.Glasser2016.annot
	mri_surf2surf --srcsubject fsaverage --trgsubject sub-$SUBJID --hemi lh --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.Glasser2016 --tval $SUBJECTS_DIR/sub-${SUBJID}/label/lh.Glasser2016.annot

	mkdir -p $SUBJECTS_DIR/sub-${SUBJID}/label/Glasser2016

	mri_annotation2label --subject sub-$SUBJID --hemi rh --annotation Glasser2016 --outdir ${SUBJECTS_DIR}/sub-${SUBJID}/label/Glasser2016
	mri_annotation2label --subject sub-$SUBJID --hemi lh --annotation Glasser2016 --outdir ${SUBJECTS_DIR}/sub-${SUBJID}/label/Glasser2016

fi

if [ "$DO_CONVERT_ATLAS" == 1 ]; then

	mkdir -p ${ROI_SAVE_DIR}

	for hemi in lh rh
	do

		# Transform label/annotation into the volume. Template defines the functional data resolution of the 3D volume 
		# Here is used the first functional run boldref scan. Note: can not be a 4D volume (i.e. timeseries)
		mri_label2vol --o ${ROI_SAVE_DIR}/${hemi}.Glasser2016.VOL.nii.gz \
		--annot ${LABEL_DIR}/${hemi}.Glasser2016.annot \
		--temp ${SUBJECTS_DIR}/sub-${SUBJID}/mri/T1.mgz \
		--fillthresh 0.5 --proj frac 0 1 .1 --subject sub-${SUBJID} --hemi $hemi \
		--regheader ${SUBJECTS_DIR}/sub-${SUBJID}/mri/T1.mgz

	done

fi