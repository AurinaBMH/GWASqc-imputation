#!/bin/bash

#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-48:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=aurina.arnatkeviciute@monash.edu
#SBATCH --mail-type=ALL

# massive usage (individual subject)
#  1. [code]$ WHERESMYSCRIPT=$(pwd); SUBJECTID=XXXXX
#  2. [....]$ sbatch --job-name="${SUBJECTID}" --output="/projects/kg98/aurina/Freesurfer/MichellesStudy/errorLog/slurm-${SUBJECTID}.out" --error="/projects/kg98/aurina/Freesurfer/MichellesStudy/errorLog/slurm-${SUBJECTID}.err" ${WHERESMYSCRIPT}/Autoreconall_AA.sh $WHERESMYSCRIPT $SUBJECTID

# massive usage (more than one subject)
#  1. [code]$ WHERESMYSCRIPT=$(pwd); SUBJECTIDS="XXXXX XXXXX ... XXXXX"
#  2. [....]$ for SUBJECTID in $SUBJECTIDS; do sbatch --job-name="${SUBJECTID}" --output="//projects/kg98/aurina/Freesurfer/MichellesStudy/errorLog/slurm-${SUBJECTID}.out" --error="/projects/kg98/aurina/Freesurfer/MichellesStudy/errorLog/slurm-${SUBJECTID}.err" ${WHERESMYSCRIPT}/Autoreconall_AA.sh $WHERESMYSCRIPT $SUBJECTID; done

WHERESMYSCRIPT=$1
SUBJECTID=$2
FILESEP="/"
echo -e "\nSETTING UP FREESURFER\n"
source ${WHERESMYSCRIPT}${FILESEP}SetupEnv_A.sh


SUBJECTDIR="${SUBJECTS_DIR}${FILESEP}${SUBJECTID}"
#Parse variables
mkdir ${SUBJECTDIR}/mri
mkdir ${SUBJECTDIR}/mri/orig
cp ${SUBJECTDIR}/t1/t1.nii ${SUBJECTDIR}/mri/orig/toct1.nii.gz
robustfov -i toct1.nii.gz -r 001.nii.gz # remove neck
mri_convert ${SUBJECTDIR}/mri/orig/001.nii.gz ${SUBJECTDIR}/mri/orig/001.mgz
# delete originals
rm  ${SUBJECTDIR}/mri/orig/001.nii.gz
rm  ${SUBJECTDIR}/mri/orig/toct1.nii.gz
WORKDIR="${SUBJECTDIR}${FILESEP}mri${FILESEP}orig"

echo -e "\nRUNNING AUTORECON ALL ON ${SUBJECTID}\n"
recon-all -subjid ${SUBJECTID} -all
