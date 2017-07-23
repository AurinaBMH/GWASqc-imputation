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
#  1. [code]$ WHERESMYSCRIPT=$(pwd); CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"; WHEREISDATA=${WHEREISIMPUTATIONQC}/preparationResults
#  2. [....]$ for CHROMOSOME in $CHROMOSOMES; do sbatch --job-name="${CHROMOSOME}" --output="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}.out" --error="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}.err" ${WHERESMYSCRIPT}/split_chromosomes.sh $WHERESMYSCRIPT $CHROMOSOME $WHEREISDATA; done

WHERESMYSCRIPT=$1
CHROMOSOME=$2
WHEREISDATA=$3

WHEREISIMPUTATIONQC="/projects/kg98/aurina/GWAS/data/DondersQC/imputationQC"
cd ${WHEREISIMPUTATIONQC}
echo -e "\nSPLITTING CHROMOSOME NR ${CHROMOSOME}\n"

./plink --bfile ${WHEREISDATA}/superclean_step5 --chr ${CHROMOSOME} --recode --out ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.${CHROMOSOME}
