#!/bin/bash


#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-96:00:00
#SBATCH --mem-per-cpu=16000
#SBATCH --mail-user=aurina.arnatkeviciute@monash.edu
#SBATCH --mail-type=ALL

# massive usage (more than one subject)
#  1. [code]$ WHERESMYSCRIPT=$(pwd); CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"; WHEREISDATA=${WHEREISIMPUTATIONQC}/phasingData
#  2. [....]$ for CHROMOSOME in $CHROMOSOMES; do sbatch --job-name="${CHROMOSOME}_imputation2run" --output="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}_imputation2run.out" --error="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}_imputation2run.err" ${WHEREISCODE}/do_imputation.sh $CHROMOSOME $WHEREISDATA; done

CHROMOSOME=$1
WHEREISDATA=$2

WHEREIS1000GENOMES="/projects/kg98/aurina/GWAS/data/enigma/1KGPref"
WHEREISIMPUTATIONQC="/projects/kg98/aurina/GWAS/data/DondersQC/imputationQC"


cd ${WHEREISDATA}
	for ((j=1; j<=4; j++)); do
		if test -f chunk${j}-ready4mach.${CHROMOSOME}.dat.gz; then
			./minimac2 --vcfReference --rounds 5 --states 200 --refHaps ${WHEREIS1000GENOMES}/chr${CHROMOSOME}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.gz --haps ${WHEREISDATA}/chunk${j}-ready4mach.${CHROMOSOME}.gz --snps ${WHEREISDATA}/chunk${j}-ready4mach.${CHROMOSOME}.dat.gz.snps --autoClip ${WHEREISDATA}/autochunk-ready4mach.${CHROMOSOME}.dat.gz --gzip --prefix ${WHEREISIMPUTATIONQC}/imputationResults/chunk${j}-ready4mach.${CHROMOSOME}.imputed > ${WHEREISIMPUTATIONQC}/imputationResults/chunk${j}-ready4mach.${CHROMOSOME}-minimac2.log
		fi
	done



#
