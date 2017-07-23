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
#  1. [code]$ WHERESMYSCRIPT=$(pwd); CHROMOSOMES="23"; WHEREISDATA=${WHEREISIMPUTATIONQC}/phasingData
#  2. [....]$ for CHROMOSOME in $CHROMOSOMES; do sbatch --job-name="${CHROMOSOME}_imputation" --output="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}_imputation.out" --error="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}_imputation.err" ${WHEREISCODE}/do_imputation23.sh $CHROMOSOME $WHEREISDATA; done


CHROMOSOME=$1
WHEREISDATA=$2

WHEREIS1000GENOMES="/projects/kg98/aurina/GWAS/data/enigma/1KGPref"
WHEREISIMPUTATIONQC="/projects/kg98/aurina/GWAS/data/DondersQC/imputationQC"

cd ${WHEREISDATA}

for ((j=1; j<=4; j++)); do
	if test -f ${WHEREISDATA}/chunk${j}-ready4mach.23.dat.gz; then
		./minimac2 --vcfReference --rounds 5 --states 200 --refHaps ${WHEREISIMPUTATIONQC}/reference_Xdata/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.fixed.gz --haps ${WHEREISDATA}/chunk${j}-ready4mach.${CHROMOSOME}.gz --snps ${WHEREISDATA}/chunk${j}-ready4mach.23.dat.gz.snps --autoClip ${WHEREISDATA}/autochunk-ready4mach.23.dat.gz --gzip --prefix ${WHEREISIMPUTATIONQC}/imputationResults/chunk${j}-ready4mach.23.female.imputed > ${WHEREISIMPUTATIONQC}/imputationResults/chunk${j}-ready4mach.23.female-minimac2.log
		cat ${WHEREISDATA}/chunk${j}-ready4mach.23.dat.gz.snps > ${WHEREISDATA}/temp.23.dat.gz.snps;
	fi
done

awk '!x[$0]++' ${WHEREISDATA}/temp.23.dat.gz.snps > ${WHEREISDATA}/ready4mach.23.dat.gz.snps
./minimac2 --vcfReference --rounds 5 --states 200 --refHaps ${WHEREISIMPUTATIONQC}/reference_Xdata/chr${CHROMOSOME}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.fixed.gz --haps ${WHEREISDATA}/ready4mach.${CHROMOSOME}.male.gz --snps ${WHEREISDATA}/ready4mach.23.dat.gz.snps --gzip --prefix ${WHEREISIMPUTATIONQC}/imputationResults/ready4mach.23.male.imputed > ${WHEREISIMPUTATIONQC}/imputationResults/ready4mach.23.male-minimac2.log
