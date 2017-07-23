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
#  1. [code]$ WHERESMYSCRIPT=$(pwd); CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23"; WHEREISDATA=${WHEREISIMPUTATIONQC}/phasingData
#  2. [....]$ for CHROMOSOME in $CHROMOSOMES; do sbatch --job-name="${CHROMOSOME}" --output="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}.out" --error="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}.err" ${WHERESMYSCRIPT}/do_phasing.sh $WHERESMYSCRIPT $CHROMOSOME $WHEREISDATA; done

WHERESMYSCRIPT=$1
CHROMOSOME=$2
WHEREISDATA=$3

WHEREISIMPUTATIONQC="/projects/kg98/aurina/GWAS/data/DondersQC/imputationQC"
cd ${WHEREISDATA}
echo -e "\nPHASING CHROMOSOME NR ${CHROMOSOME}\n"

for ((j=1; j<=4; j++)); do
    if test -f chunk${j}-ready4mach.${CHROMOSOME}.dat.gz; then
        ./mach1 -d chunk${j}-ready4mach.${CHROMOSOME}.dat.gz -p ready4mach.${CHROMOSOME}.ped.gz --prefix chunk${j}-ready4mach.${CHROMOSOME} --rounds 20 --states 200 --phase > chunk${j}-ready4mach.${CHROMOSOME}.mach.log
    fi
done
