#!/bin/bash
# copy relevant files to the directory qhere QC will be performed.
# files needed - original PedMap file, list of subjects to remove Remove.txt, updates phenotypes phenotype.txt and plink file.

module load R
module load matlab

mkdir /projects/kg98/aurina/GWAS/data/DondersQC
mkdir /projects/kg98/aurina/GWAS/data/DondersQC/MDS
mkdir /projects/kg98/aurina/GWAS/data/DondersQC/imputationQC
mkdir /projects/kg98/aurina/GWAS/data/enigma
mkdir /projects/kg98/aurina/GWAS/data/enigma/1KGPref

WHEREISORIGINAL="/projects/kg98/aurina/GWAS/backupData/originalv2"
WHEREISQC="/projects/kg98/aurina/GWAS/data/DondersQC"
WHEREISCODE="/projects/kg98/aurina/GWAS/code/dondersQC"
WHEREISENIGMA="/projects/kg98/aurina/GWAS/data/enigma"
WHEREISMDS="/projects/kg98/aurina/GWAS/data/DondersQC/MDS"
WHEREISIMPUTATIONQC="/projects/kg98/aurina/GWAS/data/DondersQC/imputationQC"
WHEREIS1000GENOMES="/projects/kg98/aurina/GWAS/data/enigma/1KGPref"


###DOWNLOAD ALL THE RELEVANT DATA and SOFTWARE###
#####################----------DATA------------#####################
##########-----HapMap3, 1000 genomes reference data------###########
###HapMap3 (save it in WHEREISENIGMA):
##http://enigma.ini.usc.edu/wp-content/uploads/2012/07/HM3.bed.gz
##http://enigma.ini.usc.edu/wp-content/uploads/2012/07/HM3.bim.gz
##http://enigma.ini.usc.edu/wp-content/uploads/2012/07/HM3.fam.gz
###1000 genomes EUR samples ## save to WHEREIS1000GENOMES
##http://enigma.ini.usc.edu/wp-content/uploads/2012/07/v3.20101123.ENIGMA2.EUR.20120719.vcf.tgz
##http://enigma.ini.usc.edu/wp-content/uploads/2012/07/v3.20101123.ENIGMA2.EUR.20120719.extras.tgz
##tar -zxvf v3.20101123.ENIGMA2.EUR.20120719.vcf.tgz
##tar -zxvf v3.20101123.ENIGMA2.EUR.20120719.extras.tgz


#####################---------SOFTWARE-----------#####################
##########-----PLINK, MACH, minimac2, ChunkChromosome------###########
### put all executable files in {WHEREISCODE} and add this directory to the path
###PLINK SOFTWARE## https://www.cog-genomics.org/plink2 (Linux 64 bit, save it in WHEREISCODE)
###MACH software## ## http://csg.sph.umich.edu/abecasis/MaCH/download/ to WHEREISCODE
## tar -zxvf mach.1.0.18.Linux.tgz
## cp ${WHEREISCODE}/executables/mach1 ${WHEREISCODE}
###ChunkCHROMOSOME software###
###new version downloaded from http://genome.sph.umich.edu/wiki/ChunkChromosome (link below)
## save it in ${WHEREISCODE}
##tar -zxvf generic-ChunkChromosome-2014-05-27.tar.gz
## install it using make install INSTALLDIR=${WHEREISCODE}
###### download minimac2 from http://genome.sph.umich.edu/wiki/Minimac2 into WHEREISCODE
##tar -zxvf minimac2.2014.9.15.tgz
#####################---------SNP id convert file-----------#####################
##https://support.illumina.com/array/downloads.html (Infinium PsychArray v1.1 Loci Name to rsID Conversion File)
## save it to ${WHEREISCODE}

PATH=${WHEREISCODE}
# copy original data and plink program
cp ${WHEREISORIGINAL}/PedMap.ped ${WHEREISQC}
cp ${WHEREISORIGINAL}/PedMap.map ${WHEREISQC}
cp ${WHEREISORIGINAL}/Remove.txt ${WHEREISQC}
cp ${WHEREISORIGINAL}/phenotype.txt ${WHEREISQC}
cp ${WHEREISORIGINAL}/plink ${WHEREISQC}

# make relevant binary plink files with updated phenotype,remove selected subjects
cd ${WHEREISQC}
./plink --file PedMap --pheno phenotype.txt --remove Remove.txt --allow-no-sex --make-bed --out rawGWAdata

# PRE-protocol filtering (remove really bad subjects and really bad SNPs)
# remove subjects with more then 10% missing data, remove subjects with missing genotype >10%,
# remove SNPs with MAF <0.01 and SNPS with genotype rate <0.1
# ./plink --bfile rawGWAdata --mind 0.1 --maf 0.01 --geno 0.1 --out rawGWAdata
# --miss (creates a table of missingness)

#*********************************
# Quality Control
# ********************************
### ------------------------------------------
### Sex check - exclude subjects that have mismatching information on sex
### ------------------------------------------

./plink --bfile ${WHEREISQC}/rawGWAdata --check-sex --allow-no-sex --out ${WHEREISQC}/rawGWAdata
# make a copy of .sexcheck file to work with it in matlab
cp ${WHEREISQC}/rawGWAdata.sexcheck ${WHEREISQC}/rawGWAdata.txt
# first create a list of subjects to remove that have initial sex information, but it missmatches with the one, calculated from SNPs using matlab
# it will make a file called fail-sex-qc.txt
matlab -nojvm < excludeSexMissmatch.m
# Exclude those individuals
./plink --bfile ${WHEREISQC}/rawGWAdata --remove ${WHEREISQC}/fail-sex-qc.txt --allow-no-sex --make-bed --out ${WHEREISQC}/rawGWAdata_rsex1
# impute sex information for the ones that are left (subjects with matching info or subjects with no infor to start with).
./plink --bfile ${WHEREISQC}/rawGWAdata_rsex1 --impute-sex --make-bed --out ${WHEREISQC}/rawGWAdata_rsex

# Make list of individuals with PROBLEM in sex-check in order to exclude
#fgrep -f <(grep "PROBLEM" ${WHEREISQC}/rawGWAdata.sexcheck | awk '{print $2}') ${WHEREISQC}/rawGWAdata.fam | awk 'BEGIN{OFS="\t"} {print $1,$2}' > ${WHEREISQC}/fail-sex-QC.txt

### ----------------------------------------
### Filter SNPs according to HWE, MAF, genotyping frequency
### ----------------------------------------

./plink --bfile ${WHEREISQC}/rawGWAdata_rsex --hwe 1e-6 --geno 0.05 --maf 0.01 --allow-no-sex --make-bed --out ${WHEREISQC}/rawGWAdata_rsex_filtered

### ----------------------------------------
### Filtersubjects that have low genotyping rate (0.02??? or 0.05???)
### ----------------------------------------

./plink --noweb --bfile ${WHEREISQC}/rawGWAdata_rsex_filtered --mind 0.05 --allow-no-sex --make-bed --out ${WHEREISQC}/rawGWAdata_rsex_filteredID

### ---------------------------------------------------------
### Check for cryptic relatedness and duplicate samples (IBD)
### ---------------------------------------------------------

# Exclude markers on sex chromosomes and chromosomes 25 and 26
grep '^23\|^24\|^25\|^26' ${WHEREISQC}/rawGWAdata_rsex_filteredID.bim  | cut -f 2 > ${WHEREISQC}/chr23_26.txt
# Exclude these markers and at the same time, prune the markers based on LD, and some other filters
./plink --bfile ${WHEREISQC}/rawGWAdata_rsex_filteredID  --exclude ${WHEREISQC}/chr23_26.txt --indep-pairwise 50 5 0.25 --allow-no-sex --out ${WHEREISQC}/common_LD_variants
#qsub bin/pruning.job
## As a result, we get 81601 variants in results/preqc/common_LD_variants.prune.in
# We use common_LD_variants.prune.in to select the variants that we want to use for making the genome file (for checking the PI-HAT value):
./plink --bfile ${WHEREISQC}/rawGWAdata_rsex_filteredID --extract ${WHEREISQC}/common_LD_variants.prune.in --genome --allow-no-sex --out ${WHEREISQC}/rawGWAdata_rsex_filteredID
#qsub bin/genome.job
# We filter the output on PIHat and keep only the Pairs and PIHAt values. Perform sorting on Pihat

awk 'BEGIN{OFS="\t"} {if ($10 >= 0.25) print $2,$4,$10}' ${WHEREISQC}/rawGWAdata_rsex_filteredID.genome | sort -k3,3gr > ${WHEREISQC}/rawGWAdata_rsex_filteredID.genome.PIHat.sorted.txt
# with matlab script select one from a pair of related subjects and save their FID and IID in a file called fail-IBD-QC.txt
# there are 7 subjects that are relatives of more than one person
matlab -nojvm < excludeIBD.m
### remove samples from this list of related subjects
./plink --bfile ${WHEREISQC}/rawGWAdata_rsex_filteredID --remove ${WHEREISQC}/fail-IBD-QC.txt --allow-no-sex --make-bed --out ${WHEREISQC}/rawGWAdata_rsex_filteredID_unrelated
### ------------------------------------------
### Check for heterozigosity deviations
### ------------------------------------------
# do this step for controls only
./plink --bfile rawGWAdata_rsex_filteredID --filter-controls --make-bed --out rawGWAdata_rsex_filteredIDcontrols
./plink --bfile rawGWAdata_rsex_filteredIDcontrols --het --allow-no-sex --out rawGWAdata_rsex_filteredIDcontrols

#make a copy of this file with txt extension for matlab
cp rawGWAdata_rsex_filteredIDcontrols.het rawGWAdata_rsex_filteredIDcontrols.txt
# plot distributios in matlab to see distributions and to choose appropriate thresholds
# this script will create fail-imisshet-qc.txt; use this file to exlucde subjects.
# 16 subjects from controls excluded - not sure if this step necessary, deciations from the mean are not huge...
matlab -nojvm < excludeHeterozig.m
# remove those subjects
./plink -bfile rawGWAdata_rsex_filteredID --remove fail-imisshet-qc.txt --make-bed --out rawGWAdata_rsex_filteredIDhet

### ------------------------------------------
### Prepare data for MDS
### ------------------------------------------
cp ${WHEREISQC}/plink ${WHEREISMDS}
./plink --noweb --bfile ${WHEREISQC}/rawGWAdata_rsex_filteredIDhet --hwe 1e-6 --geno 0.05 --maf 0.01 --make-bed --out ${WHEREISMDS}/mds_filtered
# do this on pruned set of SNPs
# From HM3 only keep SNPs in our dataset.
cut -f 2 ${WHEREISMDS}/mds_filtered.bim > ${WHEREISQC}/HM3.snplist.txt
./plink --bfile ${WHEREISQC}/HM3 --extract ${WHEREISQC}/HM3.snplist.txt --make-bed --out ${WHEREISQC}/HM3_subset

# Prune the SNPs in HM3
./plink --bfile ${WHEREISQC}/HM3_subset --indep 50 5 1.15 --out ${WHEREISQC}/HM3_subset
./plink --bfile ${WHEREISQC}/HM3_subset --extract ${WHEREISQC}/HM3_subset.prune.in --make-bed --out ${WHEREISQC}/HM3_subset_pruned

# Extract the pruned SNPs from our dataset
./plink --bfile ${WHEREISMDS}/mds_filtered --extract ${WHEREISQC}/HM3_subset.prune.in --make-bed --out ${WHEREISMDS}/subset_variants_pruned

# Create the files for updating the positions in order that both datasets contain the same and one file for updating the reference allele for the same reason
awk '{print($2"\t"$1)}' ${WHEREISQC}/HM3_subset_pruned.bim > ${WHEREISMDS}/hm3_CHR_update
awk '{print($2"\t"$4)}' ${WHEREISQC}/HM3_subset_pruned.bim > ${WHEREISMDS}/hm3_BP_update
awk '{print($2"\t"$5)}' ${WHEREISQC}/HM3_subset_pruned.bim > ${WHEREISMDS}/hm3_allele_update

# Now we will update the position of our data considering the data from hapmap
./plink --bfile ${WHEREISMDS}/subset_variants_pruned --update-map ${WHEREISMDS}/hm3_CHR_update --update-chr --make-bed --out ${WHEREISMDS}/subset_variants_pruned_chr
./plink --bfile ${WHEREISMDS}/subset_variants_pruned_chr --update-map ${WHEREISMDS}/hm3_BP_update --make-bed --out ${WHEREISMDS}/subset_variants_pruned_pos

# Force the reference allele in our data using the data from hapmap
./plink --bfile ${WHEREISMDS}/subset_variants_pruned_pos --reference-allele ${WHEREISMDS}/hm3_allele_update --make-bed --out ${WHEREISMDS}/subset_variants_pruned_forced
## 20496 could not be forced (find out with grep -c "Warning: Impossible A1 allele assignment for variant" results/mds/subset_variants_pruned_forced.log)

# Get which ones have the problem
awk 'NR==FNR{a[$2]=$5;next}{ if ($2 in a) print $0"\t"a[$2]; else print $0"\t.";}' ${WHEREISMDS}/subset_variants_pruned_forced.bim ${WHEREISQC}/HM3_subset_pruned.bim | awk '{if ($5 != $7) print $2;}' > ${WHEREISMDS}/flip_snps.txt
# Flip these
./plink --bfile ${WHEREISMDS}/subset_variants_pruned_pos --flip ${WHEREISMDS}/flip_snps.txt --make-bed --out ${WHEREISMDS}/subset_variants_pruned_flip

# Re-try to force the alleles now that we flipped
./plink --bfile ${WHEREISMDS}/subset_variants_pruned_flip --reference-allele ${WHEREISMDS}/hm3_allele_update --make-bed --out ${WHEREISMDS}/subset_variants_pruned_forced_2

# Merge our dataset with hapmap
./plink --bfile ${WHEREISMDS}/subset_variants_pruned_forced_2 --bmerge ${WHEREISQC}/HM3_subset_pruned.bed ${WHEREISQC}/HM3_subset_pruned.bim ${WHEREISQC}/HM3_subset_pruned.fam --make-bed --out ${WHEREISMDS}/rawGWAdata_HM3_merged
# DONE!

# Make genome file
./plink --bfile ${WHEREISMDS}/rawGWAdata_HM3_merged --genome --out ${WHEREISMDS}/rawGWAdata_HM3_merged

# Perfom MDS
./plink --bfile ${WHEREISMDS}/rawGWAdata_HM3_merged --read-genome ${WHEREISMDS}/rawGWAdata_HM3_merged.genome  --mind .05 --cluster --mds-plot 10 --out ${WHEREISMDS}/rawGWAdata_HM3_mds

# modify mds file to add our group as ethnicity in the first column
#sed 's/^[\t ]\+//' ${WHEREISMDS}/rawGWAdata_HM3_mds.mds | sed 's/^[0-9]\+/PBSN/' > ${WHEREISMDS}/rawGWAdata_HM3_mds.withGroups.mds

### Make plot in R
cd ${WHEREISMDS}

module load R
R
install.packages("calibrate")
library(calibrate)

mds.cluster = read.table("rawGWAdata_HM3_mds.mds", header=T);
num = length(mds.cluster$C1) - 988;
colors = rev(c(rep("red", num), rep("lightblue", 49), rep("brown", 112), rep("yellow",
84),rep("green", 85), rep("purple", 88), rep("orange", 86), rep("grey50", 90), rep("black",
50),rep("darkolivegreen", 143), rep("magenta", 88), rep("darkblue", 113)));

pdf(file="mdsplot.pdf",width=7,height=7)

plot(rev(mds.cluster$C2), rev(mds.cluster$C1), col=colors, ylab="Dimension 1",
xlab="Dimension 2",pch=20)
legend("bottomleft", c("My Sample", "ASW", "CEU", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX",
"MKK","TSI", "YRI"), fill=c("red", "lightblue", "brown", "yellow", "green", "purple", "orange",
"grey50", "black", "darkolivegreen", "magenta", "darkblue"))
dev.off()

# make a copy of .mds file to work with in matlab
cp ${WHEREISMDS}/rawGWAdata_HM3_mds.mds ${WHEREISMDS}/rawGWAdata_HM3_mds.txt
cp ${WHEREISQC}/rawGWAdata.fam ${WHEREISMDS}/rawGWAdata.txt
# make a list of subjects to exclude
cd ${WHEREISMDS}
matlab -nojvm < excludeMDS.m
# using 2SD exclusion criteria, 56 subjects areto be excluded.
# exclude them
./plink --bfile ${WHEREISMDS}/mds_filtered --remove ${WHEREISMDS}/fail-mds-qc.txt --make-bed --out ${WHEREISIMPUTATIONQC}/rawGWAdata_afterMDS
# check call rates for SNPs and people after MDS perople excluded. all in one command


### -----------
### IMPUTATION
### -----------
####### QUALITY CONTROL #######

./plink --bfile ${WHEREISIMPUTATIONQC}/rawGWAdata_afterMDS –-maf 0.01 --geno 0.05 --hwe 0.000001 --make-bed --out ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4
# qsub bin/qc_excl_markers.sh
# qsub bin/qc_filter_markers.sh

####### PREPARATION #######
cd ${WHEREISIMPUTATIONQC}

mkdir preparationData
mkdir preparationResults

### find common markers in our dataset and the reference panel and compare
ln -s ${WHEREIS1000GENOMES}/1kgp.* preparationData

## we need to couple the snps to the rs number
# copy rskey file
cp ${WHEREISCODE}/InfiniumPsychArray-24v1-1_A1_b138_rsids.txt ${WHEREISIMPUTATIONQC}/preparationData/
awk '{print $2, $1, $3, $4, $5, $6}' ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.bim > ${WHEREISIMPUTATIONQC}/preparationResults/temp.bim
awk 'FNR==NR{a[$1]=$2;next} $1 in a{print a[$1],$2,$3,$4,$5,$6}' ${WHEREISIMPUTATIONQC}/preparationData/InfiniumPsychArray-24v1-1_A1_b138_rsids.txt ${WHEREISIMPUTATIONQC}/preparationResults/temp.bim > ${WHEREISIMPUTATIONQC}/preparationResults/predatfile.bim
awk '{print $2,$1,$3,$4,$5,$6}' ${WHEREISIMPUTATIONQC}/preparationResults/predatfile.bim > ${WHEREISIMPUTATIONQC}/preparationResults/rsnumbers.bim
# rename and place with the other corresponding plink files
mv ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.bim ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4_withoutRS.bim
mv ${WHEREISIMPUTATIONQC}/preparationResults/rsnumbers.bim ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.bim

#cp ${WHEREISIMPUTATIONQC}/superclean_step1.* ${WHEREISIMPUTATIONQC}/preparationResults/
#mv ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step1.fam ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.fam
#mv ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step1.bed ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.bed
# join the two files
awk '{print $2,$1,$3,$4,$5,$6}' ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.bim > ${WHEREISIMPUTATIONQC}/preparationResults/tempQC.bim
awk 'NR==FNR{s=$1;a[s]=$0;next} a[$1]{print $0 " "a[$1]}' ${WHEREISIMPUTATIONQC}/preparationResults/tempQC.bim ${WHEREISIMPUTATIONQC}/preparationData/1kgp.alleles > ${WHEREISIMPUTATIONQC}/preparationResults/merged.alleles
# select list of snps that need flipping
awk '{ if ($2!=$8 && $2!=$9) print $1}' ${WHEREISIMPUTATIONQC}/preparationResults/merged.alleles > ${WHEREISIMPUTATIONQC}/preparationResults/flip.list
#!!! # compare alleles
#!!! awk '{ if ($2==$8 && $3==$9) print $0"\tcorrect"; else if ($2==$9 && $3==$8) print $0"\tforceRef"; else print $0"\terror"}' results/preparation_2/merged.alleles > results/preparation_2/merged.comparisons.txt

# do the flip and update chrom
# get a list of duplicate IDs
./plink --bfile ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4 --list-duplicate-vars suppress-first --out ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4
cut -f 4 ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.dupvar> ${WHEREISIMPUTATIONQC}/preparationResults/snps2remove.txt
# NOTE: add 4 SNPs manually rs7120775, rs35873579, rs1800194, rs12814239 afterwards to snps2remove.txt list - they were not labeled as duplicates in the first command (don't know why)

# remove them
./plink --bfile ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4 --exclude ${WHEREISIMPUTATIONQC}/preparationResults/snps2remove.txt --make-bed --out ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.1
cp ${WHEREISQC}/plink ${WHEREISIMPUTATIONQC}
./plink --bfile ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.1 --extract ${WHEREISIMPUTATIONQC}/preparationData/1kgp.snps --update-map ${WHEREISIMPUTATIONQC}/preparationData/1kgp.chr --flip ${WHEREISIMPUTATIONQC}/preparationResults/flip.list --update-chr --make-bed --out ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.5
#./plink --bfile ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4 --extract ${WHEREISIMPUTATIONQC}/preparationData/1kgp.snps --update-map ${WHEREISIMPUTATIONQC}/preparationData/1kgp.chr --update-chr --make-bed --out ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.5
# here we lost some markers because there were many that were duplicate

# update position

./plink --bfile ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step4.5 --update-map ${WHEREISIMPUTATIONQC}/preparationData/1kgp.bp --geno 0.05 --mind 0.05 --make-bed --out ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step5

# ### In the previous step we lost the X chromosome
# awk 'NR==FNR{s=$1;a[s]=$0;next} {if ($1 in a) print $0"\t"a[$1]}' results/preparation/superclean_step4.chr23.alleles data/reference_X/1kgp.23.alleles > results/preparation/step4_1kgphase3.chr23.alleles
# # compare alleles
# awk '{ if ($2==$5 && $3==$6) print $7"\tcorrect"; else if ($2==$6 && $3==$5) print $7"\tforceRef"; else print $7"\terror"}' results/preparation/step4_1kgphase3.chr23.alleles > results/preparation/step4_1kgphase3.chr23.comparisons.txt
# # those that have "error" we will try to flip and keep only those in common
# cut -f 7 results/preparation/step4_1kgphase3.chr23.alleles > results/preparation/common.23.snps
# grep "error" results/preparation/step4_1kgphase3.chr23.comparisons.txt | cut -f 1 > results/preparation/flip.23.list
# plink --bfile results/preparation/superclean_step4 --extract results/preparation/common.23.snps --flip results/preparation/flip.23.list --make-bed --noweb --out results/preparation/superclean_step4_flipped_chr23
# # recheck how many errors we have now
# awk '{print $1":"$4"\t"$5"\t"$6"\t"$2}' results/preparation/superclean_step4_flipped_chr23.bim > results/preparation/superclean_step4_flipped.23.alleles
# awk 'NR==FNR{s=$1;a[s]=$0;next} {if ($1 in a) print $0"\t"a[$1]}' results/preparation/superclean_step4_flipped.23.alleles data/reference_X/1kgp.23.alleles > results/preparation/step4_flipped_1kgphase3.23.alleles
# awk '{ if ($2==$5 && $3==$6) print $7"\tcorrect"; else if ($2==$6 && $3==$5) print $7"\tforceRef"; else print $7"\terror"}' results/preparation/step4_flipped_1kgphase3.23.alleles > results/preparation/step4_flipped_1kgphase3.chr23.comparisons.txt
# # now 62 have an error we will remove these
# grep  "error" results/preparation/step4_flipped_1kgphase3.chr23.comparisons.txt | cut -f 1 > results/preparation/error_markers.23.txt
# plink --noweb --bfile results/preparation/superclean_step4_flipped_chr23 --exclude results/preparation/error_markers.23.txt --make-bed --out results/preparation/superclean_step5_23

# ### make lists of males and females
awk '{ if($5==1) print $1, $2}' ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step5.fam > ${WHEREISIMPUTATIONQC}/preparationResults/male.list
awk '{ if($5==2) print $1, $2}' ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step5.fam > ${WHEREISIMPUTATIONQC}/preparationResults/female.list


####### MACH PHASING #######
mkdir phasingData
mkdir phasingResults
mkdir bin
mkdir /projects/kg98/aurina/GWAS/code/errorLog/

### Split per chromosome (submit a job per chromosome, runs about 2-4s each)
WHEREISDATA=${WHEREISIMPUTATIONQC}/preparationResults
for i in {1..23}; do
./plink --bfile ${WHEREISDATA}/superclean_step5 --chr ${i} --recode --out ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.${i}
done

#WHERESMYSCRIPT=${WHEREISCODE}; CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"; WHEREISDATA=${WHEREISIMPUTATIONQC}/preparationResults
#for CHROMOSOME in $CHROMOSOMES; do sbatch --job-name="${CHROMOSOME}" --output="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}.out" --error="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}.err" ${WHERESMYSCRIPT}/split_chromosomes.sh $WHERESMYSCRIPT $CHROMOSOME $WHEREISDATA; done

# for 23 chromosome:
./plink --bfile ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step5 --chr 23 --set-hh-missing --keep ${WHEREISIMPUTATIONQC}/preparationResults/male.list --recode --out ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.23.male
./plink --bfile ${WHEREISIMPUTATIONQC}/preparationResults/superclean_step5 --chr 23 --keep ${WHEREISIMPUTATIONQC}/preparationResults/female.list --recode --out ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.23.female

#
mv ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.23.female.map ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.23.map
rm ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.23.male.map

### Reformat dat
# create jobfiles to reformat each chrom's file
for i in {1..23}; do
echo "S dummy" > ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.${i}.dat
awk '{ print "M", $1 ":" $4}' ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.${i}.map >> ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.${i}.dat
gzip ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.${i}.dat
gzip ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.${i}.ped
done

### chunking the chromosomes and phasing
mkdir ${WHEREISCODE}/ChunkChromosome
## NOTE: ChunkChromosome old version from coocbook isn't there.
## NOTE: new version downloaded from http://genome.sph.umich.edu/wiki/ChunkChromosome (link below)
## save it in ${WHEREISCODE}/ChunkChromosome
## install it using make install INSTALLDIR=${WHEREISIMPUTATIONQC}/phasingData/
## use directory where data is located
## tar -zxvf generic-ChunkChromosome-2014-05-27.tar.gz
##cd genericChunkChromosome/
##make install INSTALLDIR=${WHEREISIMPUTATIONQC}/phasingData/
cd ${WHEREISIMPUTATIONQC}/phasingData/
### chunking the chromosomes and phasing
# must cd there for the next step to work
for i in {1..23}; do
./ChunkChromosome -d ready4mach.${i}.dat.gz -n 5000 -o 500;
done

## download Linux-mach from here http://csg.sph.umich.edu/abecasis/MaCH/download/ to /projects/kg98/aurina/GWAS/code/MACH
tar -zxvf mach.1.0.18.Linux.tgz
## copy executable file mach1 to where your data is
cp /projects/kg98/aurina/GWAS/code/MACH/executables/mach1 ${WHEREISCODE}/mach1

### MaCH phasing
# prepare the scripts run from where the script is
# NOTE!!: manually edit jobfile for chromosome 23 in order to use the correct ped file: ready4mach.23.female.ped as well as the output prefix
# submit
WHERESMYSCRIPT=$(pwd); CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23"; WHEREISDATA=${WHEREISIMPUTATIONQC}/phasingData
for CHROMOSOME in $CHROMOSOMES; do sbatch --job-name="${CHROMOSOME}" --output="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}.out" --error="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}.err" ${WHERESMYSCRIPT}/do_phasing.sh $WHERESMYSCRIPT $CHROMOSOME $WHEREISDATA; done
### HOW TO MAKE GOOD 23.female file?????
### MAYBE: remane ready4mach.23.female.ped to ready4mach.23.ped (because female .map was renamed to .map in the scripts)



### setting up the male X chromosome data for imputation
gzip ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.23.male.ped
zcat ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.23.male.ped.gz | awk '{printf "%s", $1 "->"$2 " HAPLO1 "; for (N=7; N<=NF; N+=2) printf "%s", $N; printf "\n"; printf "%s", $1 "->"$2 " HAPLO2 "; for (N=7; N<=NF; N+=2) printf "%s", $N; printf "\n";}' > ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.23.male
gzip ${WHEREISIMPUTATIONQC}/phasingData/ready4mach.23.male

#### put all software downloads at the beginning of the script and add executables to one folder and add it to path to begin with
####### IMPUTATION #######

###### download minimac2 from http://genome.sph.umich.edu/wiki/Minimac2 into code/dondersQC/Minimac2
tar -zxvf minimac2.2014.9.15.tgz
# copy file where your data is
cp /projects/kg98/aurina/GWAS/code/minimac2/minimac2/bin/minimac2 ${WHEREISCODE}

mkdir ${WHEREISIMPUTATIONQC}/imputationResults
mkdir ${WHEREISIMPUTATIONQC}/reference_Xdata
# fix the female reference vcf because it has spaces instead of tabs
zgrep "^#" ${WHEREIS1000GENOMES}/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.gz > ${WHEREISIMPUTATIONQC}/reference_Xdata/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.fixed.gz
zgrep -v "^#" ${WHEREIS1000GENOMES}/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.gz | tr ' ' '\t' >> ${WHEREISIMPUTATIONQC}/reference_Xdata/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.fixed.gz


##RUN do_imputation here for 22 chromosomes
WHERESMYSCRIPT=$(pwd); CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"; WHEREISDATA=${WHEREISIMPUTATIONQC}/phasingData
for CHROMOSOME in $CHROMOSOMES; do sbatch --job-name="${CHROMOSOME}_imputation" --output="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}_imputation.out" --error="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}_imputation.err" ${WHEREISCODE}/do_imputation.sh $CHROMOSOME $WHEREISDATA; done

##RUN do_imputation23 here for 23 chromosome (not necessary, not used in other steps anyway)
WHERESMYSCRIPT=$(pwd); CHROMOSOMES="23"; WHEREISDATA=${WHEREISIMPUTATIONQC}/phasingData
for CHROMOSOME in $CHROMOSOMES; do sbatch --job-name="${CHROMOSOME}_imputation" --output="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}_imputation.out" --error="/projects/kg98/aurina/GWAS/code/dondersQC/errorLog/slurm-${CHROMOSOME}_imputation.err" ${WHEREISCODE}/do_imputation23.sh $CHROMOSOME $WHEREISDATA; done
### all done using minimac2

rm ${WHEREISIMPUTATIONQC}/imputationResults/*draft
