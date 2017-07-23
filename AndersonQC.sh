#!/bin/bash
# copy relevant files to the directory qhere QC will be performed.
# files needed - original PedMap file, list of subjects to remove Remove.txt, updates phenotypes phenotype.txt and plink file.
module load R
module load matlab

WHEREISORIGINAL="/projects/kg98/aurina/GWAS/backupData/originalv2"
WHEREISQC="/projects/kg98/aurina/GWAS/data/andersonQC"
WHEREISCODE=""
cp ${WHEREISORIGINAL}/PedMap.ped ${WHEREISQC}
cp ${WHEREISORIGINAL}/PedMap.map ${WHEREISQC}
cp ${WHEREISORIGINAL}/Remove.txt ${WHEREISQC}
cp ${WHEREISORIGINAL}/phenotype.txt ${WHEREISQC}
cp ${WHEREISORIGINAL}/plink ${WHEREISQC}

# make relevant binary plink files with updated phenotype,remove selected subjects and allow subjects to have no sex information
cd ${WHEREISQC}

# -------------------------------------------------------------------------------------------------------
### step 0 - make binary PLINK files(replace phenotype, remove pre-defined subjects)
./plink --file PedMap --pheno phenotype.txt --remove Remove.txt --make-bed --out rawGWAdataRM
# -------------------------------------------------------------------------------------------------------
#########################################################################################################
# PRE-protocol filtering
# remove subjects with more then 10% missing data emove subjects with missig genotype >10%, SNPs with MAF <0.01 and SNPS with genotype rate <0.1
./plink --bfile rawGWAdataRM --mind 0.1 --maf 0.01 --geno 0.1 --out rawGWAdataRM

#########################################################################################################
### step 4
# identify individuals with discordant sex informaion
./plink --bfile rawGWAdataRM --check-sex --allow-no-sex --out rawGWAdataRM
# -------------------------------------------------------------------------------------------------------
### step 5
# Produce a list of individuals with discordant sex data
grep PROBLEM rawGWAdataRM.sexcheck > rawGWAdataRM.sexprobs
# -------------------------------------------------------------------------------------------------------
### step 6
# to select subjects that had initial sex information which doesn't correspond to the genotyped one, use script excludeSex.m it will create fail-sexcheck-qc.txt file.
matlab -nojvm < excludeSex.m
# -------------------------------------------------------------------------------------------------------
### step 7
# Identification of individuals with elevated missing data rates or outlying heterozygosity rate.
# create the files raw-GWA-data.imiss and raw-GWA-data.lmiss. The fourth column in the .imiss file (N_MISS) denotes the number of missing SNPs and the sixth column (F_MISS) denotes the proportion of missing SNPs per individual.
./plink --bfile rawGWAdataRM --missing --allow-no-sex --out rawGWAdataRM
# -------------------------------------------------------------------------------------------------------
### step 8
# create the file raw-GWA-data.het, in which the third column denotes the observed number of homozygous genotypes [O(Hom)] and the fifth column denotes the number of nonmissing genotypes [N(NM)] per individual.
./plink --bfile rawGWAdataRM --het --allow-no-sex --out rawGWAdataRM
# -------------------------------------------------------------------------------------------------------
### step 9
# Calculate the observed heterozygosity rate per individual using the formula (N(NM) − O(Hom))/N(NM). Create a graph in which the observed heterozygosity rate per individual is plotted on the x axis and the proportion of missing SNPs per individuals is plotted on the y axis. This can be carried out using standard software such as Excel or statistical packages such as SPSS. A script for calculating the heterozygosity rate and producing the graph using R is supplied (imiss-vs-het. Rscript). Type ‘R CMD BATCH imiss-vs-het.Rscript’ at the Unix prompt to run this script and create the graph (raw-GWA-data.imiss-vs-het.pdf).

R CMD BATCH imiss-vs-het.Rscript

# -------------------------------------------------------------------------------------------------------
### step 10
#Examine the plot to decide reasonable thresholds at which to exclude individuals based on elevated missing or extreme heterozygosity. We chose to exclude all individuals with a genotype failure rate ≥0.03 (Fig. 1, vertical dashed line) and/or a heterozygosity rate ± 3 s.d. from the mean (Fig. 1, horizontal dashed lines). Add the FID and IID of the samples failing this QC to the file named ‘fail-imisshet-qc.txt’.
# to calculate heterozigosity rate use matlab script excludeHeterozig.m it will create fail-imisshet-qc.txt
# change thresholds or other options in this script according to the raw-GWA-data.imiss-vs-het.pdf plot.
# to run it make a copy of rawGWAdataRM.imiss called rawGWAdataRM_copy.imiss and delete headlines . Do same for rawGWAdataRM.het --> rawGWAdataRM_copy.het
matlab -nojvm < excludeHeterozig.m

