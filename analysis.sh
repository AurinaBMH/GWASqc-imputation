### Author: Marc Pauper
### Date: October 2016

## Data from 4 projects will be imputed: PIB_BOA_SPI_NLS
## Finally, imputed data can be split in separate files per project.

## Imputation done according to ENIGMA2 1KG cookbook (http://enigma.ini.usc.edu/wp-content/uploads/2012/07/ENIGMA2_1KGP_cookbook_v3.pdf)

## Imputation working directory: /home/mpauper/workspace/imputation_PIB_BOA_SPI_NLS

# Copy input raw data
cp /home/multifac/PsychChip_datafeb2016_PIB_BOA_SPI_NLS/Reports.ped data/
cp /home/multifac/PsychChip_datafeb2016_PIB_BOA_SPI_NLS/Reports.map data/

# Convert to binary ped
plink --noweb --file data/Reports --make-bed --out data/Reports
#qsub bin/make_binary.job

# To save space, we remove the .ped and .map
rm data/Reports.ped data/Reports.map


#*********************************
# PsychChip cleaning
# ********************************

mkdir results/PsychChip_cleaning

python data/PsychChip_cleaning_info/standardizePlinkBfile_bugfix.py \
    --out results/PsychChip_cleaning/pchip_cleaned \
    --sloc temp/ \
    --loloc doc/ \
    --p2loc /home/multifac/bin/plink \
    --exclude data/PsychChip_cleaning_info/pchip_blackList_dec2015_stripped.txt \
    --extract data/PsychChip_cleaning_info/PsychChip_snps_singlePerfectMatches_160216.txt \
    --chr-update data/PsychChip_cleaning_info/pchip_chrupdate_160219.txt \
    --pos-update data/PsychChip_cleaning_info/pchip_posupdate_160219.txt \
    --name-update data/PsychChip_cleaning_info/PsychChip_rsid_update_160202.txt \
    --bed data/Reports.bed --fam data/Reports.fam --bim data/Reports.bim
#qsub bin/pchip_cleaning.job

# Excluded 12968 snps from file pchip_blackList_dec2015_stripped.txt
# Excluded 6913 snps not in file PsychChip_snps_singlePerfectMatches_160216.txt
# Updated the chromosome of 99 snps from file pchip_chrupdate_160219.txt
# Updated the position of 3193 snps from file pchip_posupdate_160219.txt
# Updated the name of 123 snps from file PsychChip_rsid_update_160202.txt


#*********************************
# Quality Control
# ********************************

mkdir results/preqc

## First, we filter on SNP call rate, to exclude bad markers that may affect subsequent steps
plink --noweb --bfile results/PsychChip_cleaning/pchip_cleaned --geno 0.02 --make-bed --out results/preqc/preqc_geno
# qsub bin/preqc_geno.job
# Here we lose almost all markers, so we need to check a plot of various missingness thresholds and proportion of SNPs lost

### ------------------------------------------
### We try gender check first
### ------------------------------------------


plink --noweb --bfile results/PsychChip_cleaning/pchip_cleaned --check-sex --out results/preqc/preqc_sexcheck
### For some, the F-value is not calculated ('nan'). This is probably because the data on chromosome X is missing. We manually check to make sure:
# select one inidivual with 'nan' randomly: NL93PIBUNJD100182
# make a file with FID and IID
grep NL93PIBUNJD100182 results/PsychChip_cleaning/pchip_cleaned.fam | awk 'BEGIN{OFS="\t"} {print $1,$2}' > temp/NL93PIBUNJD100182.txt
# make a file with only markers from chromosome X to extract
grep "^23" results/PsychChip_cleaning/pchip_cleaned.bim | cut -f 2 > temp/chrX.markers.txt
# extract this individual and markers and make new plink files
plink --noweb --bfile results/PsychChip_cleaning/pchip_cleaned --extract temp/chrX.markers.txt --keep temp/NL93PIBUNJD100182.txt --recode --out temp/NL93PIBUNJD100182_chrX
## then we have a look at temp/NL93PIBUNJD100182_chrX.ped and see that all genotypes are set to 0.
## This is why F-value in sex-check is not calculated.
# Make list of individuals with PROBLEM in sex-check in order to exclude
fgrep -f <(grep "PROBLEM" results/preqc/preqc_sexcheck.sexcheck | awk '{print $2}') results/PsychChip_cleaning/pchip_cleaned.fam | awk 'BEGIN{OFS="\t"} {print $1,$2}' > results/preqc/preqc_sexcheck_exclude.txt
# Exclude the individuals
plink --noweb --bfile results/PsychChip_cleaning/pchip_cleaned --remove results/preqc/preqc_sexcheck_exclude.txt --make-bed --out results/preqc/preqc_sex_excluded


### ----------------------------------------
### Let's try to exclude SNP call rate again
### ----------------------------------------

plink --noweb --bfile results/PsychChip_cleaning/preqc_sex_excluded --geno 0.02 --make-bed --out results/preqc/preqc_geno
# This time we only lose 7482 variants



### ---------------------------------------------------------
### Check for cryptic relatedness and duplicate samples (IBD)
### ---------------------------------------------------------

# Exclude markers on sex chromosomes and chromosomes 25 and 26
grep '^23\|^24\|^25\|^26' results/preqc/preqc_geno.bim  | cut -f 2 > results/preqc/chr23_26.txt
# Exclude these markers and at the same time, prune the markers based on LD, and some other filters
plink --noweb --bfile results/preqc/preqc_geno --maf 0.05 --exclude results/preqc/chr23_26.txt --indep-pairwise 50 5 0.25 --geno 0.01 --hwe 0.0001 --out results/preqc/common_LD_variants
#qsub bin/pruning.job
## As a result, we get 81601 variants in results/preqc/common_LD_variants.prune.in
# We use common_LD_variants.prune.in to select the variants that we want to use for making the genome file (for checking the PI-HAT value):
plink --noweb --bfile results/preqc/preqc_geno --extract results/preqc/common_LD_variants.prune.in --genome --out results/preqc/preqc_geno_pruned
#qsub bin/genome.job
# We filter the output on PIHat and keep only the Pairs and PIHAt values. Perform sorting on Pihat
awk 'BEGIN{OFS="\t"} {if ($10 >= 0.25) print $2,$4,$10}' results/preqc/preqc_geno_pruned.genome | sort -k3,3gr > results/preqc/preqc_geno_pruned.genome.PIHat.sorted.txt
# With a python script, I add the project of each sample in each pair
python bin/IBD_pairs_assign_project.py
# In excel I highlighted the pairs that do not belong to the same project and sent the results
### After some discussion and asking the ENIGMA2 Genetics Team, we decide that we can impute even with the related individuals
### However, I am going to remove the complete duplicates and keep only one from each pair
### samples to exclude are put in: results/preqc/preqc_geno_pruned.duplicates.remove.txt
plink --noweb --bfile results/preqc/preqc_geno --remove results/preqc/preqc_geno_pruned.duplicates.remove.txt --make-bed --out results/preqc/preqc_dups_removed



### -----------------------
### Filter sample call rate
### -----------------------

plink --noweb --bfile results/preqc/preqc_dups_removed --mind 0.05 --make-bed --out results/preqc/preqc_sample_rate
# 0 individuals removed



### -------------------------
### Multi-dimensional Scaling
### -------------------------

mkdir results/mds
# Initial filtering
plink --noweb --bfile results/preqc/preqc_sample_rate --hwe 1e-6 --geno 0.05 --maf 0.01 --make-bed --out results/mds/mds_filtered
# 0 variants removed due to missing genotype data (--geno).
# --hwe: 788 variants removed due to Hardy-Weinberg exact test.
# 271114 variants removed due to minor allele threshold(s)

# Make links to HapMap3 data
mkdir data/HM3
ln -s /home/multifac/Software_utilities/HapMap/HM3.bed data/HM3/
ln -s /home/multifac/Software_utilities/HapMap/HM3.bim data/HM3/
ln -s /home/multifac/Software_utilities/HapMap/HM3.fam data/HM3/

# From HM3 only keep SNPs in our dataset.
cut -f 2 results/mds/mds_filtered.bim > data/HM3/HM3.snplist.txt
plink --noweb --bfile data/HM3/HM3 --extract data/HM3/HM3.snplist.txt --make-bed --out data/HM3/HM3_subset

# Prune the SNPs in HM3
plink --noweb --bfile data/HM3/HM3_subset --indep 50 5 1.15 --out data/HM3/HM3_subset
plink --noweb --bfile data/HM3/HM3_subset --extract data/HM3/HM3_subset.prune.in --make-bed --out data/HM3/HM3_subset_pruned

# Extract the pruned SNPs from our dataset
plink --noweb --bfile results/mds/mds_filtered --extract data/HM3/HM3_subset.prune.in --make-bed --out results/mds/subset_variants_pruned

# Create the files for updating the positions in order that both datasets contain the same and one file for updating the reference allele for the same reason
awk '{print($2"\t"$1)}' data/HM3/HM3_subset_pruned.bim > results/mds/hm3_CHR_update
awk '{print($2"\t"$4)}' data/HM3/HM3_subset_pruned.bim > results/mds/hm3_BP_update
awk '{print($2"\t"$5)}' data/HM3/HM3_subset_pruned.bim > results/mds/hm3_allele_update

# Now we will update the position of our data considering the data from hapmap
plink --noweb --bfile results/mds/subset_variants_pruned --update-map results/mds/hm3_CHR_update --update-chr --make-bed --out results/mds/subset_variants_pruned_chr
plink --noweb --bfile results/mds/subset_variants_pruned_chr --update-map results/mds/hm3_BP_update --make-bed --out results/mds/subset_variants_pruned_pos

# Force the reference allele in our data using the data from hapmap
plink --noweb --bfile results/mds/subset_variants_pruned_pos --reference-allele results/mds/hm3_allele_update --make-bed --out results/mds/subset_variants_pruned_forced
## 20496 could not be forced (find out with grep -c "Warning: Impossible A1 allele assignment for variant" results/mds/subset_variants_pruned_forced.log)

# Get which ones have the problem
awk 'NR==FNR{a[$2]=$5;next}{ if ($2 in a) print $0"\t"a[$2]; else print $0"\t.";}' results/mds/subset_variants_pruned_forced.bim data/HM3/HM3_subset_pruned.bim | awk '{if ($5 != $7) print $2;}' > results/mds/flip_snps.txt
# Flip these
plink --noweb --bfile results/mds/subset_variants_pruned_pos --flip results/mds/flip_snps.txt --make-bed --out results/mds/subset_variants_pruned_flip

# Re-try to force the alleles now that we flipped
plink --noweb --bfile results/mds/subset_variants_pruned_flip --reference-allele results/mds/hm3_allele_update --make-bed --out results/mds/subset_variants_pruned_forced_2
# NO ERRORS :)

# Merge our dataset with hapmap
plink --noweb --bfile results/mds/subset_variants_pruned_forced_2 --bmerge data/HM3/HM3_subset_pruned.bed data/HM3/HM3_subset_pruned.bim data/HM3/HM3_subset_pruned.fam --make-bed --out results/mds/rawGWAdata_HM3_merged
# DONE!

# Make genome file
# plink --noweb --bfile results/mds/PIB_BOA_SPI_NLS_HM3_merged --genome --out results/mds/PIB_BOA_SPI_NLS_HM3_merged
qsub bin/mds_mkgenome.sh

# Perfom MDS
# plink --noweb --bfile results/mds/PIB_BOA_SPI_NLS_HM3_merged --read-genome results/mds/PIB_BOA_SPI_NLS_HM3_merged.genome  --mind .05 --cluster --mds-plot 4 --out results/mds/PIB_BOA_SPI_NLS_HM3_mds
qsub bin/mds.job

# modify mds file to add our group as ethnicity in the first column
sed 's/^[\t ]\+//' results/mds/PIB_BOA_SPI_NLS_HM3_mds.mds | sed 's/^[0-9]\+/PBSN/' > results/mds/PIB_BOA_SPI_NLS_HM3_mds.withGroups.mds

### Make plot in R
library(calibrate)
mds.cluster = read.table("results/mds/PIB_BOA_SPI_NLS_HM3_mds.withGroups.mds", header=T)
# we assign colors to the groups
colors = c("darkolivegreen", "cyan", "brown", "orange", "black", "red", "pink", "grey", "darkblue", "purple", "green", "yellow")[mds.cluster$FID]
# open pdf to save the plot
pdf(file="results/mds/mdsplot6.pdf", width=7, height=7)
# make the plot
plot(mds.cluster$C2, mds.cluster$C1, col=colors, ylab="Dimension 1", xlab="Dimention 2", pch=21, main="PIB-BOA-SPI-NLS MDS")
# make the legend. make sure colors are in the same order as groups.
legend("bottomright", c("PIB-BOA-SPI-NLS", "ASW", "CEU", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX", "MKK", "TSI", "YRI"), fill=c("purple", "darkolivegreen", "cyan", "brown", "orange", "black", "red", "pink", "grey", "darkblue", "green", "yellow"), cex=.75)
# close the pdf
dev.off()

### We can see many outliers.

# We want to make the MDS plot with separate colors for each of our 4 groups
# For that we prepare a new file that contains the separate groups
head -1 results/mds/PIB_BOA_SPI_NLS_HM3_mds.withGroups.mds > results/mds/PIB_BOA_SPI_NLS_HM3_mds.separateGroups.txt
while read -r line; do
    if [[ $line == PBSN* ]]; then
        sample=$(echo $line | cut -f 2 -d " ");
        group=$(grep -w $sample data/IDs_dnaNum_projects_gender_birthdate.txt | cut -f 6);
        echo "$line" | sed "s/^PBSN/$group/";
    else echo "$line";
    fi ;
done < <(awk 'NR>1' results/mds/PIB_BOA_SPI_NLS_HM3_mds.withGroups.mds) | sort -k1,1 >> results/mds/PIB_BOA_SPI_NLS_HM3_mds.separateGroups.txt
# make plot in R
library(calibrate)
mds.cluster = read.table("results/mds/PIB_BOA_SPI_NLS_HM3_mds.separateGroups.txt", header=T)
# we assign colors to the groups
colors = c("#6d3dc5", "#66b645", "#cb4cc0", "#aba446", "#7577d3", "#ce8538", "#5f3771", "#64b391", "#d24731", "#5c90b6", "#cb4371", "#4b6331", "#c687b5", "#77382e", "#ce8a76")[mds.cluster$FID]
# open pdf to save the plot
pdf(file="results/mds/mdsplot7.pdf", width=7, height=7)
# make the plot
plot(mds.cluster$C2, mds.cluster$C1, col=colors, ylab="Dimension 1", xlab="Dimention 2", pch=21, main="PIB-BOA-SPI-NLS MDS")
# make the legend. make sure colors are in the same order as groups.
legend(
    "bottomright",
    c("ASW", "BOA", "CEU", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX", "MKK", "NLS", "PiB", "SPI", "TSI", "YRI"),
    fill=c("#6d3dc5", "#66b645", "#cb4cc0", "#aba446", "#7577d3", "#ce8538", "#5f3771", "#64b391", "#d24731", "#5c90b6", "#cb4371", "#4b6331", "#c687b5", "#77382e", "#ce8a76"),
    cex=.75)
# close the pdf
dev.off()
# Make the same plot but labelling some of the points with their corresponding group
pdf(file="results/mds/mdsplot8.pdf", width=7, height=7)
plot(mds.cluster$C2, mds.cluster$C1, col=colors, ylab="Dimension 1", xlab="Dimention 2", pch=21, main="PIB-BOA-SPI-NLS MDS")
legend("bottomright", c("ASW", "BOA", "CEU", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX", "MKK", "NLS", "PiB", "SPI", "TSI", "YRI"), fill=c("#6d3dc5", "#66b645", "#cb4cc0", "#aba446", "#7577d3", "#ce8538", "#5f3771", "#64b391", "#d24731", "#5c90b6", "#cb4371", "#4b6331", "#c687b5", "#77382e", "#ce8a76"), cex=.75)
# label SPI with C1 > 0
textxy(subset(mds.cluster, FID=="SPI" & C1 > 0, select=c(C2)), subset(mds.cluster, FID=="SPI" & C1 > 0, select=c(C1)), subset(mds.cluster, FID=="SPI" & C1 > 0, select=c(FID)), offset=1.1)
# label BOA with C1 > 0
textxy(subset(mds.cluster, FID=="BOA" & C1 > 0, select=c(C2)), subset(mds.cluster, FID=="BOA" & C1 > 0, select=c(C1)), subset(mds.cluster, FID=="BOA" & C1 > 0, select=c(FID)), offset=1.1)
# label NLS with C1 > 0
textxy(subset(mds.cluster, FID=="NLS" & C1 > 0, select=c(C2)), subset(mds.cluster, FID=="NLS" & C1 > 0, select=c(C1)), subset(mds.cluster, FID=="NLS" & C1 > 0, select=c(FID)), offset=1.1)
# label PiB with C1 > 0
textxy(subset(mds.cluster, FID=="PiB" & C1 > 0, select=c(C2)), subset(mds.cluster, FID=="PiB" & C1 > 0, select=c(C1)), subset(mds.cluster, FID=="PiB" & C1 > 0, select=c(FID)), offset=1.1)
dev.off()

### We decided to set the threshold at the lowest Dimension 1 value corresponding to the GIH group
grep "^GIH" results/mds/PIB_BOA_SPI_NLS_HM3_mds.separateGroups.txt | sort -k4,4n | head -1 | cut -f 4
# The threshold is then 0.00266741

# Check how many we lose for each group with this threshold
grep "^PiB\|^BOA\|^SPI\|^NLS" results/mds/PIB_BOA_SPI_NLS_HM3_mds.separateGroups.txt | awk '{if ($4 >= 0.00266741) print $1;}' | sort | uniq -c
# 7 NLS
# 9 PiB
# 25 SPI
# Keep a list of the names that we need to exclude
grep -wf <(grep "^PiB\|^BOA\|^SPI\|^NLS" results/mds/PIB_BOA_SPI_NLS_HM3_mds.separateGroups.txt | awk '{if ($4 >= 0.00266741) print $2;}') results/preqc/preqc_sample_rate.fam | cut -f -2 -d " " > results/mds/mds_exclude.txt

# Use plink to exclude those individuals
plink --noweb --bfile results/preqc/preqc_sample_rate --remove results/mds/mds_exclude.txt --make-bed --out results/mds/after_mds



### -----------
### IMPUTATION
### -----------

mkdir results/final_qc
## calculate marker missing data rate
plink --noweb --bfile results/mds/after_mds --missing --out results/final_qc/after_mds
#qsub bin/missing.job

## filter SNPs with call rate less than 90%
plink --noweb --bfile results/mds/after_mds --geno 0.1 --make-bed --out results/final_qc/step1
# 0 variants removed due to missing genotype data (--geno).

## more extensive qc and hwe filter
plink --noweb --bfile results/final_qc/step1 --hwe 0.000001 --geno 0.03 --mind 0.02 --make-bed --out results/final_qc/stap2
# 718 variants removed due to Hardy-Weinberg exact test.
# 0 people removed due to missing genotype data (--mind)
# 0 variants removed due to missing genotype data (--geno)

# first remove chromosome 25 and 26, because they are not imputed anyway. (25: pseudo-autosomal, 26: Mitochondrial)
grep "^25\|^26" results/final_qc/stap2.bim | cut -f 3 > results/final_qc/exclude_chr25-26.txt
plink --noweb --bfile results/final_qc/stap2 --exclude results/final_qc/exclude_chr25-26.txt --make-bed --out results/final_qc/step3

## Remove homozygosity outliers
# remove sex chromosome first (anderson et al 2010)
grep "^23\|^24" results/final_qc/step3.bim | cut -f 3 > results/final_qc/exclude_chr23-24.txt
plink --noweb --bfile results/final_qc/step3 --exclude results/final_qc/exclude_chr23-24.txt --make-bed --out results/final_qc/step4
        # # do the pruning
        # plink --noweb --bfile results/final_qc/step4 --indep 50 5 1.25 --out results/final_qc/het_pruning
        # plink --noweb --bfile results/final_qc/step4 --extract results/final_qc/het_pruning.prune.in --make-bed --out results/final_qc/step5
# Calculate fvalues
plink --noweb --bfile results/final_qc/step4 --het --out results/final_qc/step4_het
# With the commands saved in bin/imiss-vs-het.Rscript, we can plot heterozygosity against missingness and calculate heterozygosity outliers
# In the end we did in excel, calculate the mean heterozygosity, the average of meanHet and the standard deviation
# Then the threshold was set to the average(meanHet) +- 3 times the standard deviation
# Based on this, we exclude 34 individuals that we copied from the excel to results/final_qc/step4_het.exclude.txt
plink --noweb --bfile results/final_qc/step3 --remove results/final_qc/step4_het.exclude.txt --make-bed --out results/final_qc/step5

## remove SNPs with call rate below 98%
plink --noweb --bfile results/final_qc/step5 --geno 0.02 --make-bed --out results/final_qc/step6
# 171 variants removed due to missing genotype data (--geno).

## Now start Imputation steps from ENIGMA protocol
mkdir data/prepare_imp
mkdir results/prepare_imp

## make symlinks to
ln -s /home/multifac/Software_utilities/1000_genomes_imputatie/1kgp.* data/prepare_imp

## we need to couple the snps to the rs number
# copy rskey file
cp /home/multifac/BIG/PsychChip/raw_data/InfiniumPsychArray-24v1-1_A1_b138_rsids.txt data/prepare_imp
awk '{print $2, $1, $3, $4, $5, $6}' results/final_qc/step6.bim > results/prepare_imp/temp.bim
awk 'FNR==NR{a[$1]=$2;next} $1 in a{print a[$1],$2,$3,$4,$5,$6}' data/prepare_imp/InfiniumPsychArray-24v1-1_A1_b138_rsids.txt results/prepare_imp/temp.bim > results/prepare_imp/predatfile.bim
awk '{print $2,$1,$3,$4,$5,$6}' results/prepare_imp/predatfile.bim > results/prepare_imp/rsnumbers.bim
# PROBLEM: the rsnumbers.bim is not the same size as the original .bim, because not all of the variants could be linked
# SOLUTION: get the SNPs that could be linked from the original files and then couple the rsnumbers
awk '{print $2,$1}' data/prepare_imp/InfiniumPsychArray-24v1-1_A1_b138_rsids.txt > results/prepare_imp/switched_coupling.txt
awk 'FNR==NR{a[$1]=$2;next} $1 in a{print a[$1],$2,$3,$4,$5,$6}' results/prepare_imp/switched_coupling.txt results/prepare_imp/predatfile.bim > results/prepare_imp/switched_predatfile.bim
# keep only the marker id
awk '{print $1}' results/prepare_imp/switched_predatfile.bim > results/prepare_imp/PsychChip_met_rscodes.txt
# use the marker ids to extact from original data
plink --bfile results/final_qc/step6 --extract results/prepare_imp/PsychChip_met_rscodes.txt --make-bed --noweb --out results/prepare_imp/PiB-NLS-SPI-BOA_rsMatched
# recouple those that are left to the rsnumbers
awk '{print $2,$1,$3,$4,$5,$6}' results/prepare_imp/PiB-NLS-SPI-BOA_rsMatched.bim > results/prepare_imp/final_temp.bim
awk 'FNR==NR{a[$1]=$2;next} $1 in a{print a[$1],$2,$3,$4,$5,$6}' data/prepare_imp/InfiniumPsychArray-24v1-1_A1_b138_rsids.txt results/prepare_imp/final_temp.bim > results/prepare_imp/final_predatfile.bim
awk '{print $2,$1,$3,$4,$5,$6}' results/prepare_imp/final_predatfile.bim > results/prepare_imp/final_rsnumbers.bim

# Note: there is a dot (.) remaining in the marker IDs of results/prepare_imp/final_rsnumbers.bim
# Now results/prepare_bim/final_rsnumbers.bim has the same number of variants as results/prepare_imp/PiB-NLS-SPI-BOA_rsMatched.bim
# Rename the bim file to replace results/prepare_imp/neuroimage_psychchip24v1_rsMatched.bim
mv results/prepare_imp/PiB-NLS-SPI-BOA_rsMatched.bim results/prepare_imp/PiB-NLS-SPI-BOA_rsMatched.withoutRS.bim
mv results/prepare_imp/final_rsnumbers.bim results/prepare_imp/PiB-NLS-SPI-BOA_rsMatched.bim

### STEP2 ###
# find which SNPs are strand ambiguous
awk '{ if (($5=="T" && $6=="A") || ($5=="A" && $6=="T") || ($5=="C" && $6=="G") || ($5=="G" && $6=="C")) print $2}' results/prepare_imp/PiB-NLS-SPI-BOA_rsMatched.bim > results/prepare_imp/ambig.list

# use ambig.list to exclude
plink --noweb --bfile results/prepare_imp/PiB-NLS-SPI-BOA_rsMatched --exclude results/prepare_imp/ambig.list --make-founders --out results/prepare_imp/lastQC --maf 0.01 --hwe 0.000001 --make-bed
# --hwe: 4 variants removed due to Hardy-Weinberg exact test.
# 221168 variants removed due to minor allele threshold(s)

#### STEP 3 ####
awk '{print $2,$1,$3,$4,$5,$6}' results/prepare_imp/lastQC.bim > results/prepare_imp/tempQC.bim
# join the two files
awk 'NR==FNR{s=$1;a[s]=$0;next} a[$1]{print $0 " "a[$1]}' results/prepare_imp/tempQC.bim data/prepare_imp/1kgp.alleles > results/prepare_imp/merged.alleles
# select list of snps that need flipping
awk '{ if ($2!=$8 && $2!=$9) print $1}' results/prepare_imp/merged.alleles > results/prepare_imp/flip.list

# do the flip and update chrom
plink --bfile results/prepare_imp/lastQC --extract data/prepare_imp/1kgp.snps --update-map data/prepare_imp/1kgp.chr --update-chr --flip results/prepare_imp/flip.list --make-bed --noweb --out results/prepare_imp/step7
# update position
plink --bfile results/prepare_imp/step7 --update-map data/prepare_imp/1kgp.bp --geno 0.05 --mind 0.05 --make-bed --out results/prepare_imp/lastQC_b37 --noweb

# make lists of males and females
awk '{ if($5==1) print $1, $2}' results/prepare_imp/lastQC_b37.fam > results/prepare_imp/male.list
awk '{ if($5==2) print $1, $2}' results/prepare_imp/lastQC_b37.fam > results/prepare_imp/female.list

#### STEP 4 ####
# split the file per chromosome
mkdir data/mach results/mach
# create jobfile for each chromosome
bash bin/create_jobs_split_chr.sh
# submit jobs
for job in bin/split_chr*.job; do qsub $job; done
#
mv data/mach/ready4mach.23.female.map data/mach/ready4mach.23.map
rm data/mach/ready4mach.23.male.map

### Reformat dat
# create jobfiles to reformat each chrom's file
bash bin/create_jobs_reformat.sh
# from bin/reformat_chr23.job remove the last line because there is not such ped file
# submit the reformat jobs
for job in bin/reformat_chr*.job; do qsub $job; done

### chunking the chromosomes
cd data/mach
for i in {1..23}; do /home/multifac/bin/ChunkChromosome -d ready4mach."$i".dat.gz -n 5000 -o 500; done
cd ../..

### Mach phasing
# prepare the scripts
bash bin/create_jobs_phasing.sh
# NOTE!!: manually edit jobfile for chromosome 23 in order to use the correct ped file: ready4mach.23.female.ped as well as the output prefix
# submit
for job in bin/phasing_chr*.job; do qsub $job; done

### setting up the male X chromosome data for imputation
cat data/mach/ready4mach.23.male.ped | awk '{printf "%s", $1 "->"$2 " HAPLO1 "; for (N=7; N<=NF; N+=2) printf "%s", $N; printf "\n"; printf "%s", $1 "->"$2 " HAPLO2 "; for (N=7; N<=NF; N+=2) printf "%s", $N; printf "\n";}' > data/mach/ready4mach.23.male
gzip data/mach/ready4mach.23.male

### IMPUTATION
mkdir results/imputation

# fix the female reference vcf because it has spaces instead of tabs
mkdir data/reference_X
zgrep "^#" /home/multifac/Software_utilities/1000_genomes_imputatie/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.gz > data/reference_X/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.fixed.gz
zgrep -v "^#" /home/multifac/Software_utilities/1000_genomes_imputatie/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.gz | tr ' ' '\t' >> data/reference_X/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.fixed.gz

### At this point we have to move everything to the multifac account
cp -r /home/marc/workspace/imputation_PIB_BOA_SPI_NLS /home/multifac/temp
cd /home/multifac/temp/imputation_PIB_BOA_SPI_NLS

# create job files
bash bin/create_jobs_imputation.sh

# submit all jobs
for job in bin/imputation_chr*.job; do
    qsub $job
done

## When the job is finished, we can check the logs. In case of succesful completion it should end with "Run Completed in ..."
## To check all of the logs:
tail -n 3 results/imputation/*.log | less
## Some of the jobs failed, probably because I exceeded my disk quota as mpauper
## after cleaning some of my files, I manually created a jobfile to resubmit only those chunks and chromosomes that had failed
qsub bin/resubmit_imputation.job
## NOTE: this may not be necessary in other projects.

# All done

# clean up draft files
rm results/imputation/*.draft
