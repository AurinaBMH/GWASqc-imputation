
#       #___________________________________________________________________#
#       #\                                                                   \
#       # #-------------------------------------------------------------------#
#       # #                                                                   # 
#       # #                        QUALITY CONTROL                            #
#        \#                                                                   # 
#         #-------------------------------------------------------------------#

### unzip data
unzip data/genderfixes.zip
rm data/genderfixes.zip

### re-check the gender from core_data
qsub bin/check_core_gender.sh
# they all have good sex score (F < 0.2 for female, F > 0.8 for male)
# The idea is to keep all those that Ziarth says "missing" and exlude those that Ziarth insists on opposite gender

### remove those with gender discrepancy, but keep those with missing gender (ziarth)
mkdir results/gender_fix
qsub bin/remove_gender_discr.sh
# removed 45 individuals
# output is: results/gender_fix/core_data_1

### check relatedness
mkdir results/relatedness
# make list of snps not on autosomal chromosomes
grep "^23\|^24\|^26" results/gender_fix/core_data_1.bim | cut -f 2 > results/relatedness/chr23_26.txt
qsub bin/indep.sh
qsub bin/relatedness.sh ## Did not complete...




#       *___________________________________________________________________#
#       |\                                                                   \
#       | *-------------------------------------------------------------------*
#       | |                                                                   |
#       | |                        MULTIDIMENSIONAL SCALING                   |
#        \|                                                                   | 
#         *-------------------------------------------------------------------*

### Subset variants
qsub bin/filter_preMDS.sh
# 51 markers to be excluded based on HWE test ( p <= 1e-06 )
#        26 markers failed HWE test in cases
#        51 markers failed HWE test in controls
# 4242 SNPs failed missingness test ( GENO > 0.05 )
# 1559 SNPs failed frequency test ( MAF < 0.01 )
# After frequency and genotyping pruning, there are 304518 SNPs

### make links to HapMap3 data
mkdir data/HM3
ln -s /home/multifac/Software_utilities/HapMap/HM3.bed data/HM3/
ln -s /home/multifac/Software_utilities/HapMap/HM3.bim data/HM3/
ln -s /home/multifac/Software_utilities/HapMap/HM3.fam data/HM3/

### From HM3 only keep SNPs in our dataset.
cut -f 2 results/mds/substep_variants.bim > data/HM3/HM3.snplist.txt
qsub bin/subset_HM3.sh

### Prune the SNPs
# plink --noweb --bfile data/HM3/HM3_subset --indep 50 5 1.15 --out data/HM3/HM3_subset
qsub bin/prune_HM3_job1.sh

### Exclude the prunned SNPs from hapmap data
# plink --noweb --bfile data/HM3/HM3_subset --extract data/HM3/HM3_subset.prune.in --make-bed --out data/HM3/HM3_subset_pruned
qsub bin/prune_HM3_job2.sh

### Extract prunned SNPs from our plink files
# plink --noweb --bfile results/mds/substep_variants --extract data/HM3/HM3_subset.prune.in --make-bed --out results/mds/subset_variants_pruned
qsub bin/prune_AU_job1.sh

### Create the files for updating the positions in order that both files contain the same and one file (forced_allele_from_hapmap_exomeseq) for updating the reference allele for the same reason.
awk '{print($2"\t"$1)}' data/HM3/HM3_subset_pruned.bim > results/mds/hm3_CHR_update
awk '{print($2"\t"$4)}' data/HM3/HM3_subset_pruned.bim > results/mds/hm3_BP_update
awk '{print($2"\t"$5)}' data/HM3/HM3_subset_pruned.bim > results/mds/hm3_allele_update

### Now we will update the position of our data considering the data from hapmap
# plink --noweb --bfile results/mds/subset_variants_pruned --update-map results/mds/hm3_CHR_update --update-chr --make-bed --out results/mds/subset_variants_pruned_chr
qsub bin/update_chr_mds.sh
# plink --noweb --bfile results/mds/subset_variants_pruned_chr --update-map results/mds/hm3_BP_update --make-bed --out results/mds/subset_variants_pruned_pos
qsub bin/update_pos_mds.sh

### Force the reference allele in our data using the data from hapmap
# plink --noweb --bfile results/mds/subset_variants_pruned_pos --reference-allele results/mds/hm3_allele_update --make-bed --out results/mds/subset_variants_pruned_forced
qsub bin/force_allele_mds.sh
# 20434 could not be forced
# checked a few manually and it seems that it is simply wrong strand, so will try to flip the strand for those

## find which ones have the problem
awk 'NR==FNR{a[$2]=$5;next}{ if ($2 in a) print $0"\t"a[$2]; else print $0"\t.";}' results/mds/subset_variants_pruned_forced.bim data/HM3/HM3_subset_pruned.bim | awk '{if ($5 != $7) print $2;}' > results/mds/flip_snps.txt
## use this list to flip
plink --noweb --bfile results/mds/subset_variants_pruned_pos --flip results/mds/flip_snps.txt --make-bed --out results/mds/subset_variants_pruned_flip

### Re-try to force the alleles now that we flipped
# plink --noweb --bfile results/mds/subset_variants_pruned_flip --reference-allele results/mds/hm3_allele_update --make-bed --out results/mds/subset_variants_pruned_forced_2
qsub bin/force_allele_mds.2.sh
# Set reference alleles for 40807 SNPs, 2952 different from minor allele
# PROBLEM FIXED

### Merging our files with hapmap files
# plink --noweb --bfile results/mds/subset_variants_pruned_flip --bmerge data/HM3/HM3_subset_pruned.bed data/HM3/HM3_subset_pruned.bim data/HM3/HM3_subset_pruned.fam --make-bed --out results/mds/Bellgrove_HM3_merged
qsub bin/merge_bell_hm3.sh
# DONE!!!! No warnings were obtained!

### Make genome file
# plink --noweb --bfile results/mds/Bellgrove_HM3_merged --genome --out results/mds/Bellgrove_HM3_merged
qsub bin/mkgenome.sh

### MDS
# plink --noweb --bfile results/mds/Bellgrove_HM3_merged --read-genome results/mds/Bellgrove_HM3_merged.genome  --mind .05 --cluster --mds-plot 4 --out results/mds/Bellgrove_HM3_MDS
qsub bin/multidim.sh

# make plot in R
library(calibrate)
mds.cluster = read.table("results/mds/Bellgrove_HM3_MDS.mds", header=T)
ur.num = length(mds.cluster$C1) - 988
colors = rev(c(rep("red", ur.num), rep("lightblue", 112), rep("brown", 84), rep("yellow", 113), rep("green", 88), rep("purple", 86), rep("orange", 85), rep("grey50", 50), rep("black", 88), rep("darkolivegreen", 49), rep("magenta", 90), rep("darkblue", 143)))
pdf(file="results/mds/mdsplot3.pdf")
plot(rev(mds.cluster$C2), rev(mds.cluster$C1), col=colors, ylab="Dimension 1", xlab="Dimention 2", pch=21, main="Bellgrove MDS")
legend("bottomright", c("Bellgrove", "CEU", "CHB", "YRI", "TSI", "JPT", "CHD", "MEX", "GIH", "ASW", "LWK", "MKK"), fill=c("red", "lightblue", "brown", "yellow", "green", "purple", "orange", "grey50", "black", "darkolivegreen", "magenta", "darkblue"))
#textxy(mds.cluster$C2[1:ur.num], mds.cluster$C1[1:ur.num], mds.cluster$IID[1:ur.num])  ## This is used to add the labels to identify data points
dev.off()

# 7 samples seem to be outliers so we will exclude them. Those are:
# 570044
# 5615253
# 690095
# 690213
# 5607583
# 570081
# 700026

# exclude these samples
grep -v "570044\|5615253\|690095\|690213\|5607583\|570081\|700026" results/mds/Bellgrove_HM3_MDS.mds > results/mds/Bellgrove_HM3_MDS_exclude.mds

# make new plot in R without excluded samples
library(calibrate)
mds.cluster = read.table("results/mds/Bellgrove_HM3_MDS_exclude.mds", header=T)
ur.num = length(mds.cluster$C1) - 988
colors = rev(c(rep("red", ur.num), rep("lightblue", 112), rep("brown", 84), rep("yellow", 113), rep("green", 88), rep("purple", 86), rep("orange", 85), rep("grey50", 50), rep("black", 88), rep("darkolivegreen", 49), rep("magenta", 90), rep("darkblue", 143)))
pdf(file="results/mds/mdsplot4.pdf")
plot(rev(mds.cluster$C2), rev(mds.cluster$C1), col=colors, ylab="Dimension 1", xlab="Dimention 2", pch=21, main="Bellgrove MDS")
legend("bottomright", c("Bellgrove", "CEU", "CHB", "YRI", "TSI", "JPT", "CHD", "MEX", "GIH", "ASW", "LWK", "MKK"), fill=c("red", "lightblue", "brown", "yellow", "green", "purple", "orange", "grey50", "black", "darkolivegreen", "magenta", "darkblue"))
#textxy(mds.cluster$C2[1:ur.num], mds.cluster$C1[1:ur.num], mds.cluster$IID[1:ur.num])  ## This is used to add the labels to identify data points
dev.off()




#       #___________________________________________________________________#
#       #\                                                                   \
#       # #-------------------------------------------------------------------#
#       # #                                                                   # 
#       # #                             IMPUTATION                            #
#        \#                                                                   # 
#         #-------------------------------------------------------------------#

### NOTE: This part was not completed with the 1000G phase 3 reference panel because of space and walltime problems.
###       Scroll down to "IMPUTATION (OLD Reference panel)" for the part that was used afterall

####### QUALITY CONTROL #######

# as the superclean dataset provided by Rust contains 310347 snps, which is the same as the core_data
# I am not sure that the snps have been removed
# I will perform the quality steps regardings markers again

mkdir results/qc_imputation

# first we need to remove the 7 samples that we exclude based on the MDS analysis
grep "570044\|5615253\|690095\|690213\|5607583\|570081\|700026" data/genderfixes/superclean_232true/superclean-232true.fam > results/qc_imputation/mds_exclude.txt
plink --noweb --bfile data/genderfixes/superclean_232true/superclean-232true --remove results/qc_imputation/mds_exclude.txt --make-bed --out results/qc_imputation/superclean_step1
# qsub bin/qc_excl_mds.sh

# calculate marker missing data rate
plink --noweb --bfile results/qc_imputation/superclean_step1 --missing --out results/qc_imputation/superclean_step1
# qsub bin/lmiss.sh

# test for different call rate in cases and controls
plink --noweb --bfile results/qc_imputation/superclean_step1 --test-missing --out results/qc_imputation/superclean_step1
# qsub bin/tmiss.sh

# make plot of lmiss
Rscript bin/lmiss-hist.R  # this script comes from Anderson et al. 2010

# make list of markers with significantly different missingness in cases versus controls
perl bin/run-diffmiss-qc.pl results/qc_imputation/superclean_step1  # this script comes from Anderson et al. 2010 

# based on the above steps, exclude markers
plink --noweb --bfile results/qc_imputation/superclean_step1 --exclude results/qc_imputation/fail-diffmiss-qc.txt --make-bed --out results/qc_imputation/superclean_step2
plink --noweb --bfile results/qc_imputation/superclean_step2 –-maf 0.01 --geno 0.05 --hwe 0.00001 --make-bed --out results/qc_imputation/superclean_step3
# qsub bin/qc_excl_markers.sh
# qsub bin/qc_filter_markers.sh

####### PREPARATION #######
mkdir data/preparation
mkdir results/preparation

### make file with all chr:position and alleles of the reference panel
qsub bin/make_alleles_file.job

### find which SNPs are strand ambiguous
awk '{ if (($5=="T" && $6=="A") || ($5=="A" && $6=="T") || ($5=="C" && $6=="G") || ($5=="G" && $6=="C")) print $2, "ambig"; else print $2; }' results/qc_imputation/superclean_step3.bim | grep ambig | awk '{print $1}' > results/preparation/ambig.list
## use ambig.list to exclude
plink --noweb --bfile results/qc_imputation/superclean_step3 --exclude results/preparation/ambig.list --make-founders --out results/preparation/superclean_step4 --maf 0.01 --hwe 0.000001 --make-bed
# qsub bin/excl_ambig.job

### find common markers in our dataset and the reference panel and compare
awk '{print $1":"$4"\t"$5"\t"$6"\t"$2}' results/preparation/superclean_step4.bim > results/preparation/superclean_step4.alleles
# join the two files
awk 'NR==FNR{s=$1;a[s]=$0;next} {if ($1 in a) print $0"\t"a[$1]}' results/preparation/superclean_step4.alleles data/preparation/1kgp.phase3.alleles > results/preparation/step4_1kgphase3.alleles
# this creates some duplicate ones, we will keep only one occurence
while read -r line; do
    grep -w -m 1 "$line" results/preparation/step4_1kgphase3.alleles;
done < <(cut -f 1 results/preparation/step4_1kgphase3.alleles | sort | uniq -c | sed 's/^[ \t]*//' | grep "^1 " | cut -f 2 -d " ") \
> results/preparation/step4_1kgphase3.unique.alleles

# compare alleles
awk '{ if ($2==$5 && $3==$6) print $7"\tcorrect"; else if ($2==$6 && $3==$5) print $7"\tforceRef"; else print $7"\terror"}' results/preparation/step4_1kgphase3.unique.alleles > results/preparation/step4_1kgphase3.comparisons.txt

### those that have "error" we will try to flip and keep only those in common
cut -f 7 results/preparation/step4_1kgphase3.unique.alleles > results/preparation/common.snps
grep "error" results/preparation/step4_1kgphase3.comparisons.txt | cut -f 1 > results/preparation/flip.list
plink --bfile results/preparation/superclean_step4 --extract results/preparation/common.snps --flip results/preparation/flip.list --make-bed --noweb --out results/preparation/superclean_step4_flipped
# qsub bin/flip_extract.job

### remake the comparison to see how many we fixed
awk '{print $1":"$4"\t"$5"\t"$6"\t"$2}' results/preparation/superclean_step4_flipped.bim > results/preparation/superclean_step4_flipped.alleles
# join the two files
awk 'NR==FNR{s=$1;a[s]=$0;next} {if ($1 in a) print $0"\t"a[$1]}' results/preparation/superclean_step4_flipped.alleles data/preparation/1kgp.phase3.alleles > results/preparation/step4_flipped_1kgphase3.alleles
# compare alleles
awk '{ if ($2==$5 && $3==$6) print $7"\tcorrect"; else if ($2==$6 && $3==$5) print $7"\tforceRef"; else print $7"\terror"}' results/preparation/step4_flipped_1kgphase3.alleles > results/preparation/step4_flipped_1kgphase3.comparisons.txt

## Now only 54 still have an error. We will exclude them
grep "error" results/preparation/step4_flipped_1kgphase3.comparisons.txt | cut -f 1 > results/preparation/error_markers.txt
plink --noweb --bfile results/preparation/superclean_step4_flipped --exclude results/preparation/error_markers.txt --make-bed --out results/preparation/superclean_step5 
# qsub bin/excl_error.job

### In the previous step we lost the X chromosome
awk 'NR==FNR{s=$1;a[s]=$0;next} {if ($1 in a) print $0"\t"a[$1]}' results/preparation/superclean_step4.chr23.alleles data/reference_X/1kgp.23.alleles > results/preparation/step4_1kgphase3.chr23.alleles
# compare alleles
awk '{ if ($2==$5 && $3==$6) print $7"\tcorrect"; else if ($2==$6 && $3==$5) print $7"\tforceRef"; else print $7"\terror"}' results/preparation/step4_1kgphase3.chr23.alleles > results/preparation/step4_1kgphase3.chr23.comparisons.txt
# those that have "error" we will try to flip and keep only those in common
cut -f 7 results/preparation/step4_1kgphase3.chr23.alleles > results/preparation/common.23.snps
grep "error" results/preparation/step4_1kgphase3.chr23.comparisons.txt | cut -f 1 > results/preparation/flip.23.list
plink --bfile results/preparation/superclean_step4 --extract results/preparation/common.23.snps --flip results/preparation/flip.23.list --make-bed --noweb --out results/preparation/superclean_step4_flipped_chr23
# recheck how many errors we have now
awk '{print $1":"$4"\t"$5"\t"$6"\t"$2}' results/preparation/superclean_step4_flipped_chr23.bim > results/preparation/superclean_step4_flipped.23.alleles
awk 'NR==FNR{s=$1;a[s]=$0;next} {if ($1 in a) print $0"\t"a[$1]}' results/preparation/superclean_step4_flipped.23.alleles data/reference_X/1kgp.23.alleles > results/preparation/step4_flipped_1kgphase3.23.alleles
awk '{ if ($2==$5 && $3==$6) print $7"\tcorrect"; else if ($2==$6 && $3==$5) print $7"\tforceRef"; else print $7"\terror"}' results/preparation/step4_flipped_1kgphase3.23.alleles > results/preparation/step4_flipped_1kgphase3.chr23.comparisons.txt
# now 62 have an error we will remove these
grep  "error" results/preparation/step4_flipped_1kgphase3.chr23.comparisons.txt | cut -f 1 > results/preparation/error_markers.23.txt
plink --noweb --bfile results/preparation/superclean_step4_flipped_chr23 --exclude results/preparation/error_markers.23.txt --make-bed --out results/preparation/superclean_step5_23

### make lists of males and females
awk '{ if($5==1) print $1, $2}' results/preparation/superclean_step5.fam > results/preparation/male.list
awk '{ if($5==2) print $1, $2}' results/preparation/superclean_step5.fam > results/preparation/female.list


####### MACH PHASING #######
mkdir data/Mach
mkdir results/phasing

### Split per chromosome
# create jobfiles
bash bin/create_jobs_split_chr.sh
# submit jobs
for job in bin/split_chr*.job; do qsub $job; done

#
mv data/Mach/ready4mach.23.female.map data/Mach/ready4mach.23.map
rm data/Mach/ready4mach.23.male.map

### Reformat
# create jobfiles to reformat each chrom's file
bash bin/create_jobs_reformat.sh
# submit all the jobs
for job in bin/reformat_chr*.job; do qsub $job; done

### chunking the chromosomes and phasing
cd data/Mach  # must cd there for the next step to work
for i in {1..23}; do ChunkChromosome -d ready4mach."$i".dat.gz -n 5000 -o 500; done
cd ../..

### MaCH phasing
# prepare the scripts
bash bin/create_jobs_phasing.sh
# NOTE: manually edit jobfile for chromosome 23 in order to use the correct ped file: ready4mach.23.female.ped as well as the output prefix
# submit
for job in bin/phasing_chr*.job; do qsub $job; done

### setting up the male X chromosome data for imputation
cat data/Mach/ready4mach.23.male.ped | awk '{printf "%s", $1 "->"$2 " HAPLO1 "; for (N=7; N<=NF; N+=2) printf "%s", $N; printf "\n"; printf "%s", $1 "->"$2 " HAPLO2 "; for (N=7; N<=NF; N+=2) printf "%s", $N; printf "\n";}' > data/Mach/ready4mach.23.male
gzip data/Mach/ready4mach.23.male

####### IMPUTATION #######
mkdir results/imputation

# create job files
bash bin/create_jobs_imputation.sh

# submit all jobs
for job in bin/imputation_chr*.job; do
	qsub $job
done


#       #___________________________________________________________________#
#       #\                                                                   \
#       # #-------------------------------------------------------------------#
#       # #                                                                   # 
#       # #                  IMPUTATION (OLD Reference panel)                 #
#        \#                                                                   # 
#         #-------------------------------------------------------------------#


####### PREPARATION #######
mkdir data/preparation_2
mkdir results/preparation_2

### find common markers in our dataset and the reference panel and compare
ln -s /home/multifac/Software_utilities/1000_genomes_imputatie/1kgp.* data/preparation_2/

## we need to couple the snps to the rs number
# copy rskey file
cp /home/multifac/BIG/PsychChip/raw_data/InfiniumPsychArray-24v1-1_A1_b138_rsids.txt data/preparation_2/
awk '{print $2, $1, $3, $4, $5, $6}' results/preparation/superclean_step4.bim > results/preparation_2/temp.bim
awk 'FNR==NR{a[$1]=$2;next} $1 in a{print a[$1],$2,$3,$4,$5,$6}' data/preparation_2/InfiniumPsychArray-24v1-1_A1_b138_rsids.txt results/preparation_2/temp.bim > results/preparation_2/predatfile.bim
awk '{print $2,$1,$3,$4,$5,$6}' results/preparation_2/predatfile.bim > results/preparation_2/rsnumbers.bim
# rename and place with the other corresponding plink files
mv results/preparation/superclean_step4.bim results/preparation/superclean_step4_withoutRS.bim
mv results/preparation_2/rsnumbers.bim results/preparation/superclean_step4.bim

# join the two files
awk '{print $2,$1,$3,$4,$5,$6}' results/preparation/superclean_step4.bim > results/preparation_2/tempQC.bim
awk 'NR==FNR{s=$1;a[s]=$0;next} a[$1]{print $0 " "a[$1]}' results/preparation_2/tempQC.bim data/preparation_2/1kgp.alleles > results/preparation_2/merged.alleles
# select list of snps that need flipping
awk '{ if ($2!=$8 && $2!=$9) print $1}' results/preparation_2/merged.alleles > results/preparation_2/flip.list
#!!! # compare alleles
#!!! awk '{ if ($2==$8 && $3==$9) print $0"\tcorrect"; else if ($2==$9 && $3==$8) print $0"\tforceRef"; else print $0"\terror"}' results/preparation_2/merged.alleles > results/preparation_2/merged.comparisons.txt

# do the flip and update chrom
plink --bfile results/preparation/superclean_step4 --extract data/preparation_2/1kgp.snps --update-map data/preparation_2/1kgp.chr --update-chr --flip results/preparation_2/flip.list --make-bed --noweb --out results/preparation_2/superclean_step4.5
# here we lost some markers because there were many that were duplicate

# update position
plink --bfile results/preparation_2/superclean_step4.5 --update-map data/preparation_2/1kgp.bp --geno 0.05 --mind 0.05 --make-bed --out results/preparation_2/superclean_step5 --noweb

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
# awk '{ if($5==1) print $1, $2}' results/preparation/superclean_step5.fam > results/preparation/male.list
# awk '{ if($5==2) print $1, $2}' results/preparation/superclean_step5.fam > results/preparation/female.list


####### MACH PHASING #######
mkdir data/Mach_2
mkdir results/phasing_2
mkdir bin_2

### Split per chromosome
# create jobfiles
bash bin_2/create_jobs_split_chr.sh
# submit jobs
for job in bin_2/split_chr*.job; do qsub $job; done

#
mv data/Mach_2/ready4mach.23.female.map data/Mach_2/ready4mach.23.map
rm data/Mach_2/ready4mach.23.male.map

### Reformat dat
# create jobfiles to reformat each chrom's file
bash bin_2/create_jobs_reformat.sh
# submit all the jobs
for job in bin_2/reformat_chr*.job; do qsub $job; done

### chunking the chromosomes and phasing
cd data/Mach_2  # must cd there for the next step to work
for i in {1..23}; do ChunkChromosome -d ready4mach."$i".dat.gz -n 5000 -o 500; done
cd ../..

### MaCH phasing
# prepare the scripts
bash bin_2/create_jobs_phasing.sh
# NOTE!!: manually edit jobfile for chromosome 23 in order to use the correct ped file: ready4mach.23.female.ped as well as the output prefix
# submit
for job in bin_2/phasing_chr*.job; do qsub $job; done

### setting up the male X chromosome data for imputation
zcat data/Mach_2/ready4mach.23.male.ped.gz | awk '{printf "%s", $1 "->"$2 " HAPLO1 "; for (N=7; N<=NF; N+=2) printf "%s", $N; printf "\n"; printf "%s", $1 "->"$2 " HAPLO2 "; for (N=7; N<=NF; N+=2) printf "%s", $N; printf "\n";}' > data/Mach_2/ready4mach.23.male
gzip data/Mach_2/ready4mach.23.male

####### IMPUTATION #######
#mkdir results/imputation

# fix the female reference vcf because it has spaces instead of tabs
zgrep "^#" /home/multifac/Software_utilities/1000_genomes_imputatie/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.gz > data/reference_X/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.fixed.gz
zgrep -v "^#" /home/multifac/Software_utilities/1000_genomes_imputatie/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.gz | tr ' ' '\t' >> data/reference_X/chr23.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.nosingles.vcf.fixed.gz

# create job files
bash bin_2/create_jobs_imputation.sh

# submit all jobs
for job in bin_2/imputation_chr*.job; do
    qsub $job
done

### In the end, some of the jobs were run on the multifac accout as they could not be completed on my own account. (Don't know why)

### Also, switching to Minimac from Minimac2, Minimac2 raised errors regarding autoChunk files.

rm results/imputation/*draft
