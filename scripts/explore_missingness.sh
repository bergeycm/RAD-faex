#!/bin/sh

# Call from NGS-map folder
cd /scratch/cmb433/fecalRAD/NGS-map/

module load plink
module load r

INPUT_PREFIX=baboon.pass.snp

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Missingness and filtration based on missingness
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Get stats on missingness
# -------------------------------------------------------------------------------------- #

# === Getting stats on missingness ============================================
plink --noweb --file ${INPUT_PREFIX} --missing --out reports/${INPUT_PREFIX}

# -------------------------------------------------------------------------------------- #
# --- Remove SNPs with too much missing data
# -------------------------------------------------------------------------------------- #

# --geno 0.1 - only include SNPs that are genotyped in at least 90% of individuals
# --mind 0.9 - only include individuals that have genotypes for at least 10% of SNPs

# === Removing SNPs with too much missing data ================================ #";
plink --noweb --file ${INPUT_PREFIX} \
	--geno 0.1 --mind 0.9 \
	--missing \
	--make-bed --out results/all.cleaned
cp results/all.cleaned.log reports/plink_missing_exclusion_main.log

# Load file to get final genotyping rate
plink --noweb --bfile results/all.cleaned --recode --out results/all.cleaned
cp results/all.cleaned.log reports/plink_info_full.log

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Multidimensional scaling plot
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Find SNPs in LD to prune out
# -------------------------------------------------------------------------------------- #

# Slide a sliding window of 50 SNPs across each chromosome and calculate r^2 for each 
# pair in the window. Remove one random SNP in each pair with r^2 > 0.5.
# Parameters are window size, step, and r^2 threshold.

# === Finding SNPs in LD to prune out =========================================
plink --noweb --bfile results/all.cleaned \
	--indep-pairwise 50 5 0.5 \
	--out results/all.cleaned
cp results/all.cleaned.log reports/plink_LD_pruning_part1.log

# -------------------------------------------------------------------------------------- #
# --- Prune dataset for LD
# -------------------------------------------------------------------------------------- #

# === Pruning dataset for LD ================================================== #
plink --noweb --bfile results/all.cleaned \
	--exclude results/all.cleaned.prune.out \
	--make-bed --out results/all.cleaned.LDpruned
cp results/all.cleaned.LDpruned.log reports/plink_LD_pruning_part2.log

# -------------------------------------------------------------------------------------- #
# --- Make *.genome file
# -------------------------------------------------------------------------------------- #

# === Making genome file ======================================================
plink --noweb --bfile results/all.cleaned.LDpruned \
	--genome \
	--out results/all.cleaned.LDpruned
cp results/all.cleaned.LDpruned.log reports/plink_genome_file.log

# -------------------------------------------------------------------------------------- #
# --- Do multidimensional scaling of similarities (proportions of alleles IBS)
# -------------------------------------------------------------------------------------- #

# === Doing multidimensional scaling of similarities ========================== #";
plink --noweb --bfile results/all.cleaned.LDpruned \
	--read-genome results/all.cleaned.LDpruned.genome \
	--cluster --mds-plot 2 \
	--out results/all.cleaned.LDpruned
cp results/all.cleaned.LDpruned.log reports/plink_mds.log

# -------------------------------------------------------------------------------------- #
# --- Cluster by missingness
# -------------------------------------------------------------------------------------- #

plink --noweb --bfile results/all.cleaned.LDpruned \
	--read-genome results/all.cleaned.LDpruned.genome \
	--cluster missing --mds-plot 2 \
	--out results/all.cleaned.LDpruned.missing

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Association test and missingness chi-sq test
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# ----------------------------------------------------------------------------------------
# --- Do an association test
# ----------------------------------------------------------------------------------------

grep "Tx" data/fecalRAD_individual_info.csv | \
	awk 'BEGIN { FS="," } { print $1,$1,$8 }' | \
	sed '1d' | sed -e "s/\s/.PE\t/g" > data/blood_feces_categories.txt

#	CHR     Chromosome
#	SNP     SNP ID
#	BP      Physical position (base-pair)
#	A1      Minor allele name (based on whole sample)
#	F_A     Frequency of this allele in cases
#	F_U     Frequency of this allele in controls
#	A2      Major allele name
#	CHISQ   Basic allelic test chi-square (1df)
#	P       Asymptotic p-value for this test
#	OR      Estimated odds ratio (for A1, i.e. A2 is reference)

# LD and missingness pruned dataset
plink --noweb --bfile results/all.cleaned.LDpruned \
	--make-pheno data/blood_feces_categories.txt feces \
	--assoc --allow-no-sex --out results/all.cleaned.LDpruned

# All SNPs
plink --noweb --file ${INPUT_PREFIX} \
	--make-pheno data/blood_feces_categories.txt feces \
	--assoc --allow-no-sex --out results/all

# ----------------------------------------------------------------------------------------
# --- Do a missing chi-sq test (Does SNP's missingness differ between cases and controls?)
# ----------------------------------------------------------------------------------------

# All SNPs
plink --noweb --file ${INPUT_PREFIX} \
	--make-pheno data/blood_feces_categories.txt feces \
	--test-missing \
	--allow-no-sex \
	--out results/all

# LD and missingness pruned dataset
plink --noweb --bfile results/all.cleaned.LDpruned \
	--make-pheno data/blood_feces_categories.txt feces \
	--test-missing \
	--allow-no-sex \
	--out results/all.cleaned.LDpruned

# ----------------------------------------------------------------------------------------
# --- Compute distance between samples
# ----------------------------------------------------------------------------------------

# All SNPs
plink --noweb --file ${INPUT_PREFIX} \
	--distance square flat-missing \
	--out results/all

# LD and missingness pruned dataset
plink --noweb --bfile results/all.cleaned.LDpruned \
	--distance square flat-missing \
	--out results/all.cleaned.LDpruned
