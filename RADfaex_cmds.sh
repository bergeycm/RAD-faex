#!/bin/bash

module load r

# ----------------------------------------------------------------------------------------
# --- Simulate ddRAD digestion -----------------------------------------------------------
# ----------------------------------------------------------------------------------------

# Called with:
#     qsub pbs/digest_sim.pbs
# Generates results/chr*.bed

Rscript scripts/ddRAD_sim.R $THIS_CHR_FA

# ----------------------------------------------------------------------------------------
# --- Combine BED files of RADtags
# ----------------------------------------------------------------------------------------

# Generates results/RADtags.bed

sh scripts/combine_RADtag_beds.sh

# ----------------------------------------------------------------------------------------
# --- Compute GC within all RADtags (from simulation) as well as in the larger region
# ----------------------------------------------------------------------------------------

# Called with:
#     pbs/compute_GC.pbs
# Generates results/RADtags.gc.bed, 
#           results/RADtags.slop${REGION_EXPANSION}.bed
#           results/RADtags.slop${REGION_EXPANSION}.gc.bed

sh scripts/get_RADtag_GC.sh

# ----------------------------------------------------------------------------------------
# --- Download CpG islands
# ----------------------------------------------------------------------------------------

# Generates data/cpgIslandExtUnmasked.chr.txt
#           data/cpgIslandExtUnmasked.chr.bed

sh scripts/download_CpG_islands.sh

# ----------------------------------------------------------------------------------------
# --- Find CpGs
# ----------------------------------------------------------------------------------------

# Called with:
#     qsub pbs/call_CpG_finder.pbs
# Generates results/papAnu2.CpG.bed

perl scripts/find_CpGs.pl > results/papAnu2.CpG.bed

# ----------------------------------------------------------------------------------------
# --- Find closest CpG
# ----------------------------------------------------------------------------------------

# Called with:
#     qsub pbs/call_closest_CpG_finder.pbs 
# Generates results/RADtags.nearestCpGisland.bed

sh find_closest_CpG.sh

# ========================================================================================
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Compute coverage within all RADtags (from simulation)
# ----------------------------------------------------------------------------------------

# Called with:
#     qsub -t 0-33 pbs/get_RAD_cov.pbs
# Generates results/${ID}.gt1.cov.bed

sh scripts/get_RADtag_cov.sh

# ----------------------------------------------------------------------------------------
# --- Combine RADtag coverage from all individuals
# ----------------------------------------------------------------------------------------

# Generates results/RADtag.coverage.all.txt

sh scripts/combine_coverage.sh

# ----------------------------------------------------------------------------------------
# --- Plot heatmap of coverage and RADtag length
# ----------------------------------------------------------------------------------------

# Generates reports/fecalRAD-BC*-BC*.gt1.cov.pdf

Rscript scripts/coverage_length_heatmap.R

# ----------------------------------------------------------------------------------------
# --- Combine RADtag info
# ----------------------------------------------------------------------------------------

# Called with:
#     qsub pbs/combine_RAD_info.pbs
# Generates results/RADtags.info.bed

perl scripts/combine_RADtag_info.pl > results/RADtags.info.bed

# ----------------------------------------------------------------------------------------
# --- Explore how coverage varies by sample, sample type, mapped read count, etc.
# ----------------------------------------------------------------------------------------

# Generates reports/RADtag_count_gtNreads.pdf
#           reports/RADtag_count_curves_by_total_mapped.pdf

Rscript scripts/plot_coverage.R

# ========================================================================================
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Model RADtag coverage (to explain inter-RADtag variation)
# ----------------------------------------------------------------------------------------

# See inter.RADtag.regression.R

qsub pbs/model_RADtag_cvg_1000.pbs
qsub pbs/model_RADtag_cvg_5000.pbs
qsub pbs/model_RADtag_cvg_10000.pbs
qsub pbs/model_RADtag_cvg_20000.pbs

# ========================================================================================
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Compute and explore heterozygosity
# ----------------------------------------------------------------------------------------

sh scripts/compute_het.sh
Rscript scripts/parse_heterozygosity.R

# ----------------------------------------------------------------------------------------
# --- Quantify allelic dropout by comparing blood and fecal DNA from single individual
# ----------------------------------------------------------------------------------------

perl scripts/quantify_ADO.pl 2> reports/quantify_ADO_stderr.txt

sh scripts/parse_all_discordance_matrices.sh > results/discordance.txt

Rscript scripts/explore_discordance.R

# ----------------------------------------------------------------------------------------
# --- Compute stats on missingness
# ----------------------------------------------------------------------------------------

# This uses SNPs called within one individual, not the results of multi-sample SNP calling

cd /scratch/cmb433/fecalRAD/NGS-map/
sh scripts/explore_missingness.sh
Rscript scripts/explore_missingness_further.R
cd /scratch/cmb433/fecalRAD/RAD_faex/

