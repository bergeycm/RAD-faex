#!/bin/sh

# ========================================================================================
# --- Find closest feature (CpG site, CpG island) and local density (CpG)
# ========================================================================================

module load bedtools/intel/2.25.0

# === CpG sites ==========================================================================

# --- Find closest CpG site --------------------------------------------------------------

echo "Finding closest CpG site..."

bedtools closest -d \
	-a results/RADtags.bed \
	-b results/papAnu2.CpG.bed \
	> results/RADtags.nearestCpG.bed

# --- Find CpG count in RAD tag ----------------------------------------------------------

echo "Finding CpG count in RAD tag..."

bedtools intersect -c \
	-a results/RADtags.bed \
	-b results/papAnu2.CpG.bed \
	> results/RADtags.CpGcount.bed

# --- Find local CpG density (actually count) --------------------------------------------

echo "Finding local CpG density..."

REGION_EXPANSION=5000

# Make genome file
cut -f1-2 NGS-map/genomes/papAnu2/papAnu2.fa.fai \
    > NGS-map/genomes/papAnu2/papAnu2.genome
    
bedtools slop \
    -i results/RADtags.bed \
    -b ${REGION_EXPANSION} \
    -g NGS-map/genomes/papAnu2/papAnu2.genome \
    > results/RADtags.slop${REGION_EXPANSION}.bed

bedtools intersect -c \
	-a results/RADtags.slop${REGION_EXPANSION}.bed \
	-b results/papAnu2.CpG.bed \
	> results/RADtags.slop${REGION_EXPANSION}.CpGcount.bed

# === CpG islands ========================================================================

# --- Find closest CpG island ------------------------------------------------------------

echo "Finding closest CpG island..."

bedtools closest -d \
	-a results/RADtags.bed \
	-b data/cpgIslandExtUnmasked.chr.bed \
	> results/RADtags.nearestCpGisland.bed

exit
