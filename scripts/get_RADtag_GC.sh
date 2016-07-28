#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Compute GC within all RADtags (from simulation) as well as in the larger region
# ----------------------------------------------------------------------------------------

module load bedtools/intel/2.25.0

# --- Compute GC within RAD tags
    
bedtools nuc \
    -fi NGS-map/genomes/papAnu2/papAnu2.fa \
    -bed results/RADtags.bed \
    > results/RADtags.gc.bed

# --- Compute GC around RAD tags

REGION_EXPANSION=5000

# Make genome file
cut -f1-2 NGS-map/genomes/papAnu2/papAnu2.fa.fai \
    > NGS-map/genomes/papAnu2/papAnu2.genome
    
bedtools slop \
    -i results/RADtags.bed \
    -b ${REGION_EXPANSION} \
    -g NGS-map/genomes/papAnu2/papAnu2.genome \
    > results/RADtags.slop${REGION_EXPANSION}.bed

bedtools nuc \
    -fi NGS-map/genomes/papAnu2/papAnu2.fa \
    -bed results/RADtags.slop${REGION_EXPANSION}.bed \
    > results/RADtags.slop${REGION_EXPANSION}.gc.bed

exit
