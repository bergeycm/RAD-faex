#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Compute coverage within all RADtags (from simulation)
# ----------------------------------------------------------------------------------------

module load bedtools

function compute_cov {
    BAM=$1
    ID=$(echo $BAM | sed "s:.*\/\(.*\)\.PE\.bwa.*:\\1:")
    bedtools coverage \
        -counts \
        -a results/RADtags.bed \
        -b $BAM > results/${ID}.cov.bed
    
    # Also make file with those with no coverage removed
    awk '{ if ($5 != 0) print $0}' results/${ID}.cov.bed > results/${ID}.gt1.cov.bed 
}

# Run all if no BAM file passed
if [ $# -eq 0 ]; then
    for THIS_BAM in NGS-map/results/*.PE.bwa.baboon.passed.realn.bam; do
        compute_cov $THIS_BAM
    done;
else
    # Otherwise process the BAM file passed.
    compute_cov $1
fi

exit
