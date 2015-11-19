#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Compute heterozygosity
# ----------------------------------------------------------------------------------------

module load vcftools

# Fix issue with VCF not being fully diploid (. used instead of ./.)
gunzip -c ../NGS-map/baboon.INDIV.pass.snp.vcf.gz | \
    sed -e "s:\t\.:\t./.:g" | \
    sed -e "s:\./\.:.:" | \
    gzip -c > ../NGS-map/baboon.INDIV_DIPLOID.pass.snp.vcf.gz

# Computer heterozygosity
vcftools --gzvcf ../NGS-map/baboon.INDIV_DIPLOID.pass.snp.vcf.gz --het \
    --out results/baboon.INDIV_DIPLOID.pass.snp

exit
