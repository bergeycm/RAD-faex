#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Compute heterozygosity
# ----------------------------------------------------------------------------------------

module load vcftools/intel/0.1.13

# Fix issue with VCF not being fully diploid (. used instead of ./.)
gunzip -c NGS-map/baboon_snps_indiv/baboon.INDIV.pass.snp.vcf.gz | \
    sed -e "s:\t\.:\t./.:g" | \
    sed -e "s:\./\.:.:" | \
    gzip -c > NGS-map/baboon_snps_indiv/baboon.INDIV_DIPLOID.pass.snp.vcf.gz

# Compute heterozygosity
vcftools --gzvcf NGS-map/baboon_snps_indiv/baboon.INDIV_DIPLOID.pass.snp.vcf.gz \
    --het --out results/baboon.INDIV_DIPLOID.pass.snp

# And HWE stats on a per SNP basis, including counts of heterozygotes and homozygotes
vcftools --gzvcf NGS-map/baboon_snps_indiv/baboon.INDIV_DIPLOID.pass.snp.vcf.gz \
    --hardy --out results/baboon.pass.snp

# ----------------------------------------------------------------------------------------

# Compute heterozygosity for multi-sample called SNPs
vcftools --gzvcf NGS-map/baboon_snps_multi/baboon.pass.snp.vcf.gz \
    --het --out results/baboon.pass.snp

# And HWE stats on a per SNP basis, including counts of heterozygotes and homozygotes
vcftools --gzvcf NGS-map/baboon_snps_multi/baboon.pass.snp.vcf.gz \
    --hardy --out results/baboon.pass.snp

exit
