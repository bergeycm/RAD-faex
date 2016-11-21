#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download repository and submodules
# ----------------------------------------------------------------------------------------

git clone --recursive https://github.com/bergeycm/RAD-faex.git
cd RAD-faex

# ----------------------------------------------------------------------------------------
# --- Download samples
# ----------------------------------------------------------------------------------------

sh scripts/download_samples.sh

# ----------------------------------------------------------------------------------------
# --- Demultiplex samples
# ----------------------------------------------------------------------------------------

LIBRARIES=(FecalSeq1 FecalSeq2 FecalSeq3a FecalSeq3b Kafue1 Kafue2a)

cd demultiplex/

for p in ${LIBRARIES[*]}; do
    if [ $p == ${LIBRARIES[0]} ] ; then
        job=$(qsub -v PREFIX=$p pbs/DeMultiplex.pbs)
    else
        job=$(qsub -v PREFIX=$p -Wdepend=afterok:${job} pbs/DeMultiplex.pbs)
    fi
done

cd ..

# Move samples to mapping data folder (concatenating repeat runs as needed)
sh scripts/merge_and_rename_samples.sh

# ----------------------------------------------------------------------------------------
# --- Download genome
# ----------------------------------------------------------------------------------------

cd NGS-map

sh genomes/download_papAnu2.sh

cd ..

# ----------------------------------------------------------------------------------------
# --- Change from aln to mem
# ----------------------------------------------------------------------------------------

cd NGS-map/

mv scripts/align_aln.sh scripts/align.sh

cd ..

# ----------------------------------------------------------------------------------------
# --- Map and genotype samples
# ----------------------------------------------------------------------------------------

cd NGS-map

# Edit configuration variables
sed -e "s:READ1=.*:READ1=./data/\${IND_ID}.read1.fastq:g" -i config.mk
sed -e "s:READ2=.*:READ2=./data/\${IND_ID}.read2.fastq:g" -i config.mk
sed -e "s:GENOME_FA=.*:GENOME_FA=genomes/papAnu2/papAnu2.fa:g" -i config.mk
sed -e "s:GENOME_NAME=.*:GENOME_NAME=baboon:g" -i config.mk
sed -e "s:MAPQUAL=.*:MAPQUAL=0:g" -i config.mk
sed -e "s:MARK_DUPS=.*:MARK_DUPS=FALSE:g" -i config.mk

# Create list of individuals
ls data/*.fastq | xargs -n 1 basename | sed -e "s/\.read[0-9]\.fastq//g" | \
    uniq > data/individual_list.txt

# Get number of samples
N=$(wc -l data/individual_list.txt | cut -d ' ' -f1)

# Run the mapping script
qsub -t 1-${N} pbs/call_make.pbs

# ----------------------------------------------------------------------------------------
# --- Downsample to equalize coverage in blood-feces pairs
# ----------------------------------------------------------------------------------------

module load samtools/intel/1.3

# Do downsampling
perl ../scripts/downsample_bloods.pl

# Fake the precursor files to get ready, and then call Make on these downsampled samples
perl ../scripts/prepare_to_process_downsampled.sh

# Fix the headers in the BAM files
module load picard-tools/1.129
module load samtools/intel/1.3

PICARD_TOOLS_ROOT=/share/apps/picard-tools/1.129/
SAMTOOLS_ROOT=/share/apps/samtools/1.3/intel/bin/

for BAM in results/*samp*.PE.bwa.baboon.passed.realn.bam; do
    cp ${BAM} ${BAM}.backup
    DS_ID=`echo $BAM | sed -e "s:results/::" -e "s/\.bwa.*//"`
    java -jar ${PICARD_TOOLS_ROOT}/picard.jar AddOrReplaceReadGroups \
        INPUT=$BAM.backup \
        OUTPUT=$BAM \
        RGLB=${DS_ID} \
        RGPL=Illumina \
        RGPU=Group1 \
        RGSM=${DS_ID}
    ${SAMTOOLS_ROOT}/samtools index $BAM
done

# ----------------------------------------------------------------------------------------
# --- Do multi-sample SNP calling and filtration
# ----------------------------------------------------------------------------------------

# Perform multi-sample SNP-calling, one job per chromosome (this includes X as "21")
qsub -t 1-21 pbs/call_gatk_genotyper.pbs

# ----------------------------------------------------------------------------------------

# Filter SNPs (This used to be for only autosomes. Changed to 1-21 to do X too.)
qsub -t 1-21 pbs/filter_gatk_snps.pbs

# ----------------------------------------------------------------------------------------

# Call pipeline in comparison mode to merge multi-sample SNPs, convert VCF file to PED,
# and make binary PED (BED)

# make -s -f full_analysis.mk compare

qsub -t 1-21 pbs/gatk_DoC.pbs

# ----------------------------------------------------------------------------------------

# Merge SNP calls from multiple chromsomes into autosomal, X, and full datasets

scripts/merge_snps.sh baboon

# ----------------------------------------------------------------------------------------

# Rename folder of SNPs to indicate that these were called in multi-sample SNP mode
mv baboon_snps{,_multi}

mv baboon.* baboon_snps_multi/

# ----------------------------------------------------------------------------------------
# --- Now call GATK to generate SNP sets that are NOT called in multi-sample mode
# ----------------------------------------------------------------------------------------

# Copy in GATK-individual-mode PBS script from the other repo
cp ../pbs/call_gatk_genotyper_indiv.pbs pbs/

qsub -t 1-21 pbs/call_gatk_genotyper_indiv.pbs

# Clean up
rm baboon_snps/*tmp*

# Replace QUAL score of "inf" with one more than the maximum QUAL score found
###sh ../RAD-faex/scripts/replace_inf_in_indiv_vcfs.sh

# Clean up some more
# rm baboon_snps/chr*.INDIV.raw.snps.indels.vcf_BACKUP

# And filter
cp ../pbs/filter_gatk_snps_indiv.pbs pbs/
qsub -t 1-21 pbs/filter_gatk_snps_indiv.pbs

# Fix headers
# Get rid of, e.g., chr10.raw.snps.indels.tmp12_
# And replace it with info on sample used to determine downsampling level
BAMS=(`ls results/*.PE.bwa.baboon.passed.realn.bam`)

# Downsampled samples to fix
DS_SAMPS=($(grep "^#CHROM" baboon_snps/chr20.INDIV.pass.snp.vcf | \
    tr "\t" "\n" | grep "chr" | sed -e "s/chr20/chr[0-9]\*/"))

for ((i=0; i < ${#DS_SAMPS}; i++)); do

    if [[ ${DS_SAMPS[$i]} = *[!\ ]* ]]; then
    	REPLACEE=${DS_SAMPS[$i]}
        echo "    Replacing ${DS_SAMPS[$i]}..."
        DS_IDX=`echo ${DS_SAMPS[$i]} | sed -e "s/.*tmp\([0-9]*\).*/\1/"`
        echo "    ...with item indexed ${DS_IDX}...";
        REPLACER=`echo ${BAMS[$DS_IDX - 1]} | sed -e "s/.*results\///" -e "s/\.PE.*//"`
        echo "    ...with ${REPLACER}.";

        for file in baboon_snps/chr*.INDIV.pass.snp.vcf; do
            sed -e "s/$REPLACEE/$REPLACER/g" -i $file
        done
    fi
done

# Steal steps from make to merge multi-sample SNPs, convert VCF file to PED,
# and make binary PED (BED)

# Merge non-multi-sample SNPs (autosomes only)
module load vcftools/intel/0.1.13
module load plink/intel/1.90p

vcf-concat baboon_snps/chr[0-9]*.INDIV.pass.snp.vcf | \
    gzip -c > baboon_snps/baboon.INDIV.pass.snp.vcf.gz

# Convert VCF file to PED
vcftools --gzvcf baboon_snps/baboon.INDIV.pass.snp.vcf.gz --plink \
    --out baboon_snps/baboon.INDIV.pass.snp

# Edit the MAP file (baboon.pass.snp.map) and get rid of the "chr"
# VCF uses, e.g., "chr10" whereas plink wants just "10"
sed -i -e 's/^chr//' baboon_snps/baboon.INDIV.pass.snp.map

# Make binary PED file
plink --noweb --file baboon_snps/baboon.INDIV.pass.snp --make-bed \
    --out baboon_snps/baboon.INDIV.pass.snp

# ----------------------------------------------------------------------------------------

# Rename folder of SNPs to indicate that these were called in non-multi-sample SNP mode
mv baboon_snps{,_indiv}

# ========================================================================================
# --- Commands originally in RADfaex_cmds.sh
# ========================================================================================

cd ..

module load r/intel/3.2.2

# ----------------------------------------------------------------------------------------
# --- Simulate ddRAD digestion -----------------------------------------------------------
# ----------------------------------------------------------------------------------------

# The folder data/papAnu2/ should contain fastas for each chromosome named chr[1-20,X].fa
# Generates results/chr*.bed

qsub -t 1-21 pbs/digest_sim.pbs

# ----------------------------------------------------------------------------------------
# --- Combine BED files of RADtags
# ----------------------------------------------------------------------------------------

# Generates results/RADtags.bed

scripts/combine_RADtag_beds.sh

# ----------------------------------------------------------------------------------------
# --- Compute GC within all RADtags (from simulation) as well as in the larger region
# ----------------------------------------------------------------------------------------

# Called with:
#     pbs/compute_GC.pbs
# Generates results/RADtags.gc.bed,
#           results/RADtags.slop${REGION_EXPANSION}.bed
#           results/RADtags.slop${REGION_EXPANSION}.gc.bed

qsub pbs/compute_GC.pbs

# ----------------------------------------------------------------------------------------
# --- Download CpG islands
# ----------------------------------------------------------------------------------------

# Generates data/cpgIslandExtUnmasked.chr.txt
#           data/cpgIslandExtUnmasked.chr.bed

scripts/download_CpG_islands.sh

# ----------------------------------------------------------------------------------------
# --- Find CpGs
# ----------------------------------------------------------------------------------------

# Generates results/papAnu2.CpG.bed

qsub pbs/call_CpG_finder.pbs

# ----------------------------------------------------------------------------------------
# --- Find closest CpG
# ----------------------------------------------------------------------------------------

# Generates results/RADtags.nearestCpGisland.bed

qsub pbs/call_closest_CpG_finder.pbs

# ----------------------------------------------------------------------------------------
# --- Compute coverage within all RADtags (from simulation)
# ----------------------------------------------------------------------------------------

# Generates results/${ID}.gt1.cov.bed

TOTAL_IND=$(ls NGS-map/results/*.PE.bwa.baboon.passed.realn.bam | grep -v "samp" | wc -l | cut -d' ' -f1)

qsub -t 0-$((TOTAL_IND - 1)) pbs/get_RAD_cov.pbs

# ----------------------------------------------------------------------------------------
# --- Combine RADtag coverage from all individuals
# ----------------------------------------------------------------------------------------

# Generates results/RADtag.coverage.all.txt

sh scripts/combine_coverage.sh

# ----------------------------------------------------------------------------------------
# --- Plot heatmap of coverage and RADtag length
# ----------------------------------------------------------------------------------------

# Generates reports/*.gt1.cov.pdf

Rscript scripts/coverage_length_heatmap.R

# ----------------------------------------------------------------------------------------
# --- Combine RADtag info
# ----------------------------------------------------------------------------------------

# Generates results/RADtags.info.bed

qsub pbs/combine_RAD_info.pbs

# ----------------------------------------------------------------------------------------
# --- Explore how coverage varies by sample, sample type, mapped read count, etc.
# ----------------------------------------------------------------------------------------

# Generates reports/RADtag_count_gtNreads.pdf
#           reports/RADtag_count_curves_by_total_mapped.pdf
# Note the work in progress: Plot residuals by various variables

Rscript scripts/plot_coverage.R

# ----------------------------------------------------------------------------------------
# --- Model RADtag coverage (to explain inter-RADtag variation) : SKIP
# ----------------------------------------------------------------------------------------

# See inter.RADtag.regression.R

# qsub pbs/model_RADtag_cvg_1000.pbs
# qsub pbs/model_RADtag_cvg_5000.pbs
# qsub pbs/model_RADtag_cvg_10000.pbs
# qsub pbs/model_RADtag_cvg_20000.pbs

# ----------------------------------------------------------------------------------------
# --- Compute and explore heterozygosity
# ----------------------------------------------------------------------------------------

sh scripts/compute_het.sh
Rscript scripts/parse_heterozygosity.R results/baboon.pass.snp.het

# ----------------------------------------------------------------------------------------
# --- Quantify allelic dropout by comparing blood and fecal DNA from single individual
# ----------------------------------------------------------------------------------------

module load vcftools

# Multi-sample mode
perl scripts/quantify_ADO.pl NGS-map/baboon_snps_multi/baboon.pass.snp.vcf.gz

# Individual mode
perl scripts/quantify_ADO.pl \
    NGS-map/baboon_snps_indiv/baboon.INDIV_DIPLOID.pass.snp.vcf.gz

sh scripts/parse_all_discordance_matrices.sh

Rscript scripts/explore_discordance.R results/discordance.multi.txt
Rscript scripts/explore_discordance.R results/discordance.indiv.txt

# ----------------------------------------------------------------------------------------
# --- Compute stats on missingness
# ----------------------------------------------------------------------------------------

# This uses SNPs called within one individual, not the results of multi-sample SNP calling

cd /scratch/cmb433/fecalRAD/NGS-map/
sh scripts/explore_missingness.sh
Rscript scripts/explore_missingness_further.R
cd /scratch/cmb433/fecalRAD/RAD_faex/

# ----------------------------------------------------------------------------------------
# --- Make MDS plot
# ----------------------------------------------------------------------------------------


exit
