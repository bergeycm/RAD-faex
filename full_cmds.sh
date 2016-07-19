#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download samples
# ----------------------------------------------------------------------------------------

scripts/download_samples.sh

# ----------------------------------------------------------------------------------------
# --- Demultiplex samples
# ----------------------------------------------------------------------------------------

LIBRARIES=(FecalSeq1 FecalSeq2 FecalSeq3a FecalSeq3b Kafue1 Kafue2a)

cd demultiplex/

for p in ${LIBRARIES[*]}
do
if [ $p == ${LIBRARIES[0]} ] ; then
job=$(qsub -v PREFIX=$p pbs/DeMultiplex.pbs)
else
job=$(qsub -v PREFIX=$p -Wdepend=afterok:${job} pbs/DeMultiplex.pbs)
fi
done

cd ..

# Move samples to mapping data folder (concatenating repeat runs as needed)
scripts/merge_and_rename_samples.sh

# ----------------------------------------------------------------------------------------
# --- Map and genotype samples
# ----------------------------------------------------------------------------------------

cd NGS-map/

# Edit configuration variables
sed -e "s:READ1=.*:READ1=./data/\${IND_ID}.read1.fastq:g" -i config.mk
sed -e "s:READ2=.*:READ2=./data/\${IND_ID}.read2.fastq:g" -i config.mk
sed -e "s:GENOME_FA=.*:GENOME_FA=genomes/papAnu2/papAnu2.fa:g" -i config.mk
sed -e "s:GENOME_NAME=.*:GENOME_NAME=baboon:g" -i config.mk
sed -e "s:MAPQUAL=.*:MAPQUAL=0:g" -i config.mk
sed -e "s:MARK_DUPS=.*:MARK_DUPS=FALSE:g" -i config.mk

# Create list of individuals
ls data/*.fastq | xargs -n 1 basename | sed -e "s/\.read[0-9]\.fastq//g" | uniq > data/individual_list.txt

# Get number of samples
N=$(wc -l data/individual_list.txt | cut -d ' ' -f1)

# Run the mapping script
qsub -t 1-${N} pbs/call_make.pbs

# Run the genotyping and filtering scripts
qsub -t 1-21 pbs/call_gatk_genotyper.pbs
qsub -t 1-21 pbs/filter_gatk_snps.pbs

# ----------------------------------------------------------------------------------------
# --- Profile samples with Kraken - Now skipped
# ----------------------------------------------------------------------------------------
# 
# cd data
# mkdir profiled_samples
# 
# module load metaphlan/2.0.0 
# module load bowtie2/intel/2.2.3
# 
# fec_samples=`ls *.fastq`
# 
# for f in ${fec_samples}; do
#     python ${METAPHLAN_ROOT}/metaphlan2.py \
#         --mpa_pkl ${MPADB_DIR}/mpa_v20_m200.pkl \
#         --bowtie2db ${MPADB_DIR}/mpa_v20_m200 \
#         --bt2_ps very-sensitive \
#         --input_type multifastq \
#         --bowtie2out ${f}.bt2out > profiled_samples/${f}.metaphlan.txt
# done
# 
# export PATH=$PATH:${HOME}/downloads/jellyfish-1.1.11/bin/bin/
# 
# module load kraken/intel/0.10.4
# 
# DBDIR=${SCRATCH}/kraken_db
# mkdir $DBDIR
# DBNAME=${DBDIR}/kraken_db
# 
# # Build database
# kraken-build --standard --threads 12 --db $DBNAME
# 
# kraken --db $DBNAME seqs.fa
#
# cd ..

# ----------------------------------------------------------------------------------------
# --- Downsample to equalize coverage in blood-feces pairs
# ----------------------------------------------------------------------------------------

# Do downsampling
perl ../RAD-faex/scripts/downsample_bloods.pl

# Fake the precursor files to get ready, and then call Make on these downsampled samples
perl ../RAD-faex/scripts/prepare_to_process_downsampled.sh

# Fix the headers in the BAM files
module load picard-tools/1.129
module load samtools

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
    samtools index $BAM
done










