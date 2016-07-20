#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Prepare and call make on downsampled BAMs
# ----------------------------------------------------------------------------------------

module load samtools/intel/1.3
export SAMTOOLS=/share/apps/samtools/1.3/intel/bin

for DS_BAM in results/*samp*passed.realn.bam; do

    IND_ID=`echo $DS_BAM | sed -e "s/results\///" -e "s/.PE.*//"`
    IND_ID_W_PE_SE=${IND_ID}.PE

    GENOME_NAME=baboon

    # Make fake versions of precursor files
    touch data/${IND_ID}.R1.fastq
    touch data/${IND_ID}.R2.fastq
    touch reports/${IND_ID}.read1.stats.zip
    touch reports/${IND_ID}.read2.stats.zip
    touch results/${IND_ID}.read1.bwa.${GENOME_NAME}.sai
    touch results/${IND_ID}.PE.bwa.${GENOME_NAME}.sam
    touch results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.sam.bam
    touch results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.sam.bam.sorted.bam.bai
    touch reports/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.aln_stats.txt
    touch results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.fixed.bam
    touch reports/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.aln_stats.pairsfix.txt
    touch results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.fixed.filtered.bam.bai
    touch reports/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.aln_stats.pairsfix.flt.txt
    touch results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.fixed.filtered.postdup.bam.bai 
    touch reports/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.aln_stats.pairsfix.flt.postdup.txt
    touch results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.fixed.filtered.RG.bam
    touch results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.passed.bam.bai
    touch reports/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.aln_stats.passed.txt
    touch results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.passed.bam.list

    # Update last modified time of downsampled BAM
    touch results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.passed.realn.bam
    
    # Index downsampled BAM
    scripts/index_bam.sh results/${IND_ID}.PE.bwa.${GENOME_NAME}.passed.realn.bam
    echo $IND_ID >> data/individual_list.txt

done

# Submit jobs to call make
FIRST=`grep -n "samp" data/individual_list.txt | cut -f1 -d":" | head -n1`
LAST=`grep -n "samp" data/individual_list.txt | cut -f1 -d":" | tail -n1`

qsub -t ${FIRST}-${LAST} pbs/call_make.pbs
