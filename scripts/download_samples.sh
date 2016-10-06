#!/bin/bash

aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq1_BC1_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq1_BC1_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq1_BC2_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq1_BC2_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq1_BC3_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq1_BC3_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq1_BC4_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq1_BC4_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq2_BC1_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq2_BC1_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq3a_BC1_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq3a_BC1_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq3b_BC1_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/FecalSeq/FecalSeq3b_BC1_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/KafueRAD/Kafue1_BC1_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/KafueRAD/Kafue1_BC1_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/KafueRAD/Kafue2a_BC1_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/KafueRAD/Kafue2a_BC1_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/KafueRAD/Kafue2a_BC2_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/KafueRAD/Kafue2a_BC2_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/KafueRAD/Kafue2a_BC3_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/KafueRAD/Kafue2a_BC3_R2.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/KafueRAD/Kafue2a_BC4_R1.fastq.gz demultiplex/data/
aws s3 cp s3://chiou-genomics-bucket/KafueRAD/Kafue2a_BC4_R2.fastq.gz demultiplex/data/

exit
