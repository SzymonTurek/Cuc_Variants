# Cuc_Variants
Nextflow pipeline for detecting genomic variants (using the B10 cucumber genome by default). The pipeline combines the following programs:
FastQC, MultiQC, BWA, Bowtie2, bcftools, GATK, freebayes. 

## Instruction:
In the first step, you need to complete the sample_info.csv file so that it contains information about the sample name and and file paths. 

| sample_id  |fastq1  |fastq2 |
| ------------- | ------------- || ------------- |
| Sample1  | input_fastq/Sample1_1.fastq.gz  |nput_fastq/Sample1_2.fastq.gz |
| Sample2  | input_fastq/Sample2_1.fastq.gz  |nput_fastq/Sample2_2.fastq.gz|

## Info
_The pipeline was prepared as part of the National Science Center project UMO-2020/37/B/NZ9/00586._
