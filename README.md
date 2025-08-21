# Genomics Structural Variant (SV) Pipeline

## Overview
This Snakemake pipeline processes whole genome sequencing (WGS) data to identify structural variants (SVs) in a patient’s genome.

## Input
- Paired-end FASTQ files of WGS data (example: HG002_R1_wgs_chr21.fastq.gz, HG002_R2_wgs_chr21.fastq.gz)
- Reference genome in FASTA format (example: chr21.fa)

## Output
- CSV file containing detected SVs with columns: CHROM, START, END, SIZE, QUAL, FILTER
- Optional columns: TYPE (deletion, duplication, etc.), GENE (overlapping gene), IMPACT (predicted effect)

## Steps in the pipeline
1. Align reads to reference genome using **BWA**
2. Convert and sort using **Samtools**
3. Call structural variants using **Delly, Manta, and Pindel**
4. Merge SV results into final CSV

## Usage
```bash
snakemake --cores 4

## Troubleshooting
- Check input file paths
- Ensure conda environment and tools are installed
- Increase memory if processes fail

