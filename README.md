# Genomics Structural Variant (SV) Pipeline

## Overview
This Snakemake workflow processes patient whole genome sequencing (WGS) data to identify **structural variants (SVs)**.  
The workflow is designed to be modular, reproducible, and easy to troubleshoot. It uses multiple SV callers to maximize accuracy.

**Pipeline steps:**
1. Alignment of reads to a reference genome
2. Sorting and indexing BAM files
3. Structural variant calling using Delly2, Manta, and Pindel
4. Merging SV calls into a final CSV file for downstream interpretation

---

## Pipeline Diagram

FASTQ (R1 & R2) ──► BWA MEM ──► BAM ──► samtools sort & index
│
▼
┌───────────────┐
│ SV Calling │
│ Delly / Manta │
│ Pindel │
└───────┬───────┘
▼
Merge SVs
▼
final_SV.csv

---

## Requirements

- **Operating System:** macOS / Linux
- **Software & Tools:**
  - [Conda](https://docs.conda.io/en/latest/miniconda.html)
  - bwa
  - samtools
  - delly
  - manta
  - pindel
  - pandas (Python library for CSV processing)

- **Conda environment:** `envs/sv_env_min.yaml` included

---

## Setup

1. **Clone the repository**
```bash
git clone https://github.com/Chanveshala/Genomics_sv_pipeline.git
cd Genomics_sv_pipeline

2. Create the Conda environment
conda env create -f envs/sv_env_min.yaml
conda activate sv_env_min

3. Check tools installation
bwa
samtools
delly
manta
pindel
All should respond with their version information.


Usage
Run Snakemake to process data:
snakemake --cores 4


Output directories:
aligned/ → BAM files
results/final_SV.csv → CSV file with detected structural variants
CSV Columns:
CHROM: Chromosome
START: Start position
END: End position
SIZE: SV size
QUAL: Quality score
FILTER: Filter status
(Optional) Additional columns can be added, e.g., gene affected, known pathogenicity, population frequency

Troubleshooting Tips
Missing input files: Check paths to FASTQ or reference genome
Environment issues: Ensure Conda environment is active and tools are installed
Memory issues: Increase resources in Snakemake rule definitions
Failed rules: Check .snakemake/log/ for detailed error messages



Version Control
All changes are tracked using Git. Commit history shows:
Initial pipeline setup
Added SV calling rules
Updated README with instructions


References
BWA
SAMtools
Delly2
Manta
Pindel
