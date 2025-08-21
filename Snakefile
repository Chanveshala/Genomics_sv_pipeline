# Snakefile for structural variant detection
# Uses system-installed Delly2, Manta, Pindel, BWA, and SAMtools
# Input: FASTQ files and reference genome
# Output: CSV of structural variants

import os
import pandas as pd

SAMPLES = ["HG002"]  # list your sample names here

REF = "chr21.fa"
FASTQ_DIR = "data"
ALN_DIR = "aligned"
SV_DIR = "sv"

rule all:
    input:
        expand(f"{SV_DIR}/{{sample}}_structural_variants.csv", sample=SAMPLES)

# Step 1: Align FASTQ to reference genome
rule align_bwa:
    input:
        r1=f"{FASTQ_DIR}/{{sample}}_R1_wgs_chr21.fastq.gz",
        r2=f"{FASTQ_DIR}/{{sample}}_R2_wgs_chr21.fastq.gz",
        ref=REF
    output:
        bam=f"{ALN_DIR}/{{sample}}.bam"
    shell:
        """
        mkdir -p {ALN_DIR}
        bwa mem {input.ref} {input.r1} {input.r2} | \
        samtools view -bS - | \
        samtools sort -o {output.bam}
        samtools index {output.bam}
        """

# Step 2a: Call SVs with Delly2
rule call_delly:
    input:
        bam=f"{ALN_DIR}/{{sample}}.bam"
    output:
        vcf=f"{SV_DIR}/delly/{{sample}}.vcf"
    shell:
        "mkdir -p {SV_DIR}/delly && "
        "delly call -g {REF} -o {output.vcf} {input.bam}"

# Step 2b: Call SVs with Manta
rule call_manta:
    input:
        bam=f"{ALN_DIR}/{{sample}}.bam"
    output:
        vcf=f"{SV_DIR}/manta/{{sample}}.vcf"
    shell:
        "mkdir -p {SV_DIR}/manta && "
        "configManta.py --bam {input.bam} --referenceFasta {REF} --runDir {SV_DIR}/manta/{{wildcards.sample}}_manta && "
        "{SV_DIR}/manta/{{wildcards.sample}}_manta/runWorkflow.py -m local -j 4 && "
        "cp {SV_DIR}/manta/{{wildcards.sample}}_manta/results/variants/diploidSV.vcf {output.vcf}"

# Step 2c: Call SVs with Pindel
rule call_pindel:
    input:
        bam=f"{ALN_DIR}/{{sample}}.bam"
    output:
        vcf=f"{SV_DIR}/pindel/{{sample}}.vcf"
    shell:
        """
        mkdir -p {SV_DIR}/pindel
        pindel -f {REF} -i {input.bam} -o {SV_DIR}/pindel/{{wildcards.sample}}
        pindel2vcf -p {SV_DIR}/pindel/{{wildcards.sample}} -r {REF} -R chr21 -d 2025-08-21 -v {output.vcf}
        """

# Step 3: Merge SVs into single CSV
rule merge_sv:
    input:
        delly=f"{SV_DIR}/delly/{{sample}}.vcf",
        manta=f"{SV_DIR}/manta/{{sample}}.vcf",
        pindel=f"{SV_DIR}/pindel/{{sample}}.vcf"
    output:
        csv=f"{SV_DIR}/{{sample}}_structural_variants.csv"
    run:
        def parse_vcf(vcf_file, source):
            rows = []
            with open(vcf_file) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    chrom, pos, _, ref, alt, qual, filt, info = parts[:8]
                    svtype = "NA"
                    svsize = "NA"
                    # parse info field
                    for field in info.split(";"):
                        if field.startswith("SVTYPE="):
                            svtype = field.split("=")[1]
                        if field.startswith("SVLEN="):
                            svsize = field.split("=")[1]
                    rows.append({
                        "CHROM": chrom,
                        "START": pos,
                        "END": pos,
                        "SIZE": svsize,
                        "QUAL": qual,
                        "FILTER": filt,
                        "SOURCE": source,
                        "SVTYPE": svtype
                    })
            return rows

        all_rows = []
        all_rows += parse_vcf(input.delly, "Delly")
        all_rows += parse_vcf(input.manta, "Manta")
        all_rows += parse_vcf(input.pindel, "Pindel")

        df = pd.DataFrame(all_rows)
        os.makedirs(os.path.dirname(output.csv), exist_ok=True)
        df.to_csv(output.csv, index=False)
