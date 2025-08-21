import os

configfile: "config.yaml"

SAMPLES = config["samples"]
REF = config["reference"]

rule all:
    input:
        expand("results/{sample}_sv_final.csv", sample=SAMPLES)

rule align:
    input:
        fastq1 = "data/{sample}_R1.fastq",
        fastq2 = "data/{sample}_R2.fastq",
        ref = REF
    output:
        bam = "results/{sample}.bam"
    shell:
        """
        bwa mem {input.ref} {input.fastq1} {input.fastq2} | samtools view -bS - > {output.bam}
        samtools sort -o {output.bam} {output.bam}
        samtools index {output.bam}
        """

rule manta:
    input:
        bam = "results/{sample}.bam",
        ref = REF
    output:
        vcf = "results/{sample}_manta.vcf"
    shell:
        """
        configManta.py --bam {input.bam} --referenceFasta {input.ref} --runDir results/manta_{wildcards.sample}
        results/manta_{wildcards.sample}/runWorkflow.py -m local -j 4
        cp results/manta_{wildcards.sample}/results/variants/diploidSV.vcf {output.vcf}
        """

rule delly:
    input:
        bam = "results/{sample}.bam",
        ref = REF
    output:
        vcf = "results/{sample}_delly.vcf"
    shell:
        """
        delly call -g {input.ref} -o {output.vcf} {input.bam}
        """

rule merge_vcfs_to_csv:
    input:
        manta = "results/{sample}_manta.vcf",
        delly = "results/{sample}_delly.vcf"
    output:
        csv = "results/{sample}_sv_final.csv"
    shell:
        """
        vcf-concat {input.manta} {input.delly} | vcf-sort | vcf2csv -o {output.csv}
        """
