from snakemake.utils import validate, min_version

min_version("7.0.4")

configfile: "qc-config.json"
validate(config, schema="schemas/qc-config-schema.json")

rule all:
    input:
        expand("results/qc-raw/{sample}_{n}_fastqc.html", sample=config["samples"], n=[1, 2]),
        expand("results/qc-trimmed/{sample}_{n}_{pair}_fastqc.html", sample=config["samples"], n=[1, 2], pair=config["pair"]),


wildcard_constraints:
    pair="((un)?paired)|(combined)",
    n="1|2"


rule unzip:
    input:
        "raw/{sample}/{sample}_{n}.fastq.gz"
    output:
        "raw/{sample}/{sample}_{n}.fastq"
    shell:
        "gunzip {input}"


rule raw_quality_control:
    input:
        "raw/{sample}/{sample}_{n}.fastq"
    output:
        "results/qc-raw/{sample}_{n}_fastqc.html"
    conda:
        "envs/fastqc.yaml"
    log:
        "logs/fastqc/fastqc_{sample}_{n}.log"
    shell:
        "fastqc --outdir results/qc-raw --format fastq --quiet {input} 2> {log}"


rule trimmed_quality_control:
    input:
        "trimmed/{sample}/{sample}_{n}_{pair}.fastq"
    output:
        "results/qc-trimmed/{sample}_{n}_{pair}_fastqc.html"
    conda:
        "envs/fastqc.yaml"
    log:
        "logs/fastqc/fastqc_{sample}_{n}_{pair}.log"
    shell:
        "fastqc --outdir results/qc-trimmed --format fastq --quiet {input} 2> {log}"


rule trim:
    input:
        "raw/{sample}/{sample}_1.fastq",
        "raw/{sample}/{sample}_2.fastq"
    output:
        "trimmed/{sample}/{sample}_1_paired.fastq",
        "trimmed/{sample}/{sample}_1_unpaired.fastq",
        "trimmed/{sample}/{sample}_2_paired.fastq",
        "trimmed/{sample}/{sample}_2_unpaired.fastq"
    log:
        "logs/trimmomatic_{sample}.log"
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "TrimmomaticPE {input} {output} ILLUMINACLIP:raw/trimlog.fas:2:30:10:2 MINLEN:36 &> {log}"


rule combine_paired_unpaired:
    input:
        "trimmed/{sample}/{sample}_{n}_paired.fastq",
        "trimmed/{sample}/{sample}_{n}_unpaired.fastq"
    output:
        "trimmed/{sample}/{sample}_{n}_combined.fastq"
    shell:
        "cat {input} > {output}"
