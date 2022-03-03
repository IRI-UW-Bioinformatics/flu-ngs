from snakemake.utils import validate, min_version

min_version("7.0.4")

configfile: "config.json"
validate(config, schema="schemas/config-schema.json")

rule all:
    input:
        expand("results/qc-raw/{sample}_{n}_fastqc.html", sample=config["samples"], n=[1, 2]),
        expand("results/qc-trimmed/{sample}_{n}_{pair}_fastqc.html", sample=config["samples"], n=[1, 2], pair=["paired", "unpaired"]),


wildcard_constraints:
    fastq=".fastq$",
    pair="(un)?paired",
    n="1|2"


rule unzip:
    input:
        "{fastq}.gz"
    output:
        "{fastq}"
    shell:
        "gunzip {input}"


rule raw_quality_control:
    input:
        "raw/{sample}/{sample}_{n}.fastq"
    output:
        "results/qc-raw/{sample}_{n}_fastqc.html"
    shell:
        "fastqc --outdir results/qc-raw --format fastq --quiet {input}"


rule trimmed_quality_control:
    input:
        "trimmed/{sample}/{sample}_{n}_{pair}.fastq"
    output:
        "results/qc-trimmed/{sample}_{n}_{pair}_fastqc.html"
    shell:
        "fastqc --outdir results/qc-trimmed --format fastq --quiet {input}"


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
        "trimmed/{sample}/trimmomatic.log"
    shell:
        "TrimmomaticPE {input} {output} ILLUMINACLIP:raw/trimlog.fas:2:30:10:2 MINLEN:36 &> {log}"