from snakemake.utils import validate, min_version

min_version("7.0.4")


configfile: "config.json"


validate(config, schema="schemas/config-schema.json")


rule all:
    input:
        expand(
            "results/qc-raw/{sample}_{n}_fastqc.html",
            sample=config["samples"],
            n=[1, 2],
        ),
        expand(
            "results/qc-trimmed/{sample}_{n}_{pair}_fastqc.html",
            sample=config["samples"],
            n=[1, 2],
            pair=config["pair"],
        ),


wildcard_constraints:
    pair="((un)?paired)|(combined)",
    n="1|2",


rule unzip:
    input:
        "raw/{sample}/{sample}_{n}.fastq.gz",
    output:
        "raw/{sample}/{sample}_{n}.fastq",
    shell:
        "gunzip {input}"


rule raw_quality_control:
    input:
        "raw/{sample}/{sample}_{n}.fastq",
    output:
        "results/qc-raw/{sample}_{n}_fastqc.html",
    conda:
        "envs/fastqc.yaml"
    log:
        ".logs/fastqc/fastqc_{sample}_{n}.log",
    shell:
        "fastqc --outdir results/qc-raw --format fastq --quiet {input} 2> {log}"


rule trimmed_quality_control:
    input:
        "processed_reads/{sample}/{sample}_{n}_{pair}.fastq",
    output:
        "results/qc-trimmed/{sample}_{n}_{pair}_fastqc.html",
    conda:
        "envs/fastqc.yaml"
    log:
        ".logs/fastqc/fastqc_{sample}_{n}_{pair}.log",
    shell:
        "fastqc --outdir results/qc-trimmed --format fastq --quiet {input} 2> {log}"


rule trim:
    input:
        expand("raw/{{sample}}/{{sample}}_{n}.fastq", n=(1, 2)),
    output:
        expand(
            "processed_reads/{{sample}}/{{sample}}_{n}_{pair}.fastq",
            n=(1, 2),
            pair=("paired", "unpaired"),
        ),
    log:
        ".logs/trimmomatic/trimmomatic_{sample}.log",
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "TrimmomaticPE {input} {output} ILLUMINACLIP:raw/trimlog.fas:2:30:10:2 MINLEN:36 &> {log}"


rule combine_paired_unpaired:
    input:
        expand("processed_reads/{{sample}}/{{sample}}_{n}_paired.fastq", n=(1, 2)),
    output:
        "processed_reads/{sample}/{sample}_{n}_combined.fastq",
    shell:
        "cat {input} > {output}"
