from snakemake.utils import validate, min_version

min_version("7.0.4")

configfile: "config.json"
validate(config, schema="schemas/config-schema.json")

rule all:
    input:
        expand("results/qc/{sample}_{n}_fastqc.html", sample=config["samples"], n=[1, 2])


rule gunzip:
    input:
        "{file}.fastq.gz"
    output:
        "{file}.fastq"
    shell:
        "gunzip {input}"


rule fastq_quality_control:
    input:
        "reads/{sample}/{sample}_{n}.fastq"
    output:
        "results/qc/{sample}_{n}_fastqc.html"
    shell:
        "fastqc --outdir results/qc --format fastq --quiet {input}"