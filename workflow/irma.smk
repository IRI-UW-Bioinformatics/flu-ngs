from snakemake.utils import validate, min_version

min_version("7.0.4")

configfile: "config.json"
validate(config, schema="schemas/config-schema.json")


rule all:
    input:
        expand("results/irma/{sample}_{group}", sample=config["samples"], group=["paired", "combined"])


wildcard_constraints:
    n="1|2",
    

rule combine_paired_unpaired:
    input:
        "trimmed/{sample}/{sample}_{n}_paired.fastq",
        "trimmed/{sample}/{sample}_{n}_unpaired.fastq"
    output:
        "trimmed/{sample}/{sample}_{n}_combined.fastq"
    shell:
        "cat {input} > {output}"


rule irma:
    wildcard_constraints:
        group="(paired)|(combined)"
    input:
        "trimmed/{sample}/{sample}_1_{group}.fastq",
        "trimmed/{sample}/{sample}_2_{group}.fastq",
    output:
        directory("results/irma/{sample}_{group}")
    log:
        "logs/irma_{sample}_{group}.log"
    shell:
        "IRMA FLU {input} {output} > {log}"
