from snakemake.utils import validate, min_version

min_version("7.0.4")


configfile: "config.json"


validate(config, schema="schemas/config-schema.json")


def build_targets(wildcards):
    """Generate a list of targets to build."""

    if config["platform"] == "minion":
        return expand(
            "processed_reads/{sample}/{sample}_filtered.fastq.gz",
            sample=config["samples"],
        )

    elif config["platform"] == "miseq":
        return expand(
            "results/qc-raw/{sample}_{n}_fastqc.html",
            sample=config["samples"],
            n=[1, 2],
        ) + expand(
            "results/qc-trimmed/{sample}_{n}_{pair}_fastqc.html",
            sample=config["samples"],
            n=[1, 2],
            pair=config["pair"],
        )

    else:
        raise ValueError("'platform' should be 'miseq' or 'minion'")


rule all:
    input:
        build_targets,


if config["platform"] == "minion":

    include: "rules/preprocess-minion.smk"

elif config["platform"] == "miseq":

    include: "rules/preprocess-miseq.smk"

else:
    raise ValueError("'platform' should be 'miseq' or 'minion'")
