from snakemake.utils import validate, min_version

min_version("7.0.4")


configfile: "config.json"


validate(config, schema="workflow/schemas/config-schema.json")


def build_targets(wildcards):
    """Generate a list of targets to build."""

    if config["platform"] == "minion":
        return (
            expand(
                "combined/{sample}/{sample}_combined.fastq.gz", sample=config["samples"]
            )
            + expand(
                "combined/{sample}/{sample}_filtered.fastq.gz", sample=config["samples"]
            )
            + expand(
                "combined/{sample}/{sample}_{pair}.fastq.gz",
                sample=config["samples"],
                pair=config["pair"],
            )
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

    include: "workflow/rules/filter-trim-qc.smk"

elif config["platform"] == "miseq":

    include: "workflow/rules/trim-qc.smk"

else:
    raise ValueError("'platform' should be 'miseq' or 'minion'")
