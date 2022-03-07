from snakemake.utils import validate, min_version

min_version("7.0.4")

configfile: "config.json"
validate(config, schema="schemas/config-schema.json")



rule all:
    input:
        expand(
            "results/variants/{sample}_{pair}/{sample}_{pair}.xlsx",
            sample=config["samples"],
            pair=["paired", "combined"]
        )


wildcard_constraints:
    n="1|2",
    pair="(paired)|(combined)"    


rule combine_paired_unpaired:
    input:
        "trimmed/{sample}/{sample}_{n}_paired.fastq",
        "trimmed/{sample}/{sample}_{n}_unpaired.fastq"
    output:
        "trimmed/{sample}/{sample}_{n}_combined.fastq"
    shell:
        "cat {input} > {output}"


checkpoint irma:
    input:
        "trimmed/{sample}/{sample}_1_{pair}.fastq",
        "trimmed/{sample}/{sample}_2_{pair}.fastq",
    output:
        directory("results/irma/{sample}_{pair}")
    log:
        "logs/irma_{sample}_{pair}.log"
    shell:
        "IRMA FLU {input} {output} > {log}"


rule write_gff:
    input:
        "results/irma/{sample}_{pair}/{segment}.fasta"
    output:
        gff="results/irma/{sample}_{pair}/{segment}.gff",
        gffgz="results/irma/{sample}_{pair}/{segment}.gff.gz"
    shell:
        """
        workflow/scripts/write-gff.py --fasta {input} > {output.gff}
        bgzip -c {output.gff} > {output.gffgz}
        tabix -p gff {output.gffgz}
        """


rule summarise_variants:
    input:
        vcf="results/irma/{sample}_{pair}/{segment}.vcf",
        fas="results/irma/{sample}_{pair}/{segment}.fasta",
        gff="results/irma/{sample}_{pair}/{segment}.gff.gz"
    output:
        "results/variants/{sample}_{pair}/{segment}.tsv"
    shell:
        """
        vep --input_file {input.vcf} \
            --output_file {output} \
            --gff {input.gff} \
            --fasta {input.fas} \
            --tab \
            --format vcf  # stops error with empty VCF files 
        """


def aggregate_segments(wildcards):
    """
    Returns a list of segments found by IRMA.
    """
    irma_out_dir = checkpoints.irma.get(**wildcards).output[0]
    segments = [file[:-4] for file in os.listdir(irma_out_dir) if file.endswith(".vcf")]
    return expand("results/variants/{sample}_{pair}/{segment}.tsv", segment=segments, **wildcards)


rule aggregate_summary_tsvs:
    input:
        aggregate_segments
    output:
        "results/variants/{sample}_{pair}/{sample}_{pair}.xlsx"
    shell:
        "workflow/scripts/combine-tables.py {input} --excel {output} --comment '##'"
