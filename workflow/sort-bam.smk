from pathlib import Path

def sorted_bams_to_make(wildcards):
    return [
        str(file).replace(".bam", ".sorted.bam")
        for file in Path(".").glob("results/irma/*/*.bam")
        if not str(file).endswith(".sorted.bam")
    ]

rule all:
    input:
        sorted_bams_to_make


rule sort_bam:
    input:
        "{file}.bam"
    output:
        "{file}.sorted.bam"
    shell:
        """
        samtools sort {input} -o {output}
        samtools index {output}
        """
    