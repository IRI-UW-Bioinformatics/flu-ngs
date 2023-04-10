"""
Running IRMA on MinION data.
"""


rule irma_raw:
    input:
        "processed_reads/{sample}/{sample}.fastq.gz",
    output:
        directory("results/{order}/irma-raw/{sample}_{pair}"),
    log:
        ".logs/irma-{order}-raw/{sample}_{pair}.log",
    conda:
        "envs/irma.yaml"
    shell:
        "IRMA FLU-{wildcards.order}-iri-minion {input} {output} > {log}"
