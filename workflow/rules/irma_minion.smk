"""
Running IRMA on MinION data.
"""


rule irma_raw:
    input:
        "processed_reads/{sample}/{sample}_filtered.fastq.gz",
    output:
        directory("results/{order}/irma-raw/{sample}_{pair}"),
    log:
        ".logs/irma-{order}-raw/{sample}_{pair}.log",
    conda:
        "../envs/irma.yaml"
    threads: workflow.cores  # take all threads to stop multiple instances of this rule fighting for threads
    shell:
        "IRMA FLU-{wildcards.order}-iri-minion {input} {output} > {log}"
