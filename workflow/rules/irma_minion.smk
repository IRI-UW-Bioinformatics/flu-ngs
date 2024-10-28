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
    threads:
        # take all threads to stop multiple instances of this rule fighting for threads
        workflow.cores
    shell:
        """
        IRMA \
            --external-config workflow/config/FLU-{wildcards.order}-iri-minion.sh \
            FLU {input} {output} > {log}
        """
