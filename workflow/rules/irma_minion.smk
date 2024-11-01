"""
Running IRMA on MinION data.
"""


rule irma_raw:
    input:
        "processed_reads/{sample}/{sample}_filtered.fastq.gz",
    output:
        protected(directory("results/{order}/irma-raw/{sample}_{pair}")),
    log:
        ".logs/irma-{order}-raw/{sample}_{pair}.log",
    conda:
        "../envs/irma.yaml"
    threads:
        # take all threads to stop multiple instances of this rule fighting for threads
        workflow.cores
    params:
        config="workflow/config/FLU-{order}-iri-minion.sh"
    shell:
        "IRMA --external-config {params.config} FLU {input} {output} > {log}"
