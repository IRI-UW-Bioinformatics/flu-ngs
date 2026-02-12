"""
Running IRMA on MiSeq data.
"""


rule irma_raw:
    input:
        expand("processed_reads/{{sample}}/{{sample}}_{n}_{{pair}}.fastq", n=(1, 2)),
    output:
        directory("results/{order}/irma-raw/{sample}_{pair}")
    log:
        ".logs/irma-{order}-raw/{sample}_{pair}.log",
    conda:
        "../envs/irma.yaml"
    threads:
        workflow.cores
    params:
        config="workflow/config/FLU-{order}-iri.sh"
    shell:
        "IRMA --external-config {params.config} FLU {input} {output} > {log}"
