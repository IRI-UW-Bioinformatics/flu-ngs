"""
Running IRMA on MiSeq data.
"""


rule irma_raw:
    input:
        expand("processed_reads/{{sample}}/{{sample}}_{n}_{{pair}}.fastq", n=(1, 2)),
    output:
        directory("results/{order}/irma-raw/{sample}_{pair}"),
    log:
        ".logs/irma-{order}-raw/{sample}_{pair}.log",
    conda:
        "../envs/irma.yaml"
    threads:
        4  # Feel free to set higher if you don't have many samples
    shell:
        """
        IRMA \
            --external-config workflow/config/FLU-{wildcards.order}-iri.sh \
            FLU {input} {output} > {log}
        """
