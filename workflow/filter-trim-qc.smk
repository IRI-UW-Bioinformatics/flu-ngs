from snakemake.utils import validate, min_version

min_version("7.0.4")

configfile: "trim-qc-config.json"
validate(config, schema="schemas/trim-qc-config-schema.json")

rule all:
    input: 
        expand("combined/{sample}/{sample}_combined.fastq.gz", sample=config["samples"]),
        expand("combined/{sample}/{sample}_filtered.fastq.gz", sample=config["samples"]), 

rule combine_minion_reads:            
    output:
        "combined/{sample}/{sample}_combined.fastq.gz"
    shell:         
        "cat raw/{wildcards.sample}/*fastq.gz > {output}"
     
rule filter_minion_reads:
    input:
        "combined/{sample}/{sample}_combined.fastq.gz"
    output:
        "combined/{sample}/{sample}_filtered.fastq.gz"
    shell:    
        "gunzip -c {input} | chopper -q 10 --minlength 600 --maxlength 2500 | gzip > {output}"
 
 
