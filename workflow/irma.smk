from snakemake.utils import validate, min_version

min_version("7.0.4")

configfile: "irma-config.json"
validate(config, schema="schemas/irma-config-schema.json")



rule all:
    input:
        "results/xlsx/variants-mcc-by-sample-ordered.xlsx",
        "results/xlsx/variants-mcc-by-segment-ordered.xlsx",
        "results/xlsx/variants-mcc-flat-ordered.xlsx"


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
    conda:
        "envs/irma.yaml"
    shell:
        "IRMA FLU {input} {output} > {log}"


rule trim_trailing_tabs:
    input:
        "results/irma/{sample}_{pair}/tables/{table}.txt"
    output:
        "results/irma/{sample}_{pair}/tables/{table}.tsv"
    shell:
        "sed 's/\t$//g' < {input} > {output}"


rule write_gff:
    input:
        "results/irma/{sample}_{pair}/{segment}.fasta"
    output:
        gff="results/irma/{sample}_{pair}/{segment}.gff",
        gffgz="results/irma/{sample}_{pair}/{segment}.gff.gz"
    conda:
        "envs/tabix.yaml"
    log:
        "logs/write_gff/write_gff_{sample}_{pair}_{segment}.log"
    shell:
        """
        workflow/scripts/write-gff.py \
            --fasta-in {input} \
            --segment {wildcards.segment} \
            --transcript_id {wildcards.segment} 2> {log} > {output.gff} \
            --errors {config[errors]}

        # Append any warnings about length mismatches to IRMA reference in a log
        grep 'Length of FASTA' {log} >> logs/irma-ref-length-mismatches.log
        
        bgzip -c {output.gff} > {output.gffgz}
        tabix -p gff {output.gffgz}
        """


rule summarise_variants:
    input:
        vcf="results/irma/{sample}_{pair}/{segment}.vcf",
        fas="results/irma/{sample}_{pair}/{segment}.fasta",
        gff="results/irma/{sample}_{pair}/{segment}.gff.gz"
    output:
        "results/vep/{sample}_{pair}/{segment}.tsv"
    conda:
        "envs/vep.yaml"
    shell:
        """
        vep --input_file {input.vcf} \
            --output_file {output} \
            --gff {input.gff} \
            --fasta {input.fas} \
            --tab \
            --no_stats \
            --format vcf  # stops error with empty VCF files 

        # Remove upstream and downstream variants
        egrep -v "intron_variant|downstream_gene_variant|upstream_gene_variant" {output} > {output}.tmp
        mv {output}.tmp {output}
        """


rule merge_irma_vep:
    input:
        vep="results/vep/{sample}_{pair}/{segment}.tsv",
        irma_var="results/irma/{sample}_{pair}/tables/{segment}-variants.tsv",
        irma_ins="results/irma/{sample}_{pair}/tables/{segment}-insertions.tsv",
        irma_del="results/irma/{sample}_{pair}/tables/{segment}-deletions.tsv"
    output:
        "results/variants/{sample}_{pair}/{segment}.tsv"
    shell:
        """
        workflow/scripts/merge-vep-irma.py \
            --vep {input.vep} \
            --irma-var {input.irma_var} \
            --irma-del {input.irma_del} \
            --irma-ins {input.irma_ins} \
            --fasta-consensus results/irma/{wildcards.sample}_{wildcards.pair}/{wildcards.segment}.fasta > {output}
        """

rule multiple_changes_in_codon:
    input:
        "results/variants/{sample}_{pair}/{segment}.tsv"
    output:
        "results/variants-mcc/{sample}_{pair}/{segment}.tsv"
    shell:
        "workflow/scripts/multiple-changes-in-codon.py < {input} > {output}"


def aggregate_segments(wildcards):
    irma_out_dir = checkpoints.irma.get(**wildcards).output[0]
    segments = [file[:-4] for file in os.listdir(irma_out_dir) if file.endswith(".vcf")]
    return expand("results/variants-mcc/{sample}_{pair}/{segment}.tsv", segment=segments, **wildcards)


rule concat_segements:
    input:
        aggregate_segments
    output:
        "results/variants-mcc/{sample}_{pair}/{sample}_{pair}.tsv"
    shell:
        "workflow/scripts/concat-tables.py {input} > {output}"


rule combine_samples:
    input:
        expand(
            "results/variants-mcc/{sample}_{pair}/{sample}_{pair}.tsv",
            sample=config["samples"],
            pair=config["pair"]
        )
    output:
        "results/xlsx/variants-mcc-by-sample.xlsx"
    shell:
        "workflow/scripts/combine-tables.py {input} --excel {output}"

rule by_segment_summary:
    input:
        "results/xlsx/variants-mcc-by-sample.xlsx"
    output:
        segment="results/xlsx/variants-mcc-by-segment.xlsx",
        flat="results/xlsx/variants-mcc-flat.xlsx"
    shell:
        """
        workflow/scripts/make-by-segment-summary.py \
            --in-excel {input} \
            --out-segment {output.segment} \
            --out-flat {output.flat}
        """

rule order_columns:
    input:
        "{file}.xlsx"
    output:
        "{file}-ordered.xlsx"
    shell:
        """
        workflow/scripts/alter-column-order.py \
            --input {input} \
            --output {output} \
            --order \
                Sample \
                Variant	\
                Location \
                Segment \
                cDNA_position \
                Reference_Nuc_Position \
                Protein_position \
                Consequence \
                Amino_acids \
                Consensus_Amino_Acid \
                Minority_Amino_Acid \
                Consensus_Allele \
                Minority_Allele \
                Codons \
                Total_Reads \
                Consensus_Count \
                Minority_Count \
                Consensus_Frequency \
                Minority_Frequency \
                Consensus_Average_Quality \
                Minority_Average_Quality \
                ConfidenceNotMacErr \
                PairedUB \
                QualityUB \
                Phase \
                Mutation_Type \
                Codon_Position \
                Multiple_Changes_In_Codon
        """