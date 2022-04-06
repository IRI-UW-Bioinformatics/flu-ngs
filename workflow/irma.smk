from snakemake.utils import validate, min_version

min_version("7.0.4")

configfile: "irma-config.json"
validate(config, schema="schemas/irma-config-schema.json")



rule all:
    input:
        "results/xlsx/variants-mcc-by-sample-ordered.xlsx",
        "results/xlsx/variants-mcc-by-segment-ordered.xlsx",
        "results/xlsx/variants-mcc-flat-ordered.xlsx",
        expand(
            "results/seq/{sample}_{pair}/aa.fasta",
            sample=config["samples"],
            pair=config["pair"]
        ),
        expand(
            "results/seq/{sample}_{pair}/nt.fasta",
            sample=config["samples"],
            pair=config["pair"]
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
        "logs/irma/irma_{sample}_{pair}.log"
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
        "results/irma/{sample}_{pair}/{segment}.gff"
    log:
        "logs/write_gff/write_gff_{sample}_{pair}_{segment}.log"
    shell:
        """
        workflow/scripts/write-gff.py \
            --fasta-in {input} \
            --segment {wildcards.segment} \
            --transcript_id {wildcards.segment} \
            --errors {config[errors]} 2> {log} > {output}
        """


rule make_gffgz:
    input:
        "results/irma/{sample}_{pair}/{segment}.gff"
    output:
        "results/irma/{sample}_{pair}/{segment}.gff.gz"
    log:
        "logs/write_gffgz/write_gffgz_{sample}_{pair}_{segment}.log"
    conda:
        "envs/tabix.yaml"
    shell:
        """
        bgzip -c {input} > {output} 2> {log}
        tabix -p gff {output} 2> {log}
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


def collect_segments(path):
    """
    A function that returns a function which make a list of files containing
    segment names, based on segments that IRMA has found.

    Args:
        path (str): What file names should look like. It should contain
            {segment} (which is expanded based on what segments IRMA finds), and can
            contain {sample} and {pair} (which are expanded based on wildcards).
    """

    def fun(wildcards):
        """
        Make a list of files containing segment names, based on segments that IRMA has found.
        """
        irma_dir = checkpoints.irma.get(**wildcards).output[0]
        segments = [
            file[:-4]
            for file in os.listdir(irma_dir)
            if file.endswith(".vcf")
        ]
        return expand(path, segment=segments, **wildcards)

    return fun


rule concat_segment_aa:
    input:
        collect_segments("results/seq/{sample}_{pair}/separate/{segment}-aa.fasta")
    output:
        "results/seq/{sample}_{pair}/aa.fasta"
    shell:
        "cat {input} > {output}"


rule concat_segment_nt:
    input:
        collect_segments("results/seq/{sample}_{pair}/separate/{segment}-nt.fasta")
    output:
        "results/seq/{sample}_{pair}/nt.fasta"
    shell:
        "cat {input} > {output}"


rule transcribe:
    input:
        fasta="results/irma/{sample}_{pair}/{segment}.fasta",
        gff="results/irma/{sample}_{pair}/{segment}.gff"
    output:
        "results/seq/{sample}_{pair}/separate/{segment}-nt.fasta"
    conda:
        "envs/gffread.yaml"
    log:
        "logs/gffread/gffread_{sample}_{pair}_{segment}.log"
    shell:
        "gffread -w {output} -g {input.fasta} {input.gff} > {log} 2>&1"


rule translate:
    input:
        "results/seq/{sample}_{pair}/separate/{segment}-nt.fasta"
    output:
        "results/seq/{sample}_{pair}/separate/{segment}-aa.fasta"
    conda:
        "envs/emboss.yaml"
    log:
        "logs/transeq/transeq_{sample}_{pair}_{segment}.log"
    shell:
        "transeq -sequence {input} -outseq {output} > {log} 2>&1"


rule concat_segements:
    input:
        collect_segments("results/variants-mcc/{sample}_{pair}/{segment}.tsv")
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
    log:
        "logs/make-segment-summary.log"
    shell:
        """
        workflow/scripts/make-by-segment-summary.py \
            --in-excel {input} \
            --out-segment {output.segment} \
            --out-flat {output.flat} > {log} 2>&1

        # grep's exit status is 1 if it doesn't find any matches, causing snakemake to throw an error
        # set +e, set -e prevents this
        set +e
        grep 'Length of consensus found by IRMA ' logs/write_gff/*.log > logs/incorrect-splice-vars.log
        set -e
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