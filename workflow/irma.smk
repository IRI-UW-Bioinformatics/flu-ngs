from typing import Callable
from snakemake.utils import validate, min_version

min_version("7.0.4")


configfile: "config.json"


validate(config, schema="schemas/config-schema.json")


def expand_order(path):
    "Helper function to call expand with config order."
    return expand(path, order=config["order"])


def expand_sample_pair_order(path):
    "Helper function to call expand with config sample, pair and order."
    return expand(
        path, sample=config["samples"], pair=config["pair"], order=config["order"]
    )


rule all:
    input:
        expand_order("results/{order}/xlsx/variants-mcc-by-sample-ordered.xlsx"),
        expand_order("results/{order}/xlsx/variants-mcc-by-segment-ordered.xlsx"),
        expand_order("results/{order}/xlsx/variants-mcc-flat-ordered.xlsx"),
        expand_sample_pair_order("results/{order}/seq/{sample}_{pair}/aa.fasta"),
        expand_sample_pair_order("results/{order}/seq/{sample}_{pair}/nt.fasta"),


wildcard_constraints:
    n="1|2",
    pair="(paired)|(combined)|(longread)",
    order="(primary)|(secondary)",


include: f"rules/irma_{config['platform']}.smk"


checkpoint find_irma_output:
    """
    IRMA puts output from secondary assemblies in a directory called:

        secondary_assembly

    Sometimes, however it puts it in

        secondary_assembly/secondary_assembly

    Other times it puts it in:

        secondary_assembly/secondary_assembly/secondary_assembly

    Here, find the most nested secondary_assembly path and link it in to the
    results directory so that it is easy to point to for other rules.
    """
    input:
        "results/{order}/irma-raw/{sample}_{pair}",
    output:
        directory("results/{order}/irma/{sample}_{pair}"),
    log:
        ".logs/irma-{order}/{sample}_{pair}.log",
    shell:
        """
        # Find the most nested secondary_assembly dir
        DIR="$(find {input} -name secondary_assembly | sort | tail -n 1)" > {log}

        # Some samples may not trigger secondary_assembly ($DIR will be empty for these)
        # For these cases link the regular IRMA output
        [ -z "$DIR" ] && DIR={input} >> {log}

        # Make directory if necessary
        [ ! -d results/{wildcards.order}/irma ] && mkdir results/{wildcards.order}/irma >> {log}

        # Finally, make the link
        ln -s "../../../$DIR" {output} >> {log}
        """


rule trim_trailing_tabs:
    """
    Some output tables from IRMA has trailing tabs which generates errors in
    pandas.
    """
    input:
        "results/{order}/irma/{sample}_{pair}/tables/{table}.txt",
    output:
        "results/{order}/irma/{sample}_{pair}/tables/{table}.tsv",
    shell:
        "sed 's/\t$//g' < {input} > {output}"


rule write_gff:
    input:
        "results/{order}/irma/{sample}_{pair}/{segment}.fasta",
    output:
        "results/{order}/irma/{sample}_{pair}/{segment}.gff",
    log:
        ".logs/write_gff/write_gff_{sample}_{pair}_{segment}_{order}.log",
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
        "results/{order}/irma/{sample}_{pair}/{segment}.gff",
    output:
        "results/{order}/irma/{sample}_{pair}/{segment}.gff.gz",
    log:
        ".logs/write_gffgz/write_gffgz_{sample}_{pair}_{segment}_{order}.log",
    conda:
        "envs/tabix.yaml"
    shell:
        """
        bgzip -c {input} > {output} 2> {log}
        tabix -p gff {output} 2> {log}
        """


rule summarise_variants:
    input:
        vcf="results/{order}/irma/{sample}_{pair}/{segment}.vcf",
        fas="results/{order}/irma/{sample}_{pair}/{segment}.fasta",
        gff="results/{order}/irma/{sample}_{pair}/{segment}.gff.gz",
    output:
        "results/{order}/vep/{sample}_{pair}/{segment}.tsv",
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
        vep="results/{order}/vep/{sample}_{pair}/{segment}.tsv",
        irma_var="results/{order}/irma/{sample}_{pair}/tables/{segment}-variants.tsv",
        irma_ins="results/{order}/irma/{sample}_{pair}/tables/{segment}-insertions.tsv",
        irma_del="results/{order}/irma/{sample}_{pair}/tables/{segment}-deletions.tsv",
    output:
        "results/{order}/variants/{sample}_{pair}/{segment}.tsv",
    shell:
        """
        workflow/scripts/merge-vep-irma.py \
            --vep {input.vep} \
            --irma-var {input.irma_var} \
            --irma-del {input.irma_del} \
            --irma-ins {input.irma_ins} \
            --fasta-consensus results/{wildcards.order}/irma/{wildcards.sample}_{wildcards.pair}/{wildcards.segment}.fasta > {output}
        """


rule multiple_changes_in_codon:
    input:
        "results/{order}/variants/{sample}_{pair}/{segment}.tsv",
    output:
        "results/{order}/variants-mcc/{sample}_{pair}/{segment}.tsv",
    shell:
        "workflow/scripts/multiple-changes-in-codon.py < {input} > {output}"


def collect_segments(path):
    """
    Returns a function which make a list of files containing segment names, based on
    segments that IRMA has found.

    Args:
        path (str): What file names should look like. It should contain
            {segment} (which is expanded based on what segments IRMA finds), and can
            contain {sample} and {pair} (which are expanded based on wildcards).
    """

    def fun(wildcards):
        """
        Make a list of files containing segment names, based on segments that IRMA has found.
        """
        irma_dir = checkpoints.find_irma_output.get(**wildcards).output[0]
        segments = [file[:-4] for file in os.listdir(irma_dir) if file.endswith(".vcf")]
        return expand(path, segment=segments, **wildcards)

    return fun


rule concat_segment_aa:
    input:
        collect_segments(
            "results/{order}/seq/{sample}_{pair}/separate/{segment}-aa.fasta"
        ),
    output:
        "results/{order}/seq/{sample}_{pair}/aa.fasta",
    shell:
        "cat {input} > {output}"


rule concat_segment_nt:
    input:
        collect_segments(
            "results/{order}/seq/{sample}_{pair}/separate/{segment}-nt.fasta"
        ),
    output:
        "results/{order}/seq/{sample}_{pair}/nt.fasta",
    shell:
        "cat {input} > {output}"


rule transcribe:
    input:
        fasta="results/{order}/irma/{sample}_{pair}/{segment}.fasta",
        gff="results/{order}/irma/{sample}_{pair}/{segment}.gff",
    output:
        "results/{order}/seq/{sample}_{pair}/separate/{segment}-nt.fasta",
    conda:
        "envs/gffread.yaml"
    log:
        ".logs/gffread/gffread_{sample}_{pair}_{segment}_{order}.log",
    shell:
        "gffread -w {output} -g {input.fasta} {input.gff} > {log} 2>&1"


rule translate:
    input:
        "results/{order}/seq/{sample}_{pair}/separate/{segment}-nt.fasta",
    output:
        "results/{order}/seq/{sample}_{pair}/separate/{segment}-aa.fasta",
    conda:
        "envs/emboss.yaml"
    log:
        ".logs/transeq/transeq_{sample}_{pair}_{segment}_{order}.log",
    shell:
        "transeq -sequence {input} -outseq {output} > {log} 2>&1"


rule concat_segements:
    input:
        collect_segments("results/{order}/variants-mcc/{sample}_{pair}/{segment}.tsv"),
    output:
        "results/{order}/variants-mcc/{sample}_{pair}/{sample}_{pair}.tsv",
    shell:
        "workflow/scripts/concat-tables.py {input} > {output}"


rule combine_samples:
    input:
        expand(
            "results/{{order}}/variants-mcc/{sample}_{pair}/{sample}_{pair}.tsv",
            sample=config["samples"],
            pair=config["pair"],
        ),
    output:
        temp("results/{order}/xlsx/variants-mcc-by-sample.xlsx"),
    shell:
        "workflow/scripts/combine-tables.py {input} --excel {output}"


rule by_segment_summary:
    input:
        "results/{order}/xlsx/variants-mcc-by-sample.xlsx",
    output:
        segment=temp("results/{order}/xlsx/variants-mcc-by-segment.xlsx"),
        flat=temp("results/{order}/xlsx/variants-mcc-flat.xlsx"),
    log:
        ".logs/make-segment-summary-{order}.log",
    shell:
        """
        workflow/scripts/make-by-segment-summary.py \
            --in-excel {input} \
            --out-segment {output.segment} \
            --out-flat {output.flat} > {log} 2>&1

        # grep's exit status is 1 if it doesn't find any matches, causing snakemake to throw an error
        # set +e, set -e prevents this
        set +e
        grep 'Length of consensus found by IRMA ' .logs/write_gff/*.log > .logs/incorrect-splice-vars-{wildcards.order}.log
        set -e
        """


rule order_columns:
    input:
        "{file}.xlsx",
    output:
        "{file}-ordered.xlsx",
    shell:
        """
        workflow/scripts/alter-column-order.py \
            --input {input} \
            --output {output} \
            --order \
                Sample \
                Variant    \
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


rule trim_fastq:
    input:
        "raw/{sample}/{sample}_{n}.fastq"
    output:
        [
            f".processed_reads_qsr/{{sample}}/{{sample}}_{n}_{pair}.fastq"
            for n in (1, 2)
            for pair in ("paired", "unpaired")
        ]
    log:
        ".logs/qsr/trim_{wildcards.sample}_{wildcards.segment}.txt"
    shell:
        "TrimmomaticPE {input} {output} ILLUMINACLIP:raw/trimlog.fas:2:30:10:2 MINLEN:36 2> {log}"


rule align_unfiltered_to_segment:
    input:
        "results/primary/seq/{sample}_paired/separate/{segment}-nt.fasta",
        ".processed_reads_qsr/{sample}/{sample}_1_paired.fastq",
        ".processed_reads_qsr/{sample}/{sample}_2_paired.fastq"
    output:
        "results/qsr/{sample}/{segment}/aligned.sam"
    log:
        ".logs/qsr/minimap2_{wildcards.sample}_{wildcards.segment}.txt"
    shell:
        "minimap2 -ax sr {input} > {output} 2> {log}"


rule make_abayesqr_config:
    input:
        fasta="results/primary/seq/{sample}_paired/separate/{segment}-nt.fasta",
        sam="results/qsr/{sample}/{segment}/aligned.sam"
    output:
        "results/qsr/{sample}/{segment}/abayesqr_config.txt"
    shell:
        "workflow/scripts/make-abayesqr-config.py --fasta {input.fasta} --sam {input.sam} > {output}"


rule abayes_qsr:
    input:
        make_abayesqr_config.output
    output:
        "results/qsr/{sample}/{segment}/abayesqr_ViralSeq.fasta"
    params:
        working_dir="results/qsr/{wildcards.sample}/{wildcards.segment}"
    log:
        ".logs/qsr/aBayesQr_{wildcards.sample}_{wildcards.segment}.txt"
    shell:
        """
        cd {params.working_dir}
        aBayesQR abayesqr_config.txt > {log} 2>&1

        # Make the output of aBayesQR FASTA format
        awk 'NR % 2 == 1 { print ">" $0 } NR % 2 == 0 { print $0 }' abayesqr_ViralSeq.txt \
            > abayesqr_ViralSeq.fasta
        """