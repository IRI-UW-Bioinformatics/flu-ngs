from pathlib import Path
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
        expand(
            "results/qsr/{order}/{pair}/{sample}/{sample}_{qsr_type}.fasta",
            pair=config["pair"],
            sample=config["samples"],
            qsr_type=config["qsr"],
            order=config["order"]
        )


wildcard_constraints:
    n="1|2",
    pair="(paired)|(combined)|(longread)",
    order="(primary)|(secondary)",
    qsr_type="(abayesqsr)|(tensqr)"


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


def collect_segments(path, default_wildcards=None):
    """
    Returns a function which make a list of files containing segment names, based on
    segments that IRMA has found.

    Args:
        path (str): What file names should look like. It should contain
            {segment} (which is expanded based on what segments IRMA finds), and can
            contain {sample} and {pair} (which are expanded based on wildcards).
        default_wildcards (dict): Optional default wildcardds to pass to expand. This is useful if
            the path that is passed doesn't contain wildcards that are necessary for finding the
            IRMA output.
    """
    default_wildcards = {} if default_wildcards is None else default_wildcards

    def fun(wildcards):
        """
        Make a list of files containing segment names, based on segments that IRMA has found.
        """
        irma_dir = checkpoints.find_irma_output.get(**wildcards, **default_wildcards).output[0]
        segments = [path.stem for path in Path(irma_dir).glob("*.vcf")]
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


rule combine_unpaired:
    input:
        "processed_reads/{sample}/{sample}_1_unpaired.fastq",
        "processed_reads/{sample}/{sample}_2_unpaired.fastq"
    output:
        temp("processed_reads/{sample}/{sample}_unpaired.fastq")
    shell:
        "cat {input} > {output}"


rule minimap_unpaired:
    """
    TenSQR doesn't run on .bam or .sam files generated by IRMA, so have to use minimap.

    Map the unpaired reads to the combined IRMA output (we don't generally run IRMA on just the
    unpaired reads).
    """
    input:
        "results/{order}/irma-raw/{sample}_combined/{segment}.fasta"
        "processed_reads/{sample}/{sample}_unpaired.fastq"
    output:
        temp("results/qsr/{order}/unpaired/{sample}/{segment}/aligned.sam")
    log:
        ".logs/qsr/minimap2_{order}_unpaired_{sample}_{segment}.txt"
    shell:
        "minimap2 -ax sr {input} > {output} 2> {log}"


rule minimap_paired:
    input:
        "results/{order}/irma-raw/{sample}_paired/{segment}.fasta"
        "processed_reads/{sample}/{sample}_1_paired.fastq",
        "processed_reads/{sample}/{sample}_2_paired.fastq"
    output:
        temp("results/qsr/{order}/paired/{sample}/{segment}/aligned.sam")
    log:
        ".logs/qsr/minimap2_{order}_paired_{sample}_{segment}.txt"
    shell:
        "minimap2 -ax sr {input} > {output} 2> {log}"


rule minimap_combined:
    input:
        "results/qsr/{order}/unpaired/{sample}/{segment}/aligned.sam"
        "results/qsr/{order}/paired/{sample}/{segment}/aligned.sam"
    output:
        temp("results/qsr/{order}/combined/{sample}/{segment}/aligned.sam")
    shell:
        "samtools merge -o {output} {input}"


rule make_qsr_config:
    input:
        fasta="results/{order}/irma/{sample}_{pair}/{segment}.fasta",
        sam="results/qsr/{order}/{pair}/{sample}/{segment}/aligned.sam"
    output:
        "results/qsr/{order}/{pair}/{sample}/{segment}/{qsr_type}_config.txt"
    shell:
        """
        workflow/scripts/make-qsr-config.py \
            --type {wildcards.qsr_type} \
            --fasta {input.fasta} \
            --sam aligned.sam \
            --prefix {wildcards.qsr_type} > {output}
        """


rule run_abayesqr:
    input:
        "results/qsr/{order}/{pair}/{sample}/{segment}/abayesqr_config.txt",
        "results/qsr/{order}/{pair}/{sample}/{segment}/aligned.sam"
    output:
        "results/qsr/{order}/{pair}/{sample}/{segment}/abayesqr_ViralSeq.fasta"
    params:
        working_dir="results/qsr/{order}/{pair}/{sample}/{segment}"
    shell:
        """
        cd {params.working_dir}
        aBayesQR abayesqr_config.txt > .abayesqr_log.txt 2>&1

        # Make the output of aBayesQR FASTA format
        # Add the segment this is from
        if [ -f abayesqr_ViralSeq.txt ]; then
            awk 'NR % 2 == 1 {{ print ">{wildcards.segment} " $0 }} NR % 2 == 0 {{ print $0 }}' \
                abayesqr_ViralSeq.txt > abayesqr_ViralSeq.fasta
        else
            touch abayesqr_ViralSeq.fasta
        fi
        """


rule run_tensqr:
    input:
        "results/qsr/{order}/{pair}/{sample}/{segment}/tensqr_config.txt",
        "results/qsr/{order}/{pair}/{sample}/{segment}/aligned.sam"
    output:
        "results/qsr/{order}/{pair}/{sample}/{segment}/tensqr_ViralSeq.fasta"
    params:
        working_dir="results/qsr/{order}/{pair}/{sample}/{segment}"
    shell:
        """
        cd {params.working_dir}
        ExtractMatrix tensqr_config.txt > .tensqr_extractmatrix_log.txt 2>&1

        if [ -s tensqr_SNV_matrix.txt ]; then

            # Try to run TenSQR and catch its exit code
            set +e  # Disable exiting on error
            TenSQR.py --zone_name tensqr > .tensqr_log.txt 2>&1
            exit_code=$?

            if [ $exit_code -ne 0 ]; then
                echo "TenSQR.py threw an error." >> .tensqr_log.txt
                echo "Making empty tensqr_ViralSeq.fasta." >> .tensqr_log.txt
                touch tensqr_ViralSeq.fasta

            elif [ ! -f tensqr_ViralSeq.fasta ]; then
                echo "TenSQR.py exited cleanly but didn't make tensqr_ViralSeq.fasta." >> .tensqr_log.txt
                echo "Making empty tensqr_ViralSeq.fasta." >> .tensqr_log.txt
                touch tensqr_ViralSeq.fasta

            fi

        else
            echo "No SNVs (ExtractMatrix produced empty tensqr_SNV_matrix.txt) so TenSQR wasn't run." >> .tensqr_log.txt
            echo "Making empty tensqr_ViralSeq.fasta." >> .tensqr_log.txt
            touch tensqr_ViralSeq.fasta

        fi
        """


rule collect_qsr_sequences_for_all_segments:
    input:
        collect_segments(
            "results/qsr/{order}/{pair}/{sample}/{segment}/{qsr_type}_ViralSeq.fasta",
        )
    output:
        "results/qsr/{order}/{sample}/{sample}_{qsr_type}.fasta"
    shell:
        "cat {input} > {output}"
