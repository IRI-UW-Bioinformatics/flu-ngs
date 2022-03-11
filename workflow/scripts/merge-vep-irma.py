#!/usr/bin/env python3

import sys
import argparse

from Bio import SeqIO
import pandas as pd


def lookup_deleted_nucs(df_irma_del, seq):
    """
    Lookup deleted nucleotides to add to the IRMA deletions DataFrame.
    """
    nucs = []
    for _, row in df_irma_del.iterrows():
        pos = row["Upstream_Position"]
        length = row["Length"]
        nucs.append("".join(seq[pos : pos + length]))
    df_irma_del["Consensus_Allele"] = nucs
    df_irma_del["Minority_Allele"] = "-"
    return df_irma_del


def make_index_for_irma_variants(df):
    """
    Construct an index from the IRMA {segment}-variants.txt columns.
    """
    df.index = (
        df["Reference_Name"]
        + "_"
        + df["Position"].astype(str)
        + "_"
        + df["Consensus_Allele"]
        + "/"
        + df["Minority_Allele"]
    )
    return df


def make_index_for_irma_insertions(df):
    """
    Construct an index from the IRMA {segment}-insertions.txt columns.
    """
    df.index = (
        df["Reference_Name"]
        + "_"
        + (df["Upstream_Position"] + 1).astype(str)
        + "_-/"
        + df["Insert"]
    )
    return df


def make_index_for_irma_deletions(df):
    """
    Construct an index from the IRMA {segment}-deletions.txt columns.
    """
    df.index = (
        df["Reference_Name"]
        + "_"
        + (df["Upstream_Position"] + 1).astype(str)
        + "_"
        + df["Consensus_Allele"]
        + "/-"
    )
    return df


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        "merge-vep-irma.py",
        description="Merge outputs of VEP and IRMA.",
    )
    parser.add_argument("--vep", help="VEP output file", required=True)
    parser.add_argument(
        "--irma-var",
        help="IRMA {segment}-variants.txt. Located in 'tables' IRMA output directory.",
        required=True,
    )
    parser.add_argument(
        "--irma-ins",
        help="IRMA {segment}-insertions.txt. Located in 'tables' IRMA output directory.",
        required=True,
    )
    parser.add_argument(
        "--irma-del",
        help="IRMA {segment}-deletions.txt. Located in 'tables' IRMA output directory.",
        required=True,
    )
    parser.add_argument(
        "--fasta-consensus", help="FASTA consensus sequence found by IRMA.", required=True
    )
    parser.add_argument(
        "--verbose_output",
        help="Don't remove columns from output that contain superfluous or duplicated "
        "information.",
        action="store_true",
        default=False,
    )

    args = parser.parse_args()

    with open(args.fasta_consensus) as fobj:
        record = next(SeqIO.parse(fobj, format="fasta"))

    df_vep = pd.read_table(
        args.vep,
        comment="##",
        engine="python",
        index_col=0,
    )

    df_irma_var = pd.read_table(args.irma_var).pipe(make_index_for_irma_variants)

    df_irma_ins = (
        pd.read_table(args.irma_ins)
        .eval("Consensus_Allele = '-'")
        .pipe(make_index_for_irma_insertions)
        .rename(
            columns={
                "Frequency": "Minority_Frequency",
                "Average_Quality": "Minority_Average_Quality",
                "Count": "Minority_Count",
                "Insert": "Minority_Allele",
                "Upstream_Position": "Position",
            }
        )
        .drop(columns=["Context", "Called"])
    )

    df_irma_del = (
        pd.read_table(args.irma_del)
        .pipe(lookup_deleted_nucs, seq=record)
        .pipe(make_index_for_irma_deletions)
        .rename(
            columns={
                "Frequency": "Minority_Frequency",
                "Count": "Minority_Count",
                "Upstream_Position": "Position",
            }
        )
        .drop(columns=["Context", "Called", "Length"])
    )

    df_irma = pd.concat([df_irma_var, df_irma_ins, df_irma_del])

    df_out = (
        df_vep.join(df_irma)
        .sort_values(["Feature", "Position"])
        .rename(
            columns={
                "Feature": "Segment",
                "Total": "Total_Reads",
                "Position": "Reference_Nuc_Position",
            }
        )
    )

    if not args.verbose_output:
        gz_column = next(
            column for column in df_vep.columns if column.endswith(".gff.gz")
        )
        df_out = df_out.drop(
            [
                "Allele",
                "CDS_position",
                "DISTANCE",
                "Existing_variation",
                "Feature_type",
                "FLAGS",
                "Gene",
                "IMPACT",
                "SOURCE",
                "STRAND",
                "Reference_Name",
                gz_column,
            ],
            axis=1,
        )

    df_out.to_csv(sys.stdout, sep="\t")
