#!/usr/bin/env python3

import sys
import pandas as pd
from Bio import Seq


def index_upper(chars) -> int:
    """Return index of the first uppercase character in characters"""
    for i, char in enumerate(chars):
        if char.isupper():
            return i


def aggregate_multiple_changes_in_codon(df) -> pd.DataFrame:
    """
    Take a DataFrame containing rows of merged VEP and IRMA output where alleles have the
    same "Segment", "Protein_position" and "Phase". I.e. they are from the same codon.
    Combine these changes into a single row which summarises the effect on the single
    codon.

    Returns a new DataFrame containing a single row.
    """
    name = ", ".join(df.index)

    # Check the DataFrame really only contains data from one codon
    for column in "Segment Protein_position Phase".split():
        unique_values = df[column].unique()
        if len(unique_values) > 1:
            raise ValueError(
                f"""A DataFrame was passed to aggregate_multiple_changes_in_codon which
                 has multiple values in the {column} column. A DataFrame containing
                 these alleles was passed: {name}"""
            )

    # Mini dataframe containing the consensus and minority codons
    df_samecodon = (
        df["Codons"]
        .str.split("/", expand=True)
        .rename(columns={0: "Consensus_Codon", 1: "Minority_Codon"})
    )
    df = df.join(df_samecodon)

    # Check consensus codons match
    unique_consensus_codons = df["Consensus_Codon"].str.upper().unique()
    if len(unique_consensus_codons) != 1:
        raise ValueError(
            f"Something has gone wrong. There shouldn't be more than 1 consensus codon at a "
            f"single protein position! Issue occurred at: {name}"
        )

    # Position in codon of the mutation
    df["Codon_Position"] = df["Minority_Codon"].apply(index_upper)

    # Build the minority codon. Start with the consensus, then update sites for all
    # minority alleles. Start with lower case. If a position is mutated change it to
    # uppercase in consensus and minority to be consistent with VEP output
    consensus_codon = list(unique_consensus_codons[0].lower())
    minority_codon = list(consensus_codon)
    for _, row in df.iterrows():
        i = row["Codon_Position"]
        minority_codon[i] = row["Minority_Allele"].upper()
        consensus_codon[i] = consensus_codon[i].upper()

    consensus_codon = "".join(consensus_codon)
    minority_codon = "".join(minority_codon)

    # Translate codons -> amino acids
    consensus_aa = Seq.Seq(consensus_codon).translate()
    minority_aa = Seq.Seq(minority_codon).translate()

    # 1-indexed locations for consistency with rest of data
    df["Codon_Position"] += 1

    consequence = (
        "synonymous_variant" if consensus_aa == minority_aa else "missense_variant"
    )

    def join_values_in_column(col):
        return ", ".join(map(str, sorted(df[col])))

    return pd.DataFrame(
        {
            **{
                col: join_values_in_column(col)
                for col in set(df.columns) - {"Consensus_Codon", "Minority_Codon"}
            },
            **{
                "Codons": f"{consensus_codon}/{minority_codon}",
                "Consequence": consequence,
                "Consensus_Amino_Acid": consensus_aa,
                "Minority_Amino_Acid": minority_aa,
                "Multiple_Changes_In_Codon": "Yes",
                "Consensus_Allele": consensus_codon,
                "Minority_Allele": minority_codon,
                "Phase": df["Phase"][0],
                "Segment": df["Segment"][0],
                "Protein_position": df["Protein_position"][0],
            },
        },
        index=[name],
    )


if __name__ == "__main__":

    df = pd.read_table(sys.stdin, index_col=0)

    if df.empty:
        df.to_csv(sys.stdout, sep="\t")
    else:
        rows = [
            sub if len(sub) == 1 else aggregate_multiple_changes_in_codon(sub)
            for _, sub in df.groupby(["Segment", "Protein_position", "Phase"])
        ]
        (
            pd.concat(rows)
            .eval("Protein_position = Protein_position.astype('int')")
            .sort_values(["Segment", "Protein_position"])
            .to_csv(sys.stdout, sep="\t")
        )
