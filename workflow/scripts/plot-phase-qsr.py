#!/usr/bin/env python3

from collections import defaultdict
from pathlib import Path
import argparse
import re


from Bio.SeqIO import parse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use("seaborn-v0_8-dark")

SEGMENT_TO_PROTEIN = {"MP": "M1", "NS": "NS1"}
SEGMENTS = "PB2", "PB1", "PA", "NP", "HA", "NA", "MP", "NS"


def snp_position(snp: str) -> int:
    """
    Lookup the position of a SNP.
    """
    try:
        return int(re.search("(\d+)", snp).group())
    except AttributeError as err:
        print(f"'{snp}' doesn't contain an integer")
        raise err


def infer_snps(seq_a, seq_b, ignore="-"):
    ignore = set(ignore)

    if len(seq_a) != len(seq_b):
        raise ValueError("lengths don't match")

    for site, (a, b) in enumerate(zip(seq_a, seq_b), start=1):
        if a != b and a not in ignore and b not in ignore:
            yield f"{a}{site}{b}"


def snp_columns_to_df(series: "pd.Series[list[str]]"):
    return (
        pd.DataFrame(
            {strain_i: {snp: 1 for snp in snps} for strain_i, snps in series.items()}
        )
        .T.fillna(0)
        .astype(int)
    )


def records_to_df(records: "list[dict]") -> pd.DataFrame:
    df = pd.DataFrame(records).set_index("strain_i")
    return (
        df.join(snp_columns_to_df(df["snps"]))
        .drop(columns="snps")
        .set_index("freq", append=True)
    )


def clean_segment(segment: str) -> str:
    segment = re.sub("^A_", "", segment)
    segment = re.sub("_H5$", "", segment)
    segment = re.sub("_N1", "", segment)
    return segment


def add_snp_column(df: pd.DataFrame) -> pd.DataFrame:
    df["is_snp"] = (df["Consensus_Allele"].str.len() == 1) & (
        df["Minority_Allele"].str.len() == 1
    )
    return df


def make_snp(row: pd.Series) -> str:
    c = row["Consensus_Allele"]
    m = row["Minority_Allele"]

    if not isinstance(c, str) or not isinstance(m, str):
        return None

    if len(c) != 1 or len(m) != 1:
        return None

    site = int(row["Reference_Nuc_Position"])
    return f"{m}{site}{c}"


def load_irma_tables(paths: "list[str]") -> pd.DataFrame | None:
    """
    Load IRMA variant tables for multiple segments. These are output by the workflow in:

        results/primary/variants/{sample}/{segment}.tsv

    Args:
        paths: List of paths for a single sample.

    Returns:
        pd.DataFrame
    """
    tables = [pd.read_table(table) for table in paths]
    tables = [table for table in tables if not table.empty]
    if tables:
        df = pd.concat(tables).pipe(add_snp_column).query("is_snp")

        df["Segment"] = (
            df["Segment"]
            .str.replace("^A_", "", regex=True)
            .str.replace("_H5$", "", regex=True)
            .str.replace("_N1", "", regex=True)
        )

        df["SNP"] = df.apply(make_snp, axis=1)

        return df.set_index("Segment")


def load_qsr_sequences(paths: "list[str]", ref_seq_dir: str) -> dict:
    """
    Load QSR sequences for multiple segments for a single sample.

    Args:
        paths: List of *_ViralSeq_tensqr.fasta file paths from a single sample.
        ref_seq_dir: Directory containing reference sequences used in QSR inference.

    Returns:
        dict
    """
    qsr = defaultdict(dict)

    for fasta in paths:
        *_, segment_long, _ = Path(fasta).parts

        with open(Path(ref_seq_dir, f"{segment_long}.fasta")) as fobj:
            ref = str(next(parse(fobj, "fasta")).seq)

        seg = clean_segment(segment_long)

        records = []
        with open(fasta) as fobj:
            for record in parse(fobj, "fasta"):
                strain_i, _, freq = record.description.split()
                records.append(
                    {
                        "strain_i": int(strain_i.split("_")[-1]),
                        "freq": float(freq),
                        "snps": list(infer_snps(str(record.seq), ref, ignore="-N*")),
                    }
                )

        if records:
            qs_snp_df = records_to_df(records)
            qsr[seg] = {"qs": records, "ref": ref, "qs_snp_df": qs_snp_df}

    return qsr


def snp_plot(qsr: dict[str, dict], df: pd.DataFrame):
    """
    Plots SNP comparison between IRMA variants and QSR reconstructions.

    Args:
        qsr (dict[str, dict]): A dictionary containing quasispecies sequence data for each segment.
        df (pd.DataFrame): A DataFrame containing IRMA SNPs.
    """
    # Count how many SNPs each segment has (ax widths are scaled accordingly)
    snps = defaultdict(set)
    df_irma = {}
    df_qsr = {}

    for seg in SEGMENTS:

        irma_seg = SEGMENT_TO_PROTEIN.get(seg, seg)
        if irma_seg in df.index:
            # Using a list as key here means a DataFrame is returned (rather than a series) even if
            # just a single row is present.
            df_irma[seg] = df.loc[[irma_seg]]
            snps[seg].update(df_irma[seg]["SNP"])

        if seg in qsr:
            df_qsr[seg] = qsr[seg]["qs_snp_df"]
            snps[seg].update(df_qsr[seg].columns)

        snps[seg] = sorted(snps[seg], key=snp_position)

    _, axes = plt.subplots(
        ncols=len(SEGMENTS),
        sharey=True,
        gridspec_kw=dict(
            width_ratios=[
                len(snps[seg]) if len(snps[seg]) > 1 else 1 for seg in SEGMENTS
            ]
        ),
        figsize=(15, 2.5),
    )

    for seg, ax in zip(SEGMENTS, axes):

        summary = []

        # IRMA SNPs coloured by phase
        if seg in df_irma:
            for phase, sub_irma_df in df_irma[seg].groupby("Phase"):
                x = [snps[seg].index(snp) for snp in sub_irma_df["SNP"]]
                y = sub_irma_df["Minority_Frequency"]
                ax.scatter(x, y, zorder=15, label=phase, clip_on=False)
        else:
            summary.append("No IRMA SNPs")

        # Quasispecies sequence SNP frequencies
        if seg in qsr:
            for qs_strain in qsr[seg]["qs"]:
                freq = qs_strain["freq"]
                ax.axhline(freq, linewidth=0.5, c="grey", zorder=10)
                x = [snps[seg].index(snp) for snp in qs_strain["snps"]]
                y = np.repeat(freq, len(x))
                ax.scatter(
                    x,
                    y,
                    marker="+",
                    zorder=15,
                    s=75,
                    clip_on=False,
                    c="black",
                    linewidth=0.75,
                )
        else:
            summary.append("No QSR")

        ax.text(
            0,
            1,
            "\n".join(summary),
            ha="left",
            va="top",
            c="grey",
            transform=ax.transAxes,
            fontsize=6,
        )
        ax.set_xticks(
            np.arange(len(snps[seg])), snps[seg], rotation=90, family="monospace"
        )
        ax.set(ylabel="Frequency", title=seg, ylim=(0, 1))

        for spine in "top", "bottom", "right":
            ax.spines[spine].set_visible(False)
        ax.label_outer()


if __name__ == "__main__":
    parser = argparse.ArgumentParser("plot-phase-qsr.py")
    parser.add_argument(
        "--irma-tables",
        nargs="+",
        help=(
            "Variant files generated by IRMA for a single sample. By default these are put in: "
            "results/primary/variants/{sample}/{segment}.tsv."
        ),
    )
    parser.add_argument(
        "--qsr-fastas",
        nargs="+",
        help=(
            "QSR fasta files generated by TenSQR for a single sample. By default the workflow puts "
            "these in: results/qsr/{order}/{sample}/{segment}/tensqr_ViralSeq.fasta"
        ),
    )
    parser.add_argument(
        "--ref-seq-dir",
        required=True,
        help=(
            "Directory containing reference sequences that were used in QSR inference. By default "
            "the workflow uses results/primary/irma/{sample}_combined."
        ),
    )
    parser.add_argument("--title", required=False, help="Title for the plot")
    parser.add_argument(
        "--output_name",
        help="Name of output file(s). Pass whatever extensions understood by matplotlib.",
        nargs="+",
    )
    args = parser.parse_args()

    df_irma = load_irma_tables(args.irma_tables)
    qsr = load_qsr_sequences(paths=args.qsr_fastas, ref_seq_dir=args.ref_seq_dir)

    snp_plot(qsr, df_irma)

    if args.title:
        plt.suptitle(args.title, x=0.5, y=1.1)

    for output in args.output_name:
        plt.savefig(output, bbox_inches="tight")

    plt.close()
