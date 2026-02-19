#!/usr/bin/env python3
"""
Collect FASTA sequences from irma-raw and irma directories for a given segment,
writing one combined file with headers containing sample name, segment, and source.
"""

import argparse
from pathlib import Path


def read_fasta_sequence(fasta_path):
    """Read a single-sequence FASTA file and return the sequence."""
    lines = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                lines.append(line.strip())
    return "".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--irma-raw", nargs="+", required=True,
                        help="irma-raw FASTA files")
    parser.add_argument("--irma-padded", nargs="+", required=True,
                        help="irma (padded) FASTA files")
    parser.add_argument("--segment", required=True)
    args = parser.parse_args()

    for source_label, paths in [("irma-raw", args.irma_raw),
                                ("irma-padded", args.irma_padded)]:
        for fasta_path in paths:
            fasta_path = Path(fasta_path)
            # Extract sample_pair from path: results/{order}/{irma-raw|irma}/{sample_pair}/{segment}.fasta
            sample_pair = fasta_path.parent.name
            seq = read_fasta_sequence(fasta_path)
            print(f">{sample_pair} {args.segment} {source_label}")
            print(seq)


if __name__ == "__main__":
    main()
