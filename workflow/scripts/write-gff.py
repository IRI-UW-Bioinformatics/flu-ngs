#!/usr/bin/env python3

import os
from sys import stderr, stdout
import argparse
from Bio import SeqIO


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "write-gff.py",
        description="""
            Writes a simple GFF file of a protein coding sequence in a FASTA file. The
            program assumes the FASTA file contains a single CDS which starts at the first
            nucleotide in the FASTA file, and ends at the last. Only the first sequence in
            the FASTA file is used.

            Alternatively, pass either A_MP, A_PA or A_PB1 to the --segment keyword
            argument to write a predefined GFF file. If one of these segments is passed to
            this flag then the predefined GFF file is written, and the --fasta argument is
            ignored.
            """,
    )
    parser.add_argument("--fasta", help="FASTA file.", required=False)
    parser.add_argument(
        "--segment", help="A segment with a predefined GFF file.", required=False
    )
    parser.add_argument(
        "--transcript_id",
        help="ID for the transcript. (Gets put in the 'Feature' column in VEP output.)",
        default="transcript1",
        required=False,
    )
    args = parser.parse_args()

    if args.segment in {"A_MP", "A_PA", "A_PB1"}:
        path = os.path.join("workflow", "gff", "{}.gff".format(args.segment))
        with open(path) as fobj:
            gff = fobj.read()

    else:
        with open(args.fasta) as fobj:
            record = next(SeqIO.parse(fobj, format="fasta"))

        path = os.path.join("workflow", "gff", "generic_template.txt")
        with open(path) as fobj:
            template = fobj.read()
            gff = template.format(
                name=record.description,
                start=1,
                end=len(record),
                transcript_id=args.transcript_id,
            )

    stderr.write("Used {} to write GFF\n".format(path))
    stdout.write(gff)
