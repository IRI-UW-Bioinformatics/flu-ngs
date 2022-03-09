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
            ignored for the purpose of writing the GFF file. The --fasta file is still
            required to check that the sequence is the correct length for the predefined
            GFF file.
            """,
    )
    parser.add_argument("--fasta", help="FASTA file.", required=True)
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

    with open(args.fasta) as fobj:
        record = next(SeqIO.parse(fobj, format="fasta"))

    known_length = {"A_MP": 982, "A_PA": 2151, "A_PB1": 2274}

    if args.segment in known_length:

        # Lengths of IRMA reference segments GFF files define where the splice sites are.
        # These were written w.r.t the IRMA reference sequence. So, if the consensus fasta
        # sequence differs in length to the reference, then the splice sites defined in
        # the GFF probably won't be correct. (Even if they are the same length they
        # _might_ not be correct...)

        if known_length[args.segment] != len(record):
            raise ValueError(
                """
                Length of FASTA ({}) differs from IRMA reference ({}) used to write splice
                positions defined in GFF file.
                """.format(
                    len(record), expect[args.segment]
                )
            )

        path = os.path.join("workflow", "gff", "{}.gff".format(args.segment))

        with open(path) as fobj:
            gff = fobj.read()

    else:

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
