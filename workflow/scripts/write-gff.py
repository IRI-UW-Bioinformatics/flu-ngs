#!/usr/bin/env python3

import argparse
from Bio import SeqIO


TEMPLATE = """{name}	flu-ngs	gene	{start}	{end}	.	+	.	ID=gene1;Name=GENE1
{name}	flu-ngs	transcript	{start}	{end}	.	+	.	ID=transcript1;Name=GENE1-001;Parent=gene1;biotype=protein_coding
{name}	flu-ngs	exon	{start}	{end}	.	+	.	ID=exon1;Name=GENE1-001_1;Parent=transcript1
{name}	flu-ngs	CDS	{start}	{end}	.	+	0	ID=cds1;Name=CDS0001;Parent=transcript1"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "write-gff.py",
        description="""
            Writes a simple GFF file of a protein coding sequence in a FASTA file. The
            program assumes the FASTA file contains a single CDS which starts at the first
            nucleotide in the FASTA file, and ends at the last. Only the first sequence in
            the FASTA file is used.
            """,
    )
    parser.add_argument("--fasta", help="FASTA file")
    args = parser.parse_args()

    with open(args.fasta) as fobj:
        record = next(SeqIO.parse(fobj, format="fasta"))

    print(TEMPLATE.format(name=record.description, start=1, end=len(record)))
