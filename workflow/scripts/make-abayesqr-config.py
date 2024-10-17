#!/usr/bin/env python3

import argparse
from Bio.SeqIO import SeqIO

if __name__ == "__main__":

    parser = argparse.ArgumentParser("make-abayessqr-config.py")
    parser.add_argument("--fasta", help="Path to fasta reference file.")
    parser.add_argument("--sam", help="Path to .sam format aligned reads.")
    args = parser.parse_args()

    with open(args.fasta, "r") as fobj:
        record = next(SeqIO.parse(fobj, format="fasta"))

    config = f"""filename of reference sequence (FASTA) : {args.fasta}
    filname of the aligned reads (sam format) : {args.sam}
    paired-end (1 = true, 0 = false) : 1
    SNV_thres : 0.05
    reconstruction_start : 1
    reconstruction_stop: {len(record.seq)}
    min_mapping_qual : 60
    min_read_length : 100
    max_insert_length : 10
    characteristic zone name : abayesqr
    seq_err (assumed sequencing error rate(%)) : 0.1
    MEC improvement threshold : 0.0395
    """

    print(config)
