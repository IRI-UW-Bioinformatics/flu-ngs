#!/usr/bin/env python3

from pathlib import Path
import argparse
from Bio.SeqIO import parse

if __name__ == "__main__":

    parser = argparse.ArgumentParser("make-abayesqr-config.py")
    parser.add_argument(
        "--fasta",
        help="Path to fasta reference file. N.B. only the filename will be included in the config "
        "- aBayesQR must be called from the same directory as the data.",
        required=True,
    )
    parser.add_argument(
        "--sam", help="Path to .sam format aligned reads.", required=True
    )
    args = parser.parse_args()

    with open(args.fasta, "r") as fobj:
        record = next(parse(fobj, format="fasta"))

    filename = Path(args.fasta).name

    config = f"""filename of reference sequence (FASTA) : {filename}
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
