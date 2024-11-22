#!/usr/bin/env python3

from pathlib import Path
import argparse
from Bio.SeqIO import parse

ABAYESQR_TEMPLATE = """filename of reference sequence (FASTA) : {fasta}
filname of the aligned reads (sam format) : {sam}
paired-end (1 = true, 0 = false) : {paired_end}
SNV_thres : {snv_thresh}
reconstruction_start : 1
reconstruction_stop: {length}
min_mapping_qual : {min_mapping_qual}
min_read_length : {min_read_length}
max_insert_length : {max_insert_length}
characteristic zone name : {prefix}
seq_err (assumed sequencing error rate(%)) : {seq_err}
MEC improvement threshold : 0.0395"""

TENSQR_TEMPLATE = """filename of reference sequence (FASTA) : {fasta}
filname of the aligned reads (sam format) : {sam}
SNV_thres : {snv_thresh}
reconstruction_start : 1
reconstruction_stop: {length}
min_mapping_qual : {min_mapping_qual}
min_read_length : {min_read_length}
max_insert_length : {max_insert_length}
characteristic zone name : {prefix}
seq_err (assumed sequencing error rate(%)) : {seq_err}
MEC improvement threshold : 0.0312
initial population size : 5"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "make-qsr-config.py",
        description="Make a config file for running a Quasispecies Spectrum Reconstruction (QSR) "
        "program (TenSQR or aBayesQR).",
    )
    parser.add_argument("--type", choices=["abayesqr", "tensqr"], required=True)
    parser.add_argument(
        "--fasta",
        help=(
            "Path to fasta reference file. Only the filename is included in the config (TenSQR "
            "and aBayesQR both dump output in the directory they're called from). The fasta file "
            "is required to read the length of the reference that is passed into the config. The "
            "first record in this file is used as the reference."
        ),
        required=True,
    )
    parser.add_argument(
        "--sam", help="Path to .sam format aligned reads.", required=True
    )
    parser.add_argument(
        "--prefix", help="Prefix appended to output files.", required=True
    )
    parser.add_argument(
        "--paired_end",
        help="0/1 for False/True. Only used by aBayesQR. Default=1.",
        default=1,
    )
    parser.add_argument(
        "--min_read_length", help="Default=125", default=125, required=False
    )
    parser.add_argument(
        "--max_insert_length",
        help="Maximum length of insertions allowed. Default=100",
        default=100,
        required=False,
    )
    parser.add_argument(
        "--seq_err",
        default=0.2,
        help="Assumed sequencing error rate (percent) Default=0.2.",
    )
    parser.add_argument(
        "--min_mapping_qual",
        default=30,
        help="Minimum read mapping quality. Default=30.",
    )
    parser.add_argument(
        "--snv_thresh",
        default=0.01,
        help=(
            "Single nucleotide variant threshold. I don't know exactly how this parameter is used "
            "by aBayesQR or TenSQR. Perhaps SNVs below this threshold are not considered when "
            "reconstructing haplotypes?"
        ),
    )
    args = parser.parse_args()

    with open(args.fasta, "r") as fobj:
        record = next(parse(fobj, format="fasta"))

    template = {"abayesqr": ABAYESQR_TEMPLATE, "tensqr": TENSQR_TEMPLATE}[args.type]

    kwds = vars(args)

    # Remove fasta path from kwds and pass just the filename to the template
    fasta_path = kwds.pop("fasta")
    fasta_name = Path(fasta_path).name
    print(
        template.format(length=len(record.seq), fasta=fasta_name, **kwds),
        end="",
    )
