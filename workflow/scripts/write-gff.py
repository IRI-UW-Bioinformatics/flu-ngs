#!/usr/bin/env python3

import os
from sys import stderr, stdout
import argparse
from typing import Tuple, List
import warnings

from Bio import SeqIO


def findall(sub: str, string: str) -> List[int]:
    """
    Return indexes of all substrings in a string
    """
    indexes = []
    l = len(sub)
    for i in range(len(string)):
        if string[i : i + l] == sub:
            indexes.append(i)
    return indexes


def find_ns_splice_donor(seq: str) -> int:
    """
    Find the AGGT splice donor signal. Returns an int which is the index of the 'A'

    Notes:
        In NS1 sequences Gabi has sent the AGGT is on average at position 38.2.
    """
    seq = seq.upper()
    candidates = findall("AGGT", seq)

    if not candidates:
        raise ValueError(f"No NS1 splice donor site ('AGGT') in {seq}")

    # Splice site should be in the first 100 nucs
    candidates = filter(lambda x: x < 100, candidates)

    if not candidates:
        raise ValueError(f"No NS1 splice donor sites ('AGGT') in first 100 nt of {seq}")

    # Choose the site that is closest to position 38.2
    return min(candidates, key=lambda x: abs(x - 38.2))


def find_ns_splice_acceptor(seq: str, donor_loc=None) -> int:
    """
    Find the 'AG' splice acceptor signal. Returns an int which is the index of the 'A'.

    Notes:
        Should be >350 nts downstread of the splice donor location.

    Args:
        seq (str)
        donor_loc (int): Location of the splice donor
    """
    seq = seq.upper()

    donor_loc = find_ns_splice_donor(seq) if donor_loc is None else donor_loc

    # Find all candidate acceptor sites 350 nts downstream of the donor site
    candidates = findall("AG", seq[donor_loc + 350 :])
    candidates = [c + donor_loc + 350 for c in candidates]  # Fix indexing

    # sequence around the splice site should be FQDI
    candidates = list(
        filter(
            lambda x: four_aas_around_splice_site(seq, donor_loc, x) == "FQDI",
            candidates,
        )
    )

    if not candidates:
        raise ValueError(f"No NS splice acceptor sites in {seq}")

    # Pick candidate that is closest to position 507
    return min(candidates, key=lambda x: abs(x - 507))


def splice_ns(seq: str, donor_loc: int, accept_loc: int) -> str:
    """
    Splice an NS sequence given a splice donor and acceptor locations.

    Args:
        seq (str):
        donor_loc (int): Location of the AGGT. The 'AG' remains in the
            transcript, the 'GT' is lost.
        acceptor_loc (int): Location of the 'AG'. The 'AG' is lost.
    """
    seq = seq.upper()

    if (accept_loc - donor_loc) < 350:
        raise ValueError(
            f"Splice acceptor signal location ({accept_loc}) should be at least 350 nts "
            f"downstream of the donor signal ({donor_loc}) location, but it is "
            f"{accept_loc - donor_loc}"
        )

    if seq[donor_loc : donor_loc + 4] != "AGGT":
        raise ValueError(f"No AGGT at position {donor_loc} in {seq}")

    if seq[accept_loc : accept_loc + 2] != "AG":
        raise ValueError(f"No AG at position {accept_loc} in {seq}")

    return seq[: donor_loc + 2] + seq[accept_loc + 2 :]


def four_aas_around_splice_site(seq: str, donor_loc: int, accept_loc: int) -> str:
    """
    What are the four amino acids either side of the splice site given a sequence,
    donor location and acceptor location?
    """
    spliced = splice_ns(seq, donor_loc, accept_loc)
    return extract_12_nts_around_splice_site(spliced, donor_loc).translate()


def extract_12_nts_around_splice_site(seq: str, donor_loc: int) -> str:
    """
    Start 4 nts downstream of the donor location

            Start of splice donor 'AGGT' signal
            ⌄     Splice site
            |     ⌄
        XXXXTTTCAG|GA...
        ^
        Extracts 12 nts from here
    """
    start = donor_loc - 4
    end = start + 12
    return seq[start:end]


def find_ns_splice_sites(seq: str) -> Tuple[int, int]:
    """
    Lookup the splice donor and acceptor locations for an NS1 transcript.
    """
    donor_loc = find_ns_splice_donor(seq)
    accept_loc = find_ns_splice_acceptor(seq, donor_loc)
    return donor_loc, accept_loc


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
    parser.add_argument("--fasta-in", help="FASTA file.", required=True)
    parser.add_argument(
        "--segment", help="A segment with a predefined GFF file.", required=False
    )
    parser.add_argument(
        "--transcript_id",
        help="ID for the transcript. (Gets put in the 'Feature' column in VEP output.)",
        default="transcript1",
        required=False,
    )
    parser.add_argument(
        "--errors",
        help="""
            How to handle errors. 'warn' issues a warning, but tries to carry on. 'raise'
            raises an error and stops.
            """,
        default="warn",
        choices=("warn", "raise"),
    )
    args = parser.parse_args()

    with open(args.fasta_in) as fobj:
        record = next(SeqIO.parse(fobj, format="fasta"))

    known_length = {"A_MP": 982, "A_PA": 2151, "A_PB1": 2274}

    if args.segment in known_length:

        # Lengths of IRMA reference segments GFF files define where the splice sites are.
        # These were written w.r.t the IRMA reference sequence. So, if the consensus fasta
        # sequence differs in length to the reference, then the splice sites defined in
        # the GFF probably won't be correct. (Even if they are the same length they
        # _might_ not be correct...)

        if known_length[args.segment] != len(record):

            msg = (
                "Length of consensus found by IRMA ({}) differs from length of "
                "reference ({}) used to write splice positions defined in GFF "
                "file.".format(len(record), known_length[args.segment])
            )

            if args.errors == "warn":
                warnings.warn(msg)
            elif args.errors == "raise":
                raise ValueError(msg)
            else:
                raise ValueError("'errors' should be 'warn' or 'raise'")

        path = os.path.join("workflow", "gff", "{}.gff".format(args.segment))

        with open(path) as fobj:
            gff = fobj.read()

    elif args.segment == "A_NS":
        path = os.path.join("workflow", "gff", "A_NS_template.txt")

        with open(path) as fobj:
            template = fobj.read()

        donor_loc, accept_loc = find_ns_splice_sites(record.seq)

        # NS1 and NS2 have 193 and 337 additional nts after splice acceptor (AG)
        # respectively
        ns1_end = accept_loc + 193
        ns2_end = accept_loc + 337

        # donor_loc corresponds to start of AGGT signal. 'GT' is trimmed, leaving 'AG'.
        # So, the position of the first G is the end of the first exon.
        ns2_exon1_end = donor_loc + 1

        # accept_loc is the start of the AG splice acceptor signal. In splicing the AG is
        # lost. So, need the position of the next nucleotide after the AG.
        ns2_exon2_start = accept_loc + 2

        # +1 for 1-based indexing used in GFF files
        gff = template.format(
            ns1_end=ns1_end + 1,
            ns2_end=ns2_end + 1,
            ns2_exon1_end=ns2_exon1_end + 1,
            ns2_exon2_start=ns2_exon2_start + 1,
        )

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
