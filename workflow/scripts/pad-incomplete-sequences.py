#!/usr/bin/env python3

"""
Pad incomplete IRMA consensus sequences with N's to restore full-length
coordinates, and shift VCF / table positions accordingly.

Aligns each IRMA consensus to the IRMA reference (consensus.fasta) using a
semi-global pairwise alignment (no end-gap penalties on the reference). If
the consensus is shorter than the reference, the missing positions are padded
with N's at the start and/or end.
"""

import argparse
import sys
import warnings
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def load_references(reference_path):
    """Load reference FASTA into a dict keyed by sequence id."""
    refs = {}
    with open(reference_path) as fobj:
        for record in SeqIO.parse(fobj, "fasta"):
            refs[record.id] = record
    return refs


def compute_padding(query_seq, ref_seq):
    """
    Align *query_seq* (IRMA consensus) to *ref_seq* (IRMA reference) with a
    semi-global alignment (no end-gap penalties on the reference/target) and
    return (leading_ns, trailing_ns, alignment_str).
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    aligner.target_left_open_gap_score = 0
    aligner.target_left_extend_gap_score = 0
    aligner.target_right_open_gap_score = 0
    aligner.target_right_extend_gap_score = 0

    alignments = aligner.align(ref_seq, query_seq)
    aln = alignments[0]

    # aln.aligned gives two tuples-of-intervals: one for target (ref), one
    # for query.  The first interval of the *target* tells us where the query
    # starts matching on the reference.
    target_intervals = aln.aligned[0]  # intervals on reference
    leading_ns = target_intervals[0][0]
    trailing_ns = len(ref_seq) - target_intervals[-1][1]
    has_internal_gaps = len(target_intervals) > 1

    return leading_ns, trailing_ns, has_internal_gaps, format(aln)


# ---------------------------------------------------------------------------
# File-modification helpers
# ---------------------------------------------------------------------------

def pad_fasta(fasta_path, leading_ns, trailing_ns):
    """Prepend/append N's to the sequence in *fasta_path* (single record)."""
    with open(fasta_path) as fobj:
        record = next(SeqIO.parse(fobj, "fasta"))

    header = record.description
    new_seq = "N" * leading_ns + str(record.seq) + "N" * trailing_ns

    with open(fasta_path, "w") as fobj:
        fobj.write(f">{header}\n{new_seq}\n")


def shift_vcf(vcf_path, offset):
    """Add *offset* to the POS column (column 2, 1-based) of every data line."""
    lines = Path(vcf_path).read_text().splitlines(keepends=True)
    out = []
    for line in lines:
        if line.startswith("#"):
            out.append(line)
        else:
            cols = line.split("\t")
            cols[1] = str(int(cols[1]) + offset)
            out.append("\t".join(cols))
    Path(vcf_path).write_text("".join(out))


def shift_table(table_path, column_name, offset):
    """Add *offset* to *column_name* in a tab-separated IRMA table file."""
    lines = Path(table_path).read_text().splitlines(keepends=True)
    if not lines:
        return
    header = lines[0].rstrip("\n").split("\t")
    try:
        col_idx = header.index(column_name)
    except ValueError:
        return  # column not present – nothing to do
    out = [lines[0]]
    for line in lines[1:]:
        cols = line.split("\t")
        cols[col_idx] = str(int(cols[col_idx]) + offset)
        out.append("\t".join(cols))
    Path(table_path).write_text("".join(out))


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def pad_irma_dir(irma_dir, references, errors):
    """
    For every *.fasta in *irma_dir*, align to matching reference and pad if
    needed.  Also shift positions in corresponding VCF and table files.

    Writes incomplete-sequence-padding-report.txt to *irma_dir*.
    """
    irma_dir = Path(irma_dir)

    full_length_segments = []
    padded_segments = []

    for fasta_path in sorted(irma_dir.glob("*.fasta")):
        segment = fasta_path.stem

        if segment not in references:
            print(
                f"WARNING: no reference for segment '{segment}', skipping.",
                file=sys.stderr,
            )
            continue

        with open(fasta_path) as fobj:
            record = next(SeqIO.parse(fobj, "fasta"))

        ref_record = references[segment]
        leading_ns, trailing_ns, has_internal_gaps, alignment_str = compute_padding(
            str(record.seq), str(ref_record.seq)
        )

        if has_internal_gaps:
            raise ValueError(
                f"{segment}: alignment between consensus ({len(record.seq)} nt) "
                f"and reference ({len(ref_record.seq)} nt) contains internal "
                f"gaps. Padding cannot reliably adjust coordinates.\n\n"
                f"Alignment:\n{alignment_str}"
            )

        if leading_ns == 0 and trailing_ns == 0:
            full_length_segments.append(
                (segment, len(record.seq), len(ref_record.seq))
            )
            continue

        # --- warn / raise ------------------------------------------------
        msg = (
            f"{segment}: consensus length {len(record.seq)}, reference "
            f"length {len(ref_record.seq)}, padding {leading_ns} leading "
            f"and {trailing_ns} trailing N's"
        )
        print(msg, file=sys.stderr)

        if errors == "raise":
            raise ValueError(msg)

        warnings.warn(msg)

        padded_segments.append(
            (segment, len(record.seq), len(ref_record.seq),
             leading_ns, trailing_ns, alignment_str)
        )

        # --- pad FASTA ---------------------------------------------------
        pad_fasta(fasta_path, leading_ns, trailing_ns)

        # --- shift VCF ---------------------------------------------------
        if leading_ns > 0:
            vcf_path = irma_dir / f"{segment}.vcf"
            if vcf_path.exists():
                shift_vcf(vcf_path, leading_ns)

            # --- shift tables --------------------------------------------
            tables_dir = irma_dir / "tables"
            if tables_dir.is_dir():
                var_path = tables_dir / f"{segment}-variants.txt"
                if var_path.exists():
                    shift_table(var_path, "Position", leading_ns)

                ins_path = tables_dir / f"{segment}-insertions.txt"
                if ins_path.exists():
                    shift_table(ins_path, "Upstream_Position", leading_ns)

                del_path = tables_dir / f"{segment}-deletions.txt"
                if del_path.exists():
                    shift_table(del_path, "Upstream_Position", leading_ns)

    # --- write report ----------------------------------------------------
    write_report(irma_dir, full_length_segments, padded_segments)


def write_report(irma_dir, full_length_segments, padded_segments):
    """Write incomplete-sequence-padding-report.txt to *irma_dir*."""
    report_path = irma_dir / "incomplete-sequence-padding-report.txt"
    with open(report_path, "w") as f:
        f.write("Incomplete-sequence padding report\n")
        f.write("=" * 60 + "\n\n")

        # --- full-length segments ----------------------------------------
        f.write("Segments that did not need padding\n")
        f.write("-" * 60 + "\n")
        f.write("(Consensus already matches reference length — no action taken.)\n\n")
        if full_length_segments:
            for segment, query_len, ref_len in full_length_segments:
                f.write(f"  {segment}: {query_len} nt (reference {ref_len} nt)\n")
        else:
            f.write("  (none)\n")

        f.write("\n")

        # --- padded segments ---------------------------------------------
        f.write("Segments that were padded\n")
        f.write("-" * 60 + "\n\n")
        if not padded_segments:
            f.write("  (none)\n")
        else:
            for (segment, query_len, ref_len,
                 leading, trailing, alignment_str) in padded_segments:
                f.write(f"  {segment}\n")
                f.write(f"    Consensus length: {query_len}\n")
                f.write(f"    Reference length: {ref_len}\n")
                f.write(f"    Leading N's:  {leading}\n")
                f.write(f"    Trailing N's: {trailing}\n")
                f.write(f"\n    Alignment (reference on top, consensus on bottom):\n\n")
                aln_lines = alignment_str.splitlines()
                width = 80
                for start in range(0, len(aln_lines[0]), width):
                    for line in aln_lines:
                        f.write(f"      {line[start:start+width]}\n")
                    f.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "pad-incomplete-sequences.py",
        description=(
            "Pad incomplete IRMA consensus sequences with N's and shift "
            "VCF / table positions to match."
        ),
    )
    parser.add_argument(
        "--irma-dir",
        required=True,
        help="Path to an IRMA output directory (modified in place).",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Path to reference FASTA (workflow/reference/consensus.fasta).",
    )
    parser.add_argument(
        "--errors",
        default="warn",
        choices=("warn", "raise"),
        help="How to handle incomplete sequences: 'warn' or 'raise'.",
    )
    args = parser.parse_args()

    references = load_references(args.reference)
    pad_irma_dir(args.irma_dir, references, args.errors)
