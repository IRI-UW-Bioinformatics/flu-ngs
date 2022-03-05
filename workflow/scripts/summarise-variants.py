#!/usr/bin/env python3

import sys
from typing import Dict, List

import vcf
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


def load_fasta(path) -> List:
    """Load FASTA records.

    Args:
        path (str): Path to FASTA file.
    """
    with open(path) as fobj:
        return list(SeqIO.parse(fobj, format="fasta"))


def load_vcf(path) -> List:
    """Load VCF records.

    Args:
        path (str): Path to VCF file.
    """
    return [record for record in vcf.Reader(open(path))]


class SummariseVariants:
    def __init__(self, reference, variants):
        self.reference = reference
        self.variants = variants

        for variant in self.variants:
            self.check_variant_consistent_with_reference(variant)

    def check_variant_consistent_with_reference(self, variant):
        """Check the value the variant lists as reference matches that in the reference
        sequence.

        Args:
            variant (vcf.model._Record)
            reference (str)

        Raises:
            ValueError if there is a mismatch.
        """
        index = variant.POS - 1
        ref_length = len(variant.REF)
        fasta_ref = self.reference[index : index + ref_length]
        if fasta_ref != variant.REF:
            raise ValueError(
                "VCF lists {ref} at position {pos} but fasta file has {fas}".format(
                    ref=variant.REF, pos=variant.POS, fas=fasta_ref
                )
            )

    @property
    def dataframe(self) -> pd.DataFrame:
        """
        Summarise variants in a pandas DataFrame.
        """

        data = []

        for variant in self.variants:

            if not isinstance(variant.REF, str):
                raise TypeError("Expected variant.REF to be a string: {}".format(variant))

            if not isinstance(variant.ALT, list):
                raise TypeError("Expected variant.ALT to be a list: {}".format(variant))

            if not len(variant.ALT) == 1:
                raise TypeError(
                    "Expected variant.ALT list to have one element: {}".format(variant)
                )

            # Single nucleotide mutation
            if len(variant.REF) == 1 and len(variant.ALT[0]) == 1:
                data.append(self.summarise_single_variants(variant))

            # Deletion
            elif len(variant.REF) > len(variant.ALT[0]):
                data.append(self.summarise_deletion(variant))

            # Insertion
            elif len(variant.REF) < len(variant.ALT[0]):
                data.append(self.summarise_insertion(variant))

            else:
                raise NotImplementedError(
                    "This variant has multiple changes: {}\nTalk to david.".format(
                        variant
                    )
                )

        return pd.DataFrame(data).astype(
            {"Codon pos": "Int32", "Nucleotide in codon pos": "Int32"}
        )

    def summarise_single_variants(self, variant) -> Dict:
        """Summarise a single nucleotide mutation variant.

        Args:
            variant (vcf.model._Record)
        """
        nuc_index = variant.POS - 1

        codon_index = nuc_index // 3
        codon_start = codon_index * 3
        codon_end = codon_start + 3
        codon_ref = self.reference[codon_start:codon_end]

        # Position of the variant in the codon [0, 1, 2]
        pos_in_codon = nuc_index % 3
        codon_var = "".join(
            [
                nuc if i != pos_in_codon else str(variant.ALT[0])
                for (i, nuc) in enumerate(codon_ref)
            ]
        )

        aa_ref = Seq(codon_ref).translate()
        aa_var = Seq(codon_var).translate()

        effect = "Synonymous" if aa_ref == aa_var else "Non-synonymous"

        return {
            "Effect": effect,
            "Nucleotide change": "{} -> {}".format(variant.REF, variant.ALT[0]),
            "Amino acid change": "{} -> {}".format(aa_ref, aa_var),
            "Codon change": "{} -> {}".format(codon_ref, codon_var),
            "Nucleotide pos": nuc_index + 1,
            "Codon pos": codon_index + 1,
            "Nucleotide in codon pos": pos_in_codon + 1,
            "Allele frequency": variant.INFO["AF"][0],
        }

    def summarise_deletion(self, variant) -> Dict:
        """Summarise a deletion.

        Args:
            variant (vcf.model._Record)
        """
        return {
            "Effect": "Deletion",
            "Nucleotide change": "{} -> {}".format(variant.REF, variant.ALT[0]),
            "Nucleotide pos": variant.POS,
            "Allele frequency": variant.INFO["AF"][0],
        }

    def summarise_insertion(self, variant) -> Dict:
        """Summarise an insertion.

        Args:
            variant (vcf.model._Record)
        """
        return {
            "Effect": "Insertion",
            "Nucleotide change": "{} -> {}".format(variant.REF, variant.ALT[0]),
            "Nucleotide pos": variant.POS,
            "Allele frequency": variant.INFO["AF"][0],
        }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        "summarise-variants.py",
        description="""
            Summarise variants in a .vcf file, looking up coding changes relative to a
            reference sequence. The program assumes that the reference sequence is a
            single CDS, with the first codon starting at the first position in the
            sequence. A summary table is printed to stdout.
            """,
    )
    parser.add_argument("--vcf", help="Path to .vcf file", required=True)
    parser.add_argument(
        "--reference",
        help="""
            Path to fasta format reference sequence. The first sequence in the file is
            used if multiple are present.
            """,
        required=True,
    )
    args = parser.parse_args()

    reference = str(load_fasta(args.reference)[0].seq)
    variants = load_vcf(args.vcf)
    sv = SummariseVariants(reference, variants)
    sv.dataframe.to_csv(sys.stdout, index=False)
