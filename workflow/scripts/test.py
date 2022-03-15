#!/usr/bin/env python3

import unittest
import importlib
from Bio import SeqIO

# Using importlib to handle '-' in .py names
wg = importlib.import_module("write-gff")


class TestFindAll(unittest.TestCase):
    """
    Tests for write_gff.findall.
    """

    def test_case_a(self):
        """
        Simple test case.
        """
        self.assertEqual([0, 4], wg.findall("d", "david"))

    def test_finds_overlapping(self):
        """
        Should return indexes of matches even if they overlap.
        """
        self.assertEqual([0, 2], wg.findall("ada", "adada"))


class TestFindNsSpliceDonor(unittest.TestCase):
    """
    Tests for write_gff.find_ns_splice_donor.
    """

    def test_case_a(self):
        """
        Test case with known results.
        """
        seq = """agcaaaagcagggtgacaaagacataatggattccaacacagtgtcaagctttcaggtagattgttttctttggcacattcgcaaacgatttgcagaccaaaaaatgggtgatgccccgtttcttgaccgacttcgcagagatcaaaagtccttaaaaggaagaagcagcactcttggtctagacattgaaagctcaacactagcagggaggcaaatagtaaagcggattctaaaagaagaatctgatagtgaaccaaaagggactattacctcagtacccacttcatattatttaactgacatgactcttgaagaaatgtcaagggcctggttcatgcttatacccaaccaaaaaagagtaggatcactctgcatcagaatggatcaagccataatggataaggaaatcacactgaaggcaaactttagtgtggtctttaacaaactggagactctaacacttttacgagcattcacagatgatgaagcaattattggagaaatcttaccaataccttctcttccaggacatactaacgaggatgtcaaaaatgcaattgagatcctcatcggaggacttgaatggaataataacacagttcgaatctctgagattctacagagattcacttggagaaacagtaatgagaatgggggatttttactctctccaaaacaaaaacaaaaaatggagggaacaactgggccagaagtttgaagaaataagatggctgaggtattgaagaaataaggcataaactaaaaataacagagaacagttttgaacaaataacatttattcaagcattacaactattgcttgaagtggagcaagagataagaactttctcgtttcagcttatttaatgataaaaaacacccttgtttctact"""
        self.assertEqual(54, wg.find_ns_splice_donor(seq))


class TestFindNsSpliceAcceptor(unittest.TestCase):
    """
    Tests for write_gff.find_ns_splice_acceptor.
    """

    def test_known_cases(self):
        """
        Check the correct splice acceptors are found for all sequences in
        ../test-data/ns-seqs.fasta.
        """
        with open("../test-data/ns-seqs.fas") as fobj:
            records = {
                record.description: record.seq for record in SeqIO.parse(fobj, "fasta")
            }

        # Expected values were identified by looking at an alignement of these sequences,
        # finding the splice acceptor site, and then subtracting the number of empty
        # columns (due to indels) in the alignment
        expect = {
            "A/mallard/Alberta/827/1978 (Allele B)": 526 - 18,
            "A/pintail duck/Alberta/121/1979 (Allele B)": 526 - 8,
            "A/chicken/Germany/N/49 (Allele B)": 526 - 13,
            "A/duck/Hokkaido/167/2007 (Allele B)": 526 - 8,
            "A/chicken/Italy/4746/1999 (Allele B)": 526 - 15,
            "A/equine/Prague/1/56": 526 - 0,
            "A/swine/Italy/1850/77": 526 - 7,
            "A/Brevig Mission/1/18": 526 - 26,
            "A/equine/New Market/1979": 526 - 26,
            "A/Alaska/232/2015": 526 - 26,
            "A/California/07/2009": 526 - 26,
            "A/Viet Nam/1203/2004": 526 - 26 - 15,
            "A/Indonesia/5/2005": 526 - 0 - 15,
            "A/Hubei/1/2010": 526 - 21 - 15,
        }

        for desc, seq in records.items():
            with self.subTest(desc=desc):
                self.assertEqual(expect[desc], wg.find_ns_splice_acceptor(seq))


class TestExtract12NtsAroundSpliceSite(unittest.TestCase):
    """
    Tests for write_gff.extract_12_nts_around_splice_site.
    """

    def test_case_a(self):
        """
        A simple test case.
        """
        output = wg.extract_12_nts_around_splice_site("XTTTCAGGTAAGAASDFASDF", 5)
        self.assertEqual("TTTCAGGTAAGA", output)


class TestSpliceNS(unittest.TestCase):
    """
    Tests for write_gff.splice_ns.
    """

    def test_raises_value_error_incorrect_donor_loc(self):
        """
        A ValueError should be raised if asked to splice a sequence without an "AGGT" at
        the donor location.
        """
        seq = "AGCAAAAGCAGGGTGACAAAGACATAATGGATTCCAACACAGTGTCAAGCTTTCAGGTAGATTGTTTTCTTTGGCACATTCGCAAACGATTTGCAGACCAAAAAATGGGTGATGCCCCGTTTCTTGACCGACTTCGCAGAGATCAAAAGTCCTTAAAAGGAAGAAGCAGCACTCTTGGTCTAGACATTGAAAGCTCAACACTAGCAGGGAGGCAAATAGTAAAGCGGATTCTAAAAGAAGAATCTGATAGTGAACCAAAAGGGACTATTACCTCAGTACCCACTTCATATTATTTAACTGACATGACTCTTGAAGAAATGTCAAGGGCCTGGTTCATGCTTATACCCAACCAAAAAAGAGTAGGATCACTCTGCATCAGAATGGATCAAGCCATAATGGATAAGGAAATCACACTGAAGGCAAACTTTAGTGTGGTCTTTAACAAACTGGAGACTCTAACACTTTTACGAGCATTCACAGATGATGAAGCAATTATTGGAGAAATCTTACCAATACCTTCTCTTCCAGGACATACTAACGAGGATGTCAAAAATGCAATTGAGATCCTCATCGGAGGACTTGAATGGAATAATAACACAGTTCGAATCTCTGAGATTCTACAGAGATTCACTTGGAGAAACAGTAATGAGAATGGGGGATTTTTACTCTCTCCAAAACAAAAACAAAAAATGGAGGGAACAACTGGGCCAGAAGTTTGAAGAAATAAGATGGCTGAGGTATTGAAGAAATAAGGCATAAACTAAAAATAACAGAGAACAGTTTTGAACAAATAACATTTATTCAAGCATTACAACTATTGCTTGAAGTGGAGCAAGAGATAAGAACTTTCTCGTTTCAGCTTATTTAATGATAAAAAACACCCTTGTTTCTACT"
        donor_loc = 15
        with self.assertRaisesRegex(
            ValueError, f"No AGGT at position {donor_loc} in {seq}"
        ):
            wg.splice_ns(seq, donor_loc, 400)

    def test_raises_value_error_incorrect_accept_loc(self):
        """
        A ValueError should be raised if asked to splice a sequence without an "AG" at
        the acceptor location.
        """
        seq = "AGCAAAAGCAGGGTGAGGTACAAAGACATAATGGATTCCAACACAGTGTCAAGCTTTCAGGTAGATTGTTTTCTTTGGCACATTCGCAAACGATTTGCAGACCAAAAAATGGGTGATGCCCCGTTTCTTGACCGACTTCGCAGAGATCAAAAGTCCTTAAAAGGAAGAAGCAGCACTCTTGGTCTAGACATTGAAAGCTCAACACTAGCAGGGAGGCAAATAGTAAAGCGGATTCTAAAAGAAGAATCTGATAGTGAACCAAAAGGGACTATTACCTCAGTACCCACTTCATATTATTTAACTGACATGACTCTTGAAGAAATGTCAAGGGCCTGGTTCATGCTTATACCCAACCAAAAAAGAGTAGGATCACTCTGCATCAGAATGGATCAAGCCATAATGGATAAGGAAATCACACTGAAGGCAAACTTTAGTGTGGTCTTTAACAAACTGGAGACTCTAACACTTTTACGAGCATTCACAGATGATGAAGCAATTATTGGAGAAATCTTACCAATACCTTCTCTTCCAGGACATACTAACGAGGATGTCAAAAATGCAATTGAGATCCTCATCGGAGGACTTGAATGGAATAATAACACAGTTCGAATCTCTGAGATTCTACAGAGATTCACTTGGAGAAACAGTAATGAGAATGGGGGATTTTTACTCTCTCCAAAACAAAAACAAAAAATGGAGGGAACAACTGGGCCAGAAGTTTGAAGAAATAAGATGGCTGAGGTATTGAAGAAATAAGGCATAAACTAAAAATAACAGAGAACAGTTTTGAACAAATAACATTTATTCAAGCATTACAACTATTGCTTGAAGTGGAGCAAGAGATAAGAACTTTCTCGTTTCAGCTTATTTAATGATAAAAAACACCCTTGTTTCTACT"
        donor_loc = 15
        accept_loc = 400
        with self.assertRaisesRegex(
            ValueError, f"No AG at position {accept_loc} in {seq}"
        ):
            wg.splice_ns(seq, donor_loc, accept_loc)

    def test_raises_value_error_donor_loc_gt_accept_loc(self):
        """
        Should raise a ValueError if the donor_loc is within 350 nts of the accept_loc.
        """
        seq = "CAGTACCCACTTCATATTATTTAACTGACATGACTCTTGAAGA"
        donor_loc = 15
        accept_loc = donor_loc + 349
        msg = f"Splice acceptor signal location \({accept_loc}\) should be at least 350 nts downstream of the donor signal \({donor_loc}\) location, but it is {accept_loc - donor_loc}"
        with self.assertRaisesRegex(ValueError, msg):
            wg.splice_ns(seq, donor_loc, accept_loc)


if __name__ == "__main__":
    unittest.main()
