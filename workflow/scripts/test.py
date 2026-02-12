#!/usr/bin/env python3

import unittest
import importlib
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Using importlib to handle '-' in .py names
wg = importlib.import_module("write-gff")
mvi = importlib.import_module("merge-vep-irma")
pis = importlib.import_module("pad-incomplete-sequences")


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


class TestValidateGffCoordinates(unittest.TestCase):
    """
    Tests for write_gff.validate_gff_coordinates.
    """

    def test_valid_gff_returns_empty(self):
        gff = (
            "A_MP\tflu-ngs\tgene\t1\t982\t.\t+\t.\tID=gene1\n"
            "A_MP\tflu-ngs\texon\t1\t819\t.\t+\t.\tID=exon1_1;Parent=A_M1\n"
            "A_MP\tflu-ngs\texon\t715\t982\t.\t+\t.\tID=exon3;Parent=A_M2"
        )
        self.assertFalse([], wg.validate_gff_coordinates(gff))

    def test_detects_start_greater_than_end(self):
        gff = (
            "A_MP\tflu-ngs\tgene\t1\t452\t.\t+\t.\tID=gene1\n"
            "A_MP\tflu-ngs\texon\t1\t26\t.\t+\t.\tID=exon2;Parent=A_M2\n"
            "A_MP\tflu-ngs\texon\t715\t452\t.\t+\t.\tID=exon3;Parent=A_M2"
        )
        invalid = wg.validate_gff_coordinates(gff)
        self.assertEqual(len(invalid), 1)
        self.assertEqual(invalid[0][0], 3)
        self.assertIn("715", invalid[0][1])

    def test_detects_multiple_invalid_rows(self):
        gff = (
            "A_MP\tflu-ngs\texon\t715\t452\t.\t+\t.\tID=exon3;Parent=A_M2\n"
            "A_MP\tflu-ngs\tCDS\t715\t452\t.\t+\t0\tID=cds3;Parent=A_M2"
        )
        invalid = wg.validate_gff_coordinates(gff)
        self.assertEqual(len(invalid), 2)

    def test_start_equals_end_is_valid(self):
        gff = "A_MP\tflu-ngs\texon\t500\t500\t.\t+\t.\tID=exon1"
        self.assertEqual([], wg.validate_gff_coordinates(gff))


class TestClassifyTransitionTransversion(unittest.TestCase):
    """
    Tests for classify_transition_transversion in merge-vep-irma.py.
    """

    def test_cant_handle_len2_nt_changes(self):
        """
        Check throws correct error if passed a mutation containing more than 1 nt.
        """
        msg = "Nucleotide must be one of ACGT to be classifed as transition or transversion."
        with self.assertRaisesRegex(ValueError, msg):
            mvi.classify_transition_transversion(
                {"Consensus_Allele": "AT", "Minority_Allele": "CT"}
            )

    def test_classifies_transition_correctly(self):
        """
        A simple test case.
        """
        out = mvi.classify_transition_transversion(
            {"Consensus_Allele": "G", "Minority_Allele": "A"}
        )
        self.assertEqual("transition", out)

    def test_classifies_transition_correctly(self):
        """
        A simple test case.
        """
        out = mvi.classify_transition_transversion(
            {"Consensus_Allele": "C", "Minority_Allele": "G"}
        )
        self.assertEqual("transversion", out)

    def test_returns_series_data(self):
        """
        Check that a series is returned, using all relevant tables in results.
        """
        for path in Path("../../results/irma").glob("*/tables/*-variants.tsv"):
            with self.subTest(path=path):
                df = pd.read_table(path).pipe(mvi.make_index_for_irma_variants)
                output = df.apply(mvi.classify_transition_transversion, axis=1)
                self.assertIsInstance(output, pd.Series)


class TestComputePadding(unittest.TestCase):
    """
    Tests for pad-incomplete-sequences.compute_padding.
    """

    def test_no_padding_needed(self):
        """Full-length sequence should need no padding."""
        ref = "ACGTACGTACGTACGT"
        query = "ACGTACGTACGTACGT"
        leading, trailing, *_ = pis.compute_padding(query, ref)
        self.assertEqual(0, leading)
        self.assertEqual(0, trailing)

    def test_leading_padding(self):
        """Sequence missing the start should get leading N's."""
        ref = "AACCTTGGAACCTTGG"
        query = "TTGGAACCTTGG"  # missing first 4 bases
        leading, trailing, *_ = pis.compute_padding(query, ref)
        self.assertEqual(4, leading)
        self.assertEqual(0, trailing)

    def test_trailing_padding(self):
        """Sequence missing the end should get trailing N's."""
        ref = "AACCTTGGCCAAGGTT"
        query = "AACCTTGGCCAA"  # missing last 4 bases
        leading, trailing, *_ = pis.compute_padding(query, ref)
        self.assertEqual(0, leading)
        self.assertEqual(4, trailing)

    def test_both_padding(self):
        """Sequence missing both ends should get both leading and trailing N's."""
        ref = "AACCTTGGCCAAGGTT"
        query = "CCTTGGCCAAGG"  # missing first 2 and last 2
        leading, trailing, *_ = pis.compute_padding(query, ref)
        self.assertEqual(2, leading)
        self.assertEqual(2, trailing)

    def test_no_internal_gaps(self):
        """Normal leading/trailing shortening should not flag internal gaps."""
        ref = "AACCTTGGCCAAGGTT"
        query = "CCTTGGCCAAGG"  # missing first 2 and last 2
        _, _, has_internal_gaps, _ = pis.compute_padding(query, ref)
        self.assertFalse(has_internal_gaps)

    def test_internal_gap_detected(self):
        """A query with an internal deletion should flag internal gaps."""
        ref = "AAAATTTTCCCCGGGG"
        query = "AAAACCCCGGGG"  # missing TTTT in the middle
        _, _, has_internal_gaps, _ = pis.compute_padding(query, ref)
        self.assertTrue(has_internal_gaps)


class TestPadIrmaDirInternalGaps(unittest.TestCase):
    """
    Tests that pad_irma_dir raises ValueError when internal gaps are present.
    """

    def test_raises_on_internal_gaps(self):
        """pad_irma_dir should raise ValueError for a segment with internal gaps."""
        ref_seq = "AAAATTTTCCCCGGGG"
        query_seq = "AAAACCCCGGGG"  # missing TTTT in the middle

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Write consensus FASTA with an internal deletion
            fasta_path = tmpdir / "A_HA.fasta"
            fasta_path.write_text(f">A_HA\n{query_seq}\n")

            # Build references dict matching the segment name
            references = {"A_HA": SeqRecord(Seq(ref_seq), id="A_HA")}

            with self.assertRaises(ValueError) as ctx:
                pis.pad_irma_dir(tmpdir, references, errors="warn")

            self.assertIn("A_HA", str(ctx.exception))
            self.assertIn("internal gaps", str(ctx.exception))


class TestPadFasta(unittest.TestCase):
    """
    Tests for pad-incomplete-sequences.pad_fasta.
    """

    def test_pads_fasta_file(self):
        """FASTA file should be padded with N's and header preserved."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">A_PB2 some description\nATCGATCG\n")
            path = Path(f.name)
        try:
            pis.pad_fasta(path, 3, 2)
            with open(path) as fobj:
                record = next(SeqIO.parse(fobj, "fasta"))
            self.assertEqual("NNN" + "ATCGATCG" + "NN", str(record.seq))
            self.assertEqual("A_PB2 some description", record.description)
        finally:
            path.unlink()

    def test_pads_fasta_file_correctly(self):
        """Verify exact padded sequence content."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">SEG\nACGT\n")
            path = Path(f.name)
        try:
            pis.pad_fasta(path, 2, 3)
            with open(path) as fobj:
                record = next(SeqIO.parse(fobj, "fasta"))
            self.assertEqual("NNACGTNNN", str(record.seq))
        finally:
            path.unlink()


class TestShiftVcf(unittest.TestCase):
    """
    Tests for pad-incomplete-sequences.shift_vcf.
    """

    def test_shifts_pos_column(self):
        """POS column should be shifted by the offset; header lines untouched."""
        vcf_content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "A_PB2\t10\t.\tA\tG\t.\t.\tDP=100\n"
            "A_PB2\t25\t.\tC\tT\t.\t.\tDP=200\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(vcf_content)
            path = Path(f.name)
        try:
            pis.shift_vcf(path, 100)
            lines = path.read_text().splitlines()
            # Header lines unchanged
            self.assertTrue(lines[0].startswith("##"))
            self.assertTrue(lines[1].startswith("#CHROM"))
            # Data lines shifted
            self.assertEqual("110", lines[2].split("\t")[1])
            self.assertEqual("125", lines[3].split("\t")[1])
        finally:
            path.unlink()


class TestShiftTable(unittest.TestCase):
    """
    Tests for pad-incomplete-sequences.shift_table.
    """

    def test_shifts_position_column(self):
        """Position column in variants table should be shifted."""
        content = (
            "Reference_Name\tPosition\tTotal\n"
            "A_PB2\t10\t500\n"
            "A_PB2\t20\t600\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            path = Path(f.name)
        try:
            pis.shift_table(path, "Position", 50)
            lines = path.read_text().splitlines()
            self.assertEqual("Position", lines[0].split("\t")[1])
            self.assertEqual("60", lines[1].split("\t")[1])
            self.assertEqual("70", lines[2].split("\t")[1])
        finally:
            path.unlink()

    def test_shifts_upstream_position_column(self):
        """Upstream_Position column in insertions/deletions should be shifted."""
        content = (
            "Reference_Name\tUpstream_Position\tInsert\n"
            "A_PB2\t5\tAA\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            path = Path(f.name)
        try:
            pis.shift_table(path, "Upstream_Position", 100)
            lines = path.read_text().splitlines()
            self.assertEqual("105", lines[1].split("\t")[1])
        finally:
            path.unlink()

    def test_missing_column_is_noop(self):
        """If the target column doesn't exist, do nothing."""
        content = "ColA\tColB\n1\t2\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            path = Path(f.name)
        try:
            pis.shift_table(path, "Position", 10)
            self.assertEqual(content, path.read_text())
        finally:
            path.unlink()


if __name__ == "__main__":
    unittest.main()
