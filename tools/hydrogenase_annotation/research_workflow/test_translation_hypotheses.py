import unittest

from Bio.Seq import Seq

from translation_hypotheses import circular_coding_extent, coding_extent, ligand_uga


class TranslationHypothesesTests(unittest.TestCase):
    def proposal(self, changes, frame=1):
        sites = [17, 20, 45, 48]
        codons = ["AAA"] * 70
        for p in sites:
            codons[p - 1] = "TGT"
        for p, codon in changes.items():
            codons[p - 1] = codon
        sequence = Seq("".join(codons))
        aligned = str(sequence.translate(table=4))
        dna = str(sequence if frame > 0 else sequence.reverse_complement())
        return ligand_uga(dict(qstart=1, qseq="A" * 70, sseq=aligned,
                               dna_start=1, dna_end=len(dna), frame=frame,
                               translation_table="4", evalue=1e-30), dna, sites)

    def test_dual_uga_preserves_raw_translation_and_reverse_coordinates(self):
        for frame, expected in ((1, [49, 58]), (-1, [160, 151])):
            with self.subTest(frame=frame):
                hit = self.proposal({17: "TGA", 20: "TGA"}, frame)
                self.assertEqual(hit["ligand_pattern"], "UXXU/CXXC")
                self.assertEqual(hit["uga_codon_starts_1based"], expected)
                self.assertEqual(hit["raw_translation_table_11"][16:20], "*KK*")
                self.assertEqual(hit["conditional_sec_translation"][16:20], "UKKU")

    def test_other_stop_or_nonligand_uga_cannot_be_replaced(self):
        for changes in ({17: "TAA"}, {18: "TGA"}, {17: "TGA", 35: "TGA"},
                        {17: "TGA", 20: "GCT"}, {}):
            with self.subTest(changes=changes):
                self.assertIsNone(self.proposal(changes))

    def test_verified_ordinary_cds_is_preserved_on_both_strands(self):
        coding = "GTGAAATGAAAATGAAACTAA"
        protein = {"ordinary": "MK"}
        for strand in (1, -1):
            with self.subTest(strand=strand):
                dna = coding if strand == 1 else str(Seq(coding).reverse_complement())
                proposal = (dict(start=1, end=18, frame=1, uga_codon_starts_1based=[7, 13]) if strand == 1
                            else dict(start=4, end=21, frame=-1, uga_codon_starts_1based=[13, 7]))
                gene = dict(id="ordinary", start=1 if strand == 1 else 13, end=9 if strand == 1 else 21,
                            coordinate_evidence="genomic_CDS", partial="False")
                result = coding_extent(proposal, dna, [gene], protein)
                self.assertEqual(result["status"], "conditional_start_to_stop")
                self.assertEqual((result["start"], result["end"]), (1, 21))
                self.assertEqual(result["conditional_sec_translation"], "MKUKUN")
                self.assertEqual(result["raw_translation_table_11"], "VK*K*N*")
                self.assertEqual(result["ordinary_cds"]["protein_length"], 2)
                self.assertEqual(coding_extent(proposal, dna, [gene], {"ordinary": "VK"})["status"],
                                 "no_verified_ordinary_cds_anchor")

    def test_unobserved_or_ambiguous_downstream_sequence_is_not_completed(self):
        anchor = dict(id="ordinary", start=1, end=9, coordinate_evidence="genomic_CDS", partial="False")
        proposal = dict(start=1, end=15, frame=1, uga_codon_starts_1based=[7])
        prefix = "ATGAAATGAAAATGC"
        for tail, status in (("", "contig_edge_before_terminal_stop"),
                             ("NNNTAA", "ambiguous_sequence_before_terminal_stop"),
                             ("AAA" * 100 + "TAA", "terminal_stop_outside_search_limit")):
            with self.subTest(status=status):
                result = coding_extent(proposal, prefix+tail, [anchor], {"ordinary": "MK"})
                self.assertEqual(result["status"], status)
                self.assertNotIn("conditional_sec_translation", result)

    def test_circular_hypothesis_restores_original_start_and_stop(self):
        coding = "GTGAAATGAAAATGAAACTAA"
        native = coding + "GCG" * 4
        rotated = native[10:] + native[:10]
        proposal = dict(start=24, end=41, frame=1, uga_codon_starts_1based=[30, 36])
        gene = dict(id="ordinary", start=24, end=32, coordinate_evidence="genomic_CDS",
                    partial="False", topology="circular", contig_length=len(rotated))
        result = circular_coding_extent(proposal, rotated, [gene], {"ordinary": "MK"})
        self.assertEqual(result["status"], "conditional_start_to_stop")
        self.assertEqual(result["native_segments"], [(24, 33), (1, 11)])
        self.assertEqual(result["conditional_sec_translation"], "MKUKUN")
        self.assertEqual(result["ordinary_cds"]["native_segments"], [(24, 32)])
