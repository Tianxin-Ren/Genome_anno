import unittest

from Bio.Seq import Seq

from research_workflow.start_hypotheses import adjacent_starts, _neighbor_overlaps


def gene(dna, start, end, strand=1, topology="linear"):
    return dict(start=start, end=end, strand=strand, topology=topology, contig_length=len(dna),
                partial="False", coordinate_evidence="genomic_CDS")


def owned_gene(dna, gene_id, contig, start, end, strand=1, topology="linear", partial="False", evidence="genomic_CDS"):
    row = gene(dna, start, end, strand, topology)
    row.update(id=gene_id, contig=contig, partial=partial, coordinate_evidence=evidence)
    return row


class StartHypothesisTests(unittest.TestCase):
    def test_upstream_growth_reports_baseline_and_added_neighbor_overlap(self):
        called_coding = "ATGATGATGATGAAATAA"
        dna = "ATGATGATGATG" + called_coding
        selected = owned_gene(dna, "selected", "c", 13, 30)
        neighbor = owned_gene(dna, "neighbor", "c", 10, 20)
        rows = adjacent_starts(dna, selected, str(Seq(called_coding).translate(table=11, cds=True)), 11)
        alt = next(row for row in rows[1:] if row["start_offset_nt"] == -12)
        hits = _neighbor_overlaps({"c": [selected, neighbor]}, selected, rows[0], alt, len(dna))
        self.assertEqual(hits[0]["protein_id"], "neighbor")
        self.assertGreater(hits[0]["called_overlap_nt"], 0)
        self.assertGreater(hits[0]["alternative_overlap_nt"], hits[0]["called_overlap_nt"])
        self.assertGreater(hits[0]["added_overlap_nt"], 0)

    def test_negative_strand_growth_reports_opposite_neighbor(self):
        coding = "ATGATGATGATG" + "ATGATGAAATAA"
        dna = str(Seq(coding).reverse_complement())
        selected = owned_gene(dna, "selected", "c", 1, 12, strand=-1)
        neighbor = owned_gene(dna, "neighbor", "c", 13, 24, strand=1, partial="True")
        rows = adjacent_starts(dna, selected, str(Seq(coding[12:]).translate(table=11, cds=True)), 11)
        alt = next(row for row in rows[1:] if row["start_offset_nt"] == -12)
        hits = _neighbor_overlaps({"c": [neighbor, selected]}, selected, rows[0], alt, len(dna))
        self.assertEqual((hits[0]["strand_relation"], hits[0]["partial"]), ("opposite", True))

    def test_overlap_is_circular_and_sorted_independent_of_gene_order(self):
        called = "ATGATGATGATGAAATAA"; dna = "ATGATGATGATG" + called; dna = dna[1:] + dna[:1]
        selected = owned_gene(dna, "selected", "c", 12, 29, topology="circular")
        zed = owned_gene(dna, "z", "c", 1, 5)
        alpha = owned_gene(dna, "a", "c", 1, 3)
        rows = adjacent_starts(dna, selected, str(Seq(called).translate(table=11, cds=True)), 11)
        alt = next(row for row in rows[1:] if row["start_offset_nt"] == -12)
        hits = _neighbor_overlaps({"c": [zed, selected, alpha]}, selected, rows[0], alt, len(dna))
        self.assertEqual([hit["protein_id"] for hit in hits], ["a", "z"])
        self.assertTrue(all(hit["added_overlap_nt"] > 0 for hit in hits))

    def test_non_genomic_or_other_contig_and_removed_neighbor_are_ignored(self):
        dna = "ATGATGATGATG" + "ATGATGAAATAA"
        selected = owned_gene(dna, "selected", "c", 13, 24)
        alt = next(row for row in adjacent_starts(dna, selected, "MMK", 11)[1:] if row["start_offset_nt"] == -12)
        invalid = owned_gene(dna, "invalid", "c", 1, 12, evidence="predicted")
        invalid_strand = owned_gene(dna, "invalid_strand", "c", 1, 12, strand=0)
        other = owned_gene(dna, "other", "d", 1, 12)
        self.assertEqual(_neighbor_overlaps({"c": [selected, invalid, invalid_strand], "d": [other]}, selected,
                                            adjacent_starts(dna, selected, "MMK", 11)[0], alt, len(dna)), [])
        self.assertEqual(_neighbor_overlaps({"c": [selected]}, selected,
                                            adjacent_starts(dna, selected, "MMK", 11)[0], alt, len(dna)), [])

    def test_minus_strand_accepts_a_twelve_nt_extension(self):
        coding = "ATGATGATGATG" + "ATGATGAAATAA"
        dna = str(Seq(coding).reverse_complement())
        called_protein = str(Seq(coding[12:]).translate(table=11, cds=True))
        rows = adjacent_starts(dna, gene(dna, 1, 12, strand=-1), called_protein, 11)
        self.assertEqual(rows[0]["translation"], called_protein)
        self.assertEqual(next(r for r in rows[1:] if r["start_offset_nt"] == -12)["translation"], str(Seq(coding).translate(table=11, cds=True)))

    def test_table_valid_non_atg_is_caller_and_alternative(self):
        dna = "GTGGTGAAATAA"
        protein = str(Seq(dna).translate(table=11, cds=True))
        rows = adjacent_starts(dna, gene(dna, 1, len(dna)), protein, 11)
        self.assertEqual(rows[0]["start_codon"], "GTG")
        self.assertEqual(next(r for r in rows[1:] if r["start_offset_nt"] == 3)["start_codon"], "GTG")

    def test_intervening_stop_blocks_a_bounded_extension(self):
        dna = "ATGTAAGATGAAATAA"
        self.assertEqual(adjacent_starts(dna, gene(dna, 4, len(dna)), "MK", 11), [])

    def test_adjacent_start_keeps_the_stop_on_both_strands(self):
        dna = "ATGATGAAATAA"
        for seq, strand, expected in ((dna, 1, (4, 12)), (str(Seq(dna).reverse_complement()), -1, (1, 9))):
            with self.subTest(strand=strand):
                called, alternative = adjacent_starts(seq, gene(seq, 1, 12, strand), "MMK", 11)
                self.assertEqual((alternative["start"], alternative["end"]), expected)
                self.assertEqual(alternative["start_offset_nt"], 3)
                self.assertEqual(alternative["stop_codon"], called["stop_codon"])
                self.assertEqual(alternative["translation"], "MK")

    def test_upstream_atg_is_conditional_and_outside_linear_contig_is_excluded(self):
        dna = "AAAATGATGAAATAACC"
        called, alternative = adjacent_starts(dna, gene(dna, 7, 15), "MK", 11)
        self.assertEqual(called["translation"], "MK")
        self.assertEqual(alternative["translation"], "MMK")
        self.assertEqual(alternative["start_offset_nt"], -3)
        clipped = "ATGAAATAA"
        self.assertEqual(adjacent_starts(clipped, gene(clipped, 1, 9), "MK", 11), [])

    def test_wrap_and_reject_unproven_calls(self):
        dna = "ATGATGAAATAA"
        rotated = dna[5:] + dna[:5]
        row = gene(rotated, 8, 19, topology="circular")
        called, alternative = adjacent_starts(rotated, row, "MMK", 11)
        self.assertEqual(called["native_segments"], [(8, 12), (1, 7)])
        self.assertEqual(alternative["native_segments"], [(11, 12), (1, 7)])
        self.assertEqual(adjacent_starts(rotated, row, "MQK", 11), [])
        self.assertEqual(adjacent_starts(rotated, dict(row, partial="True"), "MMK", 11), [])
        self.assertEqual(adjacent_starts(rotated, row, "MMK", None), [])
        self.assertEqual(adjacent_starts(rotated, row, "MMK", 999), [])
        internal_stop = "ATGATGTAATAA"
        self.assertEqual(adjacent_starts(internal_stop, gene(internal_stop, 1, 12), "MM", 11), [])

    def test_circular_twelve_nt_extension_is_rotation_invariant(self):
        coding = "ATGATGATGATGAAATAA"
        contig = "ATGATGATGATG" + coding
        rotated = contig[1:] + contig[:1]
        row = gene(rotated, 12, 29, topology="circular")
        protein = str(Seq(coding).translate(table=11, cds=True))
        rows = adjacent_starts(rotated, row, protein, 11)
        self.assertEqual(rows[0]["translation"], protein)
        self.assertIn(-12, [r["start_offset_nt"] for r in rows[1:]])


if __name__ == "__main__":
    unittest.main()
