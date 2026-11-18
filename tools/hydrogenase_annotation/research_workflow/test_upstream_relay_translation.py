"""Regression tests for bounded two-site relay translation."""
import unittest

from Bio.Seq import Seq

from .dna_rescue import ROOT, _fasta, _load_genes, load
from .split_component_translation import _upstream_relay, audit


class UpstreamRelayTranslationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        archived = ROOT / "experiments/round201_sec_replay/runs/original"
        cls.manifest = load(archived / "dna_rescue/dna_rescue_manifest.json")
        cls.contigs = _fasta(cls.manifest["input"])
        cls.contig, = cls.contigs
        cls.dna = cls.contigs[cls.contig]
        cls.genes = _load_genes(archived / "annotation")
        cls.anchor = cls.genes["CP026606.1_1109"]
        cls.anchor = {**cls.anchor, "strand": 1, "protein_sha256": "", "mapped_sites": []}

    def search(self, dna, *anchors):
        return _upstream_relay(list(anchors or (self.anchor,)), {self.contig: dna}, {}, ROOT)

    def test_bound_anchor_recovers_one_locus_without_making_claims(self):
        evidence = audit(self.manifest, self.contigs)
        relay, = evidence["upstream_relay_translation"]["hypotheses"]
        self.assertEqual(relay["tga_genomic_starts"], [1050058, 1050211])
        self.assertEqual(relay["terminal_stop_genomic_start"], 1050421)
        self.assertEqual(relay["component_membership"], "unresolved")
        self.assertEqual(relay["catalytic_identity"], "unconfirmed")
        self.assertEqual(relay["Sec_incorporation"], "unresolved")
        self.assertEqual(relay["query_activity"], "not_tested")
        self.assertEqual([(r["start"], r["end"]) for r in relay["coding_alternatives"]],
                         [(1049992, 1050423), (1050016, 1050423), (1050019, 1050423)])
        for row in relay["coding_alternatives"]:
            self.assertEqual(row["translation_conditional"].count("X"), 2)
            self.assertNotIn("U", row["translation_conditional"])
            self.assertEqual(row["upstream_gap_bp"], 1007)
        self.assertEqual(len(evidence["hypotheses"]), 1)

    def test_either_internal_taa_or_lost_relay_cys_withdraws(self):
        for start, codon in [(1050058, "TAA"), (1050211, "TAA"), (1050049, "TCC")]:
            with self.subTest(start=start):
                mutated = self.dna[:start - 1] + codon + self.dna[start + 2:]
                self.assertEqual(self.search(mutated)["hypotheses"], [])

    def test_internal_rotation_and_reverse_strand_keep_exact_nt(self):
        length = len(self.dna)
        start, end = (int(self.anchor[k]) for k in ("start", "end"))
        cut = 1050100
        rotated = self.dna[cut:] + self.dna[:cut]
        anchored = dict(self.anchor, start=start - cut, end=end - cut)
        shifted, = self.search(rotated, anchored)["hypotheses"]
        self.assertEqual([(i + cut - 1) % length + 1 for i in shifted["tga_genomic_starts"]],
                         [1050058, 1050211])
        self.assertTrue(all(r["end"] > length for r in shifted["coding_alternatives"]))
        reverse = str(Seq(self.dna).reverse_complement())
        mirrored = dict(self.anchor, start=length + 1 - end, end=length + 1 - start, strand=-1)
        reflected, = self.search(reverse, mirrored)["hypotheses"]
        self.assertEqual([length - i - 1 for i in reflected["tga_genomic_starts"]], [1050058, 1050211])
        self.assertEqual({r["DNA_sha256"] for r in shifted["coding_alternatives"]},
                         {r["DNA_sha256"] for r in reflected["coding_alternatives"]})

    def test_competing_loci_or_anchors_withhold(self):
        original = self.dna[1050018:1050423]
        offset = 1049400
        duplicate = self.dna[:offset] + original + self.dna[offset + len(original):]
        local = self.search(duplicate)
        self.assertEqual(local["hypotheses"], [])
        self.assertEqual(local["searches"][0]["locus_count"], 2)
        other = dict(self.anchor, id="competing_anchor", start=1051600, end=1052856)
        shared = self.search(self.dna, self.anchor, other)
        self.assertEqual(shared["hypotheses"], [])
        self.assertEqual({r["status"] for r in shared["searches"]},
                         {"ambiguous_conditional_locus_or_anchor"})


if __name__ == "__main__":
    unittest.main()
