"""Nucleotide-only split-component tests against a reviewed development locus."""
import unittest

from Bio import SeqIO
from Bio.Seq import Seq

from .dna_rescue import ROOT, _fasta, _load_genes, load
from .split_component_scan import _complete_orf_window, scan
from .split_component_translation import audit


class SplitComponentTranslationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        archived = ROOT / "experiments/round201_sec_replay/runs/original"
        cls.manifest = load(archived / "dna_rescue/dna_rescue_manifest.json")
        cls.contigs = _fasta(cls.manifest["input"])
        cls.name, = cls.contigs
        cls.dna = cls.contigs[cls.name]
        cls.anchor = _load_genes(archived / "annotation")["CP026606.1_1109"]
        cls.anchor = {**cls.anchor, **{key: int(cls.anchor[key]) for key in ("start", "end", "strand")}}
        model_dir = ROOT / "references/processed/split_centers"
        cls.model = load(model_dir / "models.json")[1]
        cls.parameters = load(model_dir / "manifest.json")["parameters"]

    def search(self, dna, **kwargs):
        gene = {**self.anchor, **kwargs}
        result = scan([gene], {self.name: dna}, self.model, self.parameters, ROOT)
        return result.get(gene["id"], [])

    def test_real_jj_proposal_keeps_identity_and_membership_unconfirmed(self):
        evidence = audit(self.manifest, self.contigs)
        proposal, = evidence["hypotheses"]
        self.assertEqual(proposal["coding_start_status"], "multiple_compatible_starts")
        self.assertEqual(proposal["tga_genomic_start"], 1052760)
        self.assertEqual({(r["start"], r["end"]) for r in proposal["coding_alternatives"]},
                         {(1052703, 1052837), (1052694, 1052837)})
        self.assertEqual(proposal["translation_status"], "conditional_internal_TGA_no_U_assigned")
        self.assertEqual(proposal["catalytic_identity"], "unconfirmed")
        self.assertEqual(proposal["component_membership"], "unresolved")
        self.assertEqual(proposal["Sec_incorporation"], "unresolved")
        self.assertEqual(proposal["query_activity"], "not_tested")

    def test_real_rotated_called_fragment_preserves_one_conditional_locus(self):
        case = ROOT / "experiments/round203_vhuu_controls/runs/origin_rotated_fna"
        manifest = load(case / "dna_rescue/dna_rescue_manifest.json")
        evidence = audit(manifest, _fasta(manifest["input"]))
        proposal, = evidence["hypotheses"]
        self.assertEqual(proposal["called_component_fragment"]["terminal_TGA_start"],
                         proposal["tga_genomic_start"])
        length = len(self.dna)
        cut = 1052694
        loci = {((r["start"] + cut - 1) % length + 1, r["end"] - r["start"] + 1)
                for r in proposal["coding_alternatives"]}
        self.assertEqual(loci, {(1052694, 144), (1052703, 135)})

    def test_real_internal_taa_or_wrong_cys_abstains(self):
        for coordinate, codon in ((1052760, "TAA"), (1052769, "TCC")):
            with self.subTest(codon=codon):
                i = coordinate - 1
                mutated = self.dna[:i] + codon + self.dna[i + 3:]
                self.assertEqual(self.search(mutated), [])

    def test_strand_and_origin_transform_keep_exact_candidate_span(self):
        length = len(self.dna)
        complement = str(Seq(self.dna).reverse_complement())
        mirrored = dict(start=length + 1 - self.anchor["end"],
                        end=length + 1 - self.anchor["start"], strand=-1)
        loci = {(r["start"], r["end"]) for r in self.search(complement, **mirrored)}
        self.assertEqual(loci, {(length + 1 - 1052837, length + 1 - start) for start in (1052694, 1052703)})
        cut = 1052758
        start = (self.anchor["start"] - cut - 1) % length + 1
        rotated = self.dna[cut:] + self.dna[:cut]
        records = self.search(rotated, start=start, end=start + self.anchor["end"] - self.anchor["start"])
        self.assertEqual({((r["start"] + cut - 1) % length + 1, r["end"] - r["start"] + 1)
                          for r in records}, {(1052703, 135), (1052694, 144)})
        self.assertTrue(all(r["end"] > length for r in records))

    def test_two_local_candidates_remain_distinct(self):
        original = self.dna[1052702:1052837]
        offset = 1053000
        duplicate = self.dna[:offset] + original + self.dna[offset + len(original):]
        loci = {(r["start"], r["end"]) for r in self.search(duplicate)}
        self.assertTrue({(1052703, 1052837), (offset + 1, offset + 135)} <= loci)
        self.assertEqual({end for _, end in loci}, {1052837, offset + 135})


class CompleteOrfWindowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        record = SeqIO.read(ROOT / "experiments/round213_fn46_triage/inputs/NC_003551.1.gb", "genbank")
        cls.name, cls.dna = record.id, str(record.seq)
        features = {f.qualifiers["locus_tag"][0]: f for f in record.features
                    if f.type == "CDS" and f.qualifiers.get("locus_tag")}
        gene = features["MK_RS00945"]
        cls.anchor = dict(id="renamable_anchor", contig=cls.name, start=int(gene.location.start) + 1,
                          end=int(gene.location.end), strand=-1, topology="circular", contig_length=len(cls.dna))
        cls.small = features["MK_RS09355"]
        cls.nt = str(cls.small.extract(record.seq))
        model_dir = ROOT / "references/processed/split_centers"
        cls.model = load(model_dir / "models.json")[1]
        cls.parameters = load(model_dir / "manifest.json")["parameters"]

    def search(self, dna, **kwargs):
        anchor = dict(self.anchor, **kwargs)
        return scan([anchor], {anchor["contig"]: dna}, self.model, self.parameters, ROOT)[anchor["id"]]

    def edited(self, position, codon):
        offset = (position - 1) * 3
        nt = self.nt[:offset] + codon + self.nt[offset + 3:]
        return self.dna[:int(self.small.location.start)] + str(Seq(nt).reverse_complement()) + self.dna[int(self.small.location.end):]

    def test_native_complete_orf_recovers_only_a_conditional_translation(self):
        proposal, = self.search(self.dna)
        self.assertEqual((proposal["start"], proposal["end"], proposal["strand"]), (172723, 172851, -1))
        self.assertEqual(proposal["nt_sequence"], self.nt)
        self.assertEqual(proposal["window_alignment_method"], "complete_ORF_global_uncalibrated_v1")
        self.assertEqual(proposal["hsp_window_coverage"], 21 / 27)
        self.assertEqual(proposal["window_coverage"], 1)
        self.assertNotIn("U", proposal["translation_conditional"])
        self.assertEqual(proposal["translation_raw"].count("*"), 2)

    def test_true_termination_and_missing_ligand_do_not_extend(self):
        self.assertEqual(self.nt[60:63], "TGA")
        self.assertEqual(str(Seq(self.nt[69:72]).translate()), "C")
        # Position 7 lies before the unchanged HSP; its downstream start cannot cover the window.
        for position, codon in ((21, "TAA"), (24, "TCC"), (7, "TAA")):
            with self.subTest(position=position):
                self.assertEqual(self.search(self.edited(position, codon)), [])

    def test_same_dna_retains_its_span_after_reverse_rotation_and_rename(self):
        original, = self.search(self.dna)
        length = len(self.dna)
        reverse = str(Seq(self.dna).reverse_complement())
        reflected, = self.search(reverse, start=length + 1 - self.anchor["end"],
                                 end=length + 1 - self.anchor["start"], strand=1)
        self.assertEqual(reflected["DNA_sha256"], original["DNA_sha256"])
        cut = 172790
        start = (self.anchor["start"] - cut - 1) % length + 1
        rotated, = self.search(self.dna[cut:] + self.dna[:cut], start=start,
                              end=start + self.anchor["end"] - self.anchor["start"],
                              id="other_anchor_name", contig="other_contig_name")
        self.assertEqual(rotated["DNA_sha256"], original["DNA_sha256"])
        self.assertEqual((rotated["start"] + cut - 1) % length + 1, original["start"])
        self.assertEqual(rotated["end"] - rotated["start"], original["end"] - original["start"])

    def test_global_map_must_preserve_the_seed_sites_and_reject_ties(self):
        sites, offsets = {2: "U", 3: "C"}, {2: 1, 3: 2}
        nt = "ATGTGATGTGCTGCTTAA"
        raw, conditional = "M*CAA*", "MXCAA"
        result = _complete_orf_window("MXCAA", conditional, raw, nt, sites, offsets, 1, 5)
        self.assertEqual(result, (1, 80))
        self.assertIsNone(_complete_orf_window("MXCAA", conditional, raw, nt, sites, {2: 2, 3: 3}, 1, 5))
        self.assertIsNone(_complete_orf_window("MXCAAA", conditional, raw, nt, sites, offsets, 1, 6))
        self.assertIsNone(_complete_orf_window("MXCAA", conditional, raw, nt.replace("TGA", "TAA"), sites, offsets, 1, 5))


if __name__ == "__main__":
    unittest.main()
