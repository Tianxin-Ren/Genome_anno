import os
import shlex
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import pandas as pd

from scripts import rescue_and_classify_gene_state as rescue
from scripts import extract_anchor_loci
from scripts import prepare_clinker_inputs
from scripts.select_and_validate_refs import (clustered_hits, select_colocalized_hits,
                                             reference_gene_coordinates, reference_locus_threshold)
from scripts.validate_syntenic_absence import find_deletion_junction, parse_paf_alignments, read_aligns_across_breakpoint
from scripts.infer_evolutionary_loss import load_independent_species_tree, run_fitch_parsimony


class TestSyntenyEvidence(unittest.TestCase):
    def test_clinker_inputs_crop_unlinked_tails_but_retain_anchor_interval(self):
        from Bio import SeqIO, SeqFeature
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        with tempfile.TemporaryDirectory() as directory:
            refs_dir = os.path.join(directory, "refs")
            out_dir = os.path.join(directory, "focused")
            os.makedirs(refs_dir)
            target_path = os.path.join(directory, "target_locus.gbk")
            reference_path = os.path.join(refs_dir, "GCA_1_locus.gbk")
            topology_path = os.path.join(directory, "topology.tsv")

            def make_record(record_id, genes):
                record = SeqRecord(Seq("A" * 100000), id=record_id, name=record_id, description="synthetic")
                record.annotations["molecule_type"] = "DNA"
                record.features = [SeqFeature.SeqFeature(SeqFeature.FeatureLocation(0, 100000), type="source")]
                for protein_id, start, end, label in genes:
                    record.features.append(SeqFeature.SeqFeature(
                        SeqFeature.FeatureLocation(start, end, strand=1), type="CDS",
                        qualifiers={"protein_id": [protein_id], "locus_tag": [protein_id], "gene": [label]}
                    ))
                return record

            target = make_record("target", [
                ("far_left", 1000, 1200, "hypothetical"),
                ("left_anchor", 30000, 30300, "K02111: atpA"),
                ("module_gene", 50000, 50400, "K02112: atpB"),
                ("right_anchor", 70000, 70300, "K02113: atpC"),
                ("far_right", 99000, 99200, "hypothetical"),
            ])
            reference = make_record("GCA_1", [
                ("far_left_ref", 1000, 1200, "hypothetical"),
                ("core_a", 49000, 49300, "K02111: atpA"),
                ("core_b", 50000, 50400, "K02112: atpB"),
                ("far_right_ref", 99000, 99200, "hypothetical"),
            ])
            with open(target_path, "w") as handle:
                SeqIO.write(target, handle, "genbank")
            with open(reference_path, "w") as handle:
                SeqIO.write(reference, handle, "genbank")
            with open(topology_path, "w") as handle:
                handle.write("left_matched_genes\tright_matched_genes\nleft_anchor\tright_anchor\n")
            with open(os.path.join(refs_dir, "ordered_clinker_files.txt"), "w") as handle:
                handle.write(reference_path + "\n")

            focused_target_path = prepare_clinker_inputs.prepare_inputs(
                target_path, topology_path, refs_dir, out_dir, "M00157",
                "K02111,K02112,K02113", 5000, flank_genes=0
            )
            focused_target = SeqIO.read(focused_target_path, "genbank")
            focused_reference = SeqIO.read(os.path.join(out_dir, "GCA_1_locus.focused.gbk"), "genbank")
            target_pids = {f.qualifiers.get("protein_id", [""])[0] for f in focused_target.features if f.type == "CDS"}
            ref_pids = {f.qualifiers.get("protein_id", [""])[0] for f in focused_reference.features if f.type == "CDS"}

            self.assertLess(len(focused_target.seq), len(target.seq))
            self.assertLess(len(focused_reference.seq), len(reference.seq))
            self.assertTrue({"left_anchor", "right_anchor", "module_gene"}.issubset(target_pids))
            self.assertNotIn("far_left", target_pids)
            self.assertNotIn("far_right", target_pids)
            self.assertTrue({"core_a", "core_b"}.issubset(ref_pids))
            self.assertNotIn("far_left_ref", ref_pids)
            self.assertNotIn("far_right_ref", ref_pids)
            with open(os.path.join(out_dir, "ordered_clinker_files.txt")) as handle:
                self.assertEqual(handle.read().strip(), os.path.join(out_dir, "GCA_1_locus.focused.gbk"))

    def test_reference_requires_distinct_colocalized_genes(self):
        coords = {"a": ("ctg", 100, 200, "+"), "b": ("ctg", 500, 600, "+"),
                  "c": ("other", 700, 800, "+")}
        self.assertTrue(clustered_hits({"q1": "a", "q2": "b"}, coords, 2))
        self.assertFalse(clustered_hits({"q1": "a", "q2": "c"}, coords, 2))
        self.assertFalse(clustered_hits({"q1": "a"}, coords, 2))
        self.assertEqual(select_colocalized_hits({"q1": "a", "q2": "b", "q3": "c"}, coords, 2),
                         {"q1": "a", "q2": "b"})
        self.assertTrue(clustered_hits({"q1": "a", "q2": "b"},
                                       {**coords, "interrupt": ("ctg", 300, 400, "+")}, 2))
        self.assertTrue(clustered_hits({"q1": "a", "q2": "b"},
                                       {**coords, "b": ("ctg", 500, 600, "-")}, 2))
        self.assertEqual(reference_locus_threshold(9), 4)
        self.assertEqual(reference_locus_threshold(2), 2)

    def test_stage1_reference_sort_uses_dict_coordinates(self):
        records = [types.SimpleNamespace(id="p2"), types.SimpleNamespace(id="p1")]
        coords = {"p1": {"contig": "ctg", "start": 100, "end": 200, "strand": "+"},
                  "p2": {"contig": "ctg", "start": 10, "end": 90, "strand": "+"}}
        grouped = extract_anchor_loci.group_reference_proteins_by_coordinate(records, coords)
        self.assertEqual([record.id for _, record in grouped["ctg"]], ["p2", "p1"])

    def test_reference_coordinates_require_observed_positions(self):
        with tempfile.TemporaryDirectory() as directory:
            faa = os.path.join(directory, "ref.faa")
            with open(faa, "w") as handle:
                handle.write(">ctg_1 # 101 # 201 # 1 # gene\nMKK\n>unknown\nMKK\n")
            coords = reference_gene_coordinates(faa)
            self.assertEqual(coords["ctg_1"], ("ctg", 101, 201, "+"))
            self.assertNotIn("unknown", coords)

    def test_operon_structure_uses_order_strand_spacing_and_interruption(self):
        genes = [{"query_gene": "q1", "gene_state": "PRESENT", "hit_id": "a"},
                 {"query_gene": "q2", "gene_state": "PRESENT", "hit_id": "b"}]
        coords = {"a": ("ctg", 100, 200, "+"), "b": ("ctg", 250, 350, "+")}
        result = rescue.classify_operon_structure(genes, coords, ["q1", "q2"], ["a", "b"])
        self.assertEqual(result["operon_prediction"], "CANDIDATE_OPERON")
        self.assertEqual(result["max_intergenic_bp"], 50)
        self.assertEqual(rescue.classify_operon_structure(genes, coords, ["q1", "q2"], ["a", "b", "other"])["operon_structure_state"], "INTERRUPTED_CLUSTER")
        coords["b"] = ("ctg", 250, 350, "-")
        self.assertEqual(rescue.classify_operon_structure(genes, coords, ["q1", "q2"], ["a", "b"])["operon_structure_state"], "DISRUPTED_ORIENTATION")

    def test_deletion_requires_two_collinear_alignment_blocks(self):
        left = {"qname": "ref", "tname": "target", "strand": "+", "qstart": 0,
                "qend": 1000, "tstart": 0, "tend": 1000, "aln_len": 1000, "mapq": 60}
        right = {**left, "qstart": 1000, "qend": 2000, "tstart": 4000, "tend": 5000}
        self.assertEqual(find_deletion_junction([left, right])[0], 3000)
        self.assertIsNone(find_deletion_junction([left]))
        self.assertIsNone(find_deletion_junction([left, {**right, "tname": "other"}]))

    def test_single_alignment_cigar_deletion_supports_sample_loss(self):
        block = {"qname": "sample", "tname": "positive_ref", "strand": "+",
                 "qstart": 0, "qend": 2000, "tstart": 0, "tend": 5000,
                 "aln_len": 5000, "mapq": 60, "cigar": "1000M3000D1000M"}
        self.assertEqual(find_deletion_junction([block]), (3000, 0, 1000))
        self.assertIsNone(find_deletion_junction([block], locus_start=1200, locus_end=1500))
        self.assertIsNone(find_deletion_junction([{**block, "cigar": "1000M3000I1000M"}]))
        with tempfile.TemporaryDirectory() as directory:
            paf = os.path.join(directory, "alignment.paf")
            with open(paf, "w") as handle:
                handle.write("sample\t2000\t0\t2000\t+\tpositive_ref\t5000\t0\t5000\t1900\t5000\t60\tcg:Z:1000M3000D1000M\n")
            self.assertEqual(find_deletion_junction(parse_paf_alignments(paf)), (3000, 0, 1000))

    def test_junction_read_must_align_on_both_sides_without_deletion(self):
        self.assertTrue(read_aligns_across_breakpoint(951, "100M", 1000))
        self.assertFalse(read_aligns_across_breakpoint(941, "40M100D40M", 1000))
        self.assertFalse(read_aligns_across_breakpoint(941, "40M100N40M", 1000))

    def test_no_tree_does_not_fabricate_topology(self):
        with tempfile.TemporaryDirectory() as directory:
            output = os.path.join(directory, "tree.nwk")
            root, count, support = load_independent_species_tree("", "target", [], [], output)
            self.assertIsNone(root)
            self.assertEqual(count, 0)
            self.assertEqual(support, 0)
            with open(output) as handle:
                self.assertIn("No independently inferred", handle.read())

    def test_rooted_independent_tree_supports_nested_loss(self):
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "species.nwk")
            output = os.path.join(directory, "used.nwk")
            with open(source, "w") as handle:
                handle.write("(((target:1,Pos1:1)95:1,Pos2:1)95:1,Outgroup:1);\n")
            positive = [{"accession": "Pos1", "taxonomy": "g__G1"},
                        {"accession": "Pos2", "taxonomy": "g__G2"}]
            negative = [{"accession": "Outgroup", "taxonomy": "g__Out"}]
            root, lineages, support = load_independent_species_tree(
                source, "target", positive, negative, output, "Outgroup")
            self.assertEqual(lineages, 2)
            self.assertEqual(support, 95.0)
            self.assertEqual(run_fitch_parsimony(root)[0], "1")

    def test_dna_rescue_uses_dna_query_and_protein_database(self):
        with tempfile.TemporaryDirectory() as directory:
            query = os.path.join(directory, "query.faa")
            target_faa = os.path.join(directory, "target.faa")
            locus = os.path.join(directory, "locus.fa")
            topology = os.path.join(directory, "topology.tsv")
            gene_out = os.path.join(directory, "genes.tsv")
            cluster_out = os.path.join(directory, "cluster.tsv")
            with open(query, "w") as handle:
                handle.write(">q1\nMKKKKKKKKK\n")
            open(target_faa, "w").close()
            with open(locus, "w") as handle:
                handle.write(">target_locus\n" + "ATG" * 100 + "\n")
            pd.DataFrame([{"target_slice_start": 0, "ancestral_interval_start": 0,
                           "ancestral_interval_end": 300, "target_contig": "ctg",
                           "same_contig": "YES", "intervening_genes": ""}]).to_csv(topology, sep="\t", index=False)
            commands = []
            def fake_command(command):
                commands.append(command)
                tokens = shlex.split(command)
                if "blastx" in tokens:
                    out = tokens[tokens.index("-o") + 1]
                    with open(out, "w") as handle:
                        handle.write("target_locus\tq1\t90\t10\t0\t0\t1\t30\t1\t10\t1e-20\t70\t300\t10\tMKKKKKKKKK\tMKKKKKKKKK\n")
                return types.SimpleNamespace(returncode=0, stdout="", stderr="")
            args = ["rescue", "--target_name", "T", "--target_locus_fa", locus,
                    "--target_faa", target_faa, "--query_operon_faa", query,
                    "--topology_tsv", topology, "--out_gene_tsv", gene_out,
                    "--out_cluster_tsv", cluster_out]
            with patch.object(sys, "argv", args), patch.object(rescue, "run_cmd", side_effect=fake_command):
                rescue.main()
            self.assertIn(f"blastx -q {locus}", "\n".join(commands))
            self.assertEqual(pd.read_csv(gene_out, sep="\t").iloc[0]["gene_state"], "ANNOTATION_MISS")

    def test_old_diamond_field_falls_back_without_overclaiming_intact_orf(self):
        with tempfile.TemporaryDirectory() as directory:
            query = os.path.join(directory, "query.faa")
            target_faa = os.path.join(directory, "target.faa")
            locus = os.path.join(directory, "locus.fa")
            topology = os.path.join(directory, "topology.tsv")
            gene_out = os.path.join(directory, "genes.tsv")
            cluster_out = os.path.join(directory, "cluster.tsv")
            with open(query, "w") as handle:
                handle.write(">q1\nMKKKKKKKKK\n")
            open(target_faa, "w").close()
            with open(locus, "w") as handle:
                handle.write(">target_locus\n" + "ATG" * 100 + "\n")
            pd.DataFrame([{"target_slice_start": 0, "ancestral_interval_start": 0,
                           "ancestral_interval_end": 300, "target_contig": "ctg",
                           "same_contig": "YES", "intervening_genes": ""}]).to_csv(topology, sep="\t", index=False)
            commands = []
            def fake_command(command):
                commands.append(command)
                if "diamond blastx" in command and "qseq_translated" in command:
                    raise RuntimeError("Error: Invalid output field: qseq_translated")
                if "diamond blastx" in command:
                    tokens = shlex.split(command)
                    with open(tokens[tokens.index("-o") + 1], "w") as handle:
                        handle.write("target_locus\tq1\t90\t10\t0\t0\t1\t30\t1\t10\t1e-20\t70\t300\t10\tMKKKKKKKKK\n")
                return types.SimpleNamespace(returncode=0, stdout="", stderr="")
            args = ["rescue", "--target_name", "T", "--target_locus_fa", locus,
                    "--target_faa", target_faa, "--query_operon_faa", query,
                    "--topology_tsv", topology, "--out_gene_tsv", gene_out,
                    "--out_cluster_tsv", cluster_out]
            with patch.object(sys, "argv", args), patch.object(rescue, "run_cmd", side_effect=fake_command):
                rescue.main()
            gene = pd.read_csv(gene_out, sep="\t").iloc[0]
            self.assertEqual(gene["gene_state"], "AMBIGUOUS")
            self.assertEqual(gene["mutation_evidence_status"], "VERSION_LIMITED")
            self.assertEqual(sum("qseq_translated" not in cmd for cmd in commands if "diamond blastx" in cmd), 1)
            cluster = pd.read_csv(cluster_out, sep="\t").iloc[0]
            self.assertEqual(cluster["mutation_evidence_status"], "VERSION_LIMITED")
            self.assertEqual(cluster["A_local"], 0)
            self.assertEqual(cluster["A_global"], 0)

    def test_unknown_protein_location_is_global_presence(self):
        with tempfile.TemporaryDirectory() as directory:
            query = os.path.join(directory, "query.faa")
            target_faa = os.path.join(directory, "target.faa")
            locus = os.path.join(directory, "locus.fa")
            topology = os.path.join(directory, "topology.tsv")
            gene_out = os.path.join(directory, "genes.tsv")
            cluster_out = os.path.join(directory, "cluster.tsv")
            with open(query, "w") as handle:
                handle.write(">q1\nMKKKKKKKKK\n")
            with open(target_faa, "w") as handle:
                handle.write(">unknown_id\nMKKKKKKKKK\n")
            open(locus, "w").close()
            pd.DataFrame([{"target_slice_start": 0, "ancestral_interval_start": 0,
                           "ancestral_interval_end": 300, "target_contig": "ctg",
                           "same_contig": "YES", "intervening_genes": ""}]).to_csv(topology, sep="\t", index=False)
            def fake_command(command):
                tokens = shlex.split(command)
                if "blastp" in tokens:
                    with open(tokens[tokens.index("-o") + 1], "w") as handle:
                        handle.write("q1\tunknown_id\t100\t10\t0\t0\t1\t10\t1\t10\t1e-20\t70\t10\t10\n")
                return types.SimpleNamespace(returncode=0, stdout="", stderr="")
            args = ["rescue", "--target_name", "T", "--target_locus_fa", locus,
                    "--target_faa", target_faa, "--query_operon_faa", query,
                    "--topology_tsv", topology, "--out_gene_tsv", gene_out,
                    "--out_cluster_tsv", cluster_out]
            with patch.object(sys, "argv", args), patch.object(rescue, "run_cmd", side_effect=fake_command):
                rescue.main()
            result = pd.read_csv(cluster_out, sep="\t").iloc[0]
            self.assertEqual(result["protein_detected_count"], 1)
            self.assertEqual(result["T_global"], 1.0)
            self.assertEqual(result["A_global"], 0)

    def test_no_positive_reference_emits_unresolved_stage1_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            refs = os.path.join(directory, "refs.tsv")
            anchors = os.path.join(directory, "anchors.tsv")
            topology = os.path.join(directory, "stage1", "topology.tsv")
            locus = os.path.join(directory, "stage1", "locus.fa")
            gbk = os.path.join(directory, "stage1", "locus.gbk")
            extended = os.path.join(directory, "stage1", "extended.tsv")
            ref_dir = os.path.join(directory, "stage1", "refs")
            pd.DataFrame([{"validation_status": "UNRESOLVED_NO_POSITIVE_REFERENCE"}]).to_csv(refs, sep="\t", index=False)
            pd.DataFrame(columns=["anchor_id"]).to_csv(anchors, sep="\t", index=False)
            args = ["stage1", "--target_name", "T", "--target_fa", "missing.fa",
                    "--target_faa", "missing.faa", "--anchors_tsv", anchors,
                    "--refs_tsv", refs, "--out_topology_tsv", topology,
                    "--out_target_locus_fa", locus, "--out_target_locus_gbk", gbk,
                    "--out_refs_locus_dir", ref_dir, "--out_extended_proteins_tsv", extended]
            with patch.object(sys, "argv", args), patch.object(extract_anchor_loci, "load_kegg_ko_definitions", return_value={}):
                extract_anchor_loci.main()
            self.assertEqual(pd.read_csv(topology, sep="\t").iloc[0]["localization_status"], "NO_VALIDATED_POSITIVE_REF")
            self.assertEqual(os.path.getsize(locus), 0)
            self.assertTrue(os.path.isfile(extended))


if __name__ == "__main__":
    unittest.main()
