from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from map_gene_counts import load_gene_kos, process_summary, write_output
from scripts.deduplicate_proteins import build_catalog
from scripts.expand_dedup_hmm import expand


class HmmPostprocessTests(unittest.TestCase):
    def test_exact_dedup_and_hmm_expansion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_faa = root / "all.faa"
            unique_faa = root / "unique.faa"
            mapping = root / "map.tsv"
            stats = root / "stats.tsv"
            database = root / "catalog.sqlite"
            raw_tbl = root / "raw.tbl"
            expanded_tbl = root / "expanded.tbl"
            input_faa.write_text(">g1|p1\nMKK*\n>g2|p2\nMKK\n>g3|p3\nMNN\n")
            build_catalog(str(input_faa), str(unique_faa), str(mapping), str(stats), str(database))
            self.assertEqual(unique_faa.read_text().count(">"), 2)
            raw_tbl.write_text("# hmmsearch\nHMM\t-\tg1|p1\t-\t1e-20\t50\t0\n")
            expand(str(raw_tbl), str(expanded_tbl), str(database), str(mapping))
            self.assertIn("g1|p1", expanded_tbl.read_text())
            self.assertIn("g2|p2", expanded_tbl.read_text())

    def test_gene_counts_retain_genomes_without_target_hits(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gene_map = root / "genes.tsv"
            summary = root / "summary.tsv"
            output = root / "nested" / "gene_counts.tsv"

            gene_map.write_text(
                "KO\tdefinition\tgene_name\tec_number\n"
                "K00330\tdescription\tnuoA; subunit A\t\n"
            )
            summary.write_text(
                "Name\tProtein_ID\tScore_Annotation\tEvalue_Annotation\n"
                "genome_hit\tp1\tK00330\tK00330\n"
                "genome_zero\tp2\tK99999\tK99999\n"
            )

            genes, ko_to_genes = load_gene_kos(str(gene_map))
            counts = process_summary(str(summary), ko_to_genes)
            write_output(str(output), genes, counts)

            rows = output.read_text().splitlines()
            self.assertEqual(rows[0], "Genome_Name\tnuoA; subunit A")
            self.assertIn("genome_hit\t1", rows)
            self.assertIn("genome_zero\t0", rows)


if __name__ == "__main__":
    unittest.main()
