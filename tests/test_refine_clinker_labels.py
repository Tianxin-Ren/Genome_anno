#!/usr/bin/env python3
import unittest
import json
import tempfile
import os
import re
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from scripts.refine_clinker_labels import refine_clinker_html, extract_contig_from_genes, clean_id_spaces

class TestRefineClinker(unittest.TestCase):
    def test_extract_contig(self):
        # Case 1: orig_coord in names
        genes1 = [{"names": {"orig_coord": "CAZOPV010000010.1:100..200(+)"}}]
        self.assertEqual(extract_contig_from_genes(genes1), "CAZOPV010000010.1")

        # Case 2: orig_coord inside note
        genes2 = [{"names": {"note": "KO: K15888 | orig_coord=ctg_99:50..500(-) | score=100"}}]
        self.assertEqual(extract_contig_from_genes(genes2), "ctg_99")

        # Case 3: triple underscore in protein_id
        genes3 = [{"names": {"protein_id": "GCA_964582085.1___contig_alpha___prot00018"}}]
        self.assertEqual(extract_contig_from_genes(genes3), "contig_alpha")

    def test_clean_id_spaces(self):
        # Spaces in long ID should be cleaned
        self.assertEqual(clean_id_spaces("CA ZOPV010000010.1"), "CAZOPV010000010.1")
        # English spaces in functional annotations should be preserved
        self.assertEqual(clean_id_spaces("K15888: hypothetical protein"), "K15888: hypothetical protein")
        self.assertEqual(clean_id_spaces("membrane protein"), "membrane protein")

    def test_refine_clinker_html(self):
        mock_data = {
            "clusters": [
                {
                    "name": "GCA_964582085_1",
                    "loci": [
                        {
                            "name": "locus_1",
                            "start": 1000,
                            "end": 5000,
                            "genes": [
                                {
                                    "names": {
                                        "locus_tag": "prot00018",
                                        "gene": "K15888: putative membrane protein (sc: 329.5/202.9, e: 9.2e-97)",
                                        "product": "putative membrane protein",
                                        "note": "KO: K15888 | orig_coord=CAZOPV010000010.1:1000..5000(+)"
                                    },
                                    "label": "prot00018"
                                }
                            ]
                        }
                    ]
                },
                {
                    "name": "GCA_963760905.1_locus",
                    "loci": [
                        {
                            "name": "locus_2",
                            "start": 500,
                            "end": 3000,
                            "genes": [
                                {
                                    "names": {
                                        "locus_tag": "CAZOPV010000010.1_18",
                                        "gene": "K15888: ppaX; inorganic pyrophosphatase (sc: 215.0/180.0, e: 1.0e-50)",
                                        "product": "inorganic pyrophosphatase",
                                        "orig_coord": "CAZOPV010000010.1:500..3000(+)",
                                        "note": "KO: K15888"
                                    },
                                    "label": "CAZOPV010000010.1_18"
                                }
                            ]
                        }
                    ]
                }
            ]
        }
        
        raw_html = f"""<!DOCTYPE html><html><body><script>
const data = {json.dumps(mock_data)};
</script></body></html>"""

        with tempfile.NamedTemporaryFile(suffix=".html", mode="w+", delete=False) as tf:
            tf.write(raw_html)
            tf_path = tf.name

        try:
            refine_clinker_html(tf_path)
            with open(tf_path, "r", encoding="utf-8") as f:
                content = f.read()

            match = re.search(r"const data\s*=\s*(\{.*?\});", content, re.DOTALL)
            self.assertTrue(match is not None)
            res_data = json.loads(match.group(1))

            # Verify Target Cluster
            c0 = res_data["clusters"][0]
            self.assertEqual(c0["name"], "GCA_964582085.1 (CAZOPV010000010.1)")
            self.assertEqual(c0["loci"][0]["name"], "CAZOPV010000010.1:1000-5000")
            g0 = c0["loci"][0]["genes"][0]
            self.assertEqual(g0["label"], "K15888: putative membrane protein (sc: 329.5/202.9, e: 9.2e-97)")

            # Verify Ref Cluster
            c1 = res_data["clusters"][1]
            self.assertEqual(c1["name"], "GCA_963760905.1 (CAZOPV010000010.1)")
            self.assertEqual(c1["loci"][0]["name"], "CAZOPV010000010.1:500-3000")
            g1 = c1["loci"][0]["genes"][0]
            self.assertEqual(g1["label"], "K15888: ppaX; inorganic pyrophosphatase (sc: 215.0/180.0, e: 1.0e-50)")

        finally:
            if os.path.exists(tf_path):
                os.remove(tf_path)

if __name__ == "__main__":
    unittest.main()
