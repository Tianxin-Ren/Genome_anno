#!/usr/bin/env python3
import unittest
import tempfile
import os
import shutil
import subprocess
import sys


class TestPoolMagGenesAndProteins(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_pool_")
        self.script = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "scripts", "pool_mag_genes_and_proteins.py")
        )

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_pool_batches(self):
        batches_dir = os.path.join(self.test_dir, "batches")
        b1 = os.path.join(batches_dir, "batch_0001")
        b2 = os.path.join(batches_dir, "batch_0002")
        for b in [b1, b2]:
            os.makedirs(os.path.join(b, "faa"), exist_ok=True)
            os.makedirs(os.path.join(b, "ffn"), exist_ok=True)
            os.makedirs(os.path.join(b, "gff"), exist_ok=True)

        with open(os.path.join(b1, "faa", "mag1.faa"), "w") as f:
            f.write(">mag1_p1\nMK\n")
        with open(os.path.join(b1, "ffn", "mag1.ffn"), "w") as f:
            f.write(">mag1_g1\nATG\n")
        with open(os.path.join(b1, "gff", "mag1.gff"), "w") as f:
            f.write("##gff-version 3\n")

        with open(os.path.join(b2, "faa", "mag2.faa"), "w") as f:
            f.write(">mag2_p1\nVL\n")
        with open(os.path.join(b2, "ffn", "mag2.ffn"), "w") as f:
            f.write(">mag2_g1\nGTG\n")
        with open(os.path.join(b2, "gff", "mag2.gff"), "w") as f:
            f.write("##gff-version 3\n")

        batch_list = os.path.join(batches_dir, "batch_list.txt")
        with open(batch_list, "w") as f:
            f.write("batch_0001\nbatch_0002\n")

        proteins_dir = os.path.join(self.test_dir, "proteins", "all")
        genes_dir = os.path.join(self.test_dir, "genes")
        done_file = os.path.join(self.test_dir, "proteins", ".pool.done")

        cmd = [
            sys.executable,
            self.script,
            "--batches_dir", batches_dir,
            "--batch_list", batch_list,
            "--proteins_dir", proteins_dir,
            "--genes_dir", genes_dir,
            "--done_file", done_file,
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Script failed: {res.stderr}")

        # Check proteins
        self.assertTrue(os.path.islink(os.path.join(proteins_dir, "mag1.faa")))
        self.assertTrue(os.path.islink(os.path.join(proteins_dir, "mag2.faa")))

        # Check genes
        self.assertTrue(os.path.islink(os.path.join(genes_dir, "mag1.ffn")))
        self.assertTrue(os.path.islink(os.path.join(genes_dir, "mag2.ffn")))
        self.assertTrue(os.path.islink(os.path.join(genes_dir, "mag1.gff")))
        self.assertTrue(os.path.islink(os.path.join(genes_dir, "mag2.gff")))

        # Check done file
        self.assertTrue(os.path.isfile(done_file))

        # Test idempotency (re-run should succeed without errors)
        res2 = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res2.returncode, 0)


if __name__ == "__main__":
    unittest.main()
