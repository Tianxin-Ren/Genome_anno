#!/usr/bin/env python3
import unittest
import tempfile
import os
import shutil
import gzip
import csv
import subprocess
import sys
import time


class TestPrepareDrepInputs(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_drep_prep_")
        self.script = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "scripts", "prepare_drep_inputs.py")
        )

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_basic_staging_and_genome_info(self):
        # Create dummy fastas
        fa1 = os.path.join(self.test_dir, "bin.1.fa")
        with open(fa1, "w") as f:
            f.write(">c1\nATCG\n")

        fa2_gz = os.path.join(self.test_dir, "bin.2.fasta.gz")
        with gzip.open(fa2_gz, "wb") as f:
            f.write(b">c2\nGGCC\n")

        manifest = os.path.join(self.test_dir, "manifest.tsv")
        with open(manifest, "w") as f:
            f.write(f"bin.1\t{fa1}\n")
            f.write(f"bin.2\t{fa2_gz}\n")

        ginfo_src = os.path.join(self.test_dir, "raw_info.csv")
        with open(ginfo_src, "w") as f:
            f.write("Genome,Completeness,Contamination\n")
            f.write("bin.1.fa,98.5,1.2\n")
            f.write("bin.2,85.0,3.5\n")

        stage_dir = os.path.join(self.test_dir, "staged")
        out_list = os.path.join(self.test_dir, "genomes_list.txt")
        ginfo_dst = os.path.join(self.test_dir, "aligned_info.csv")

        cmd = [
            sys.executable,
            self.script,
            "--manifest", manifest,
            "--stage-dir", stage_dir,
            "--out-list", out_list,
            "--genome-info-src", ginfo_src,
            "--genome-info-dst", ginfo_dst,
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Script failed: {res.stderr}")

        # Check staged files
        self.assertTrue(os.path.isfile(os.path.join(stage_dir, "bin.1.fa")))
        self.assertTrue(os.path.isfile(os.path.join(stage_dir, "bin.2.fa")))

        # Check list
        with open(out_list) as f:
            lines = [line.strip() for line in f if line.strip()]
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[0], os.path.abspath(os.path.join(stage_dir, "bin.1.fa")))
        self.assertEqual(lines[1], os.path.abspath(os.path.join(stage_dir, "bin.2.fa")))

        # Check output genomeInfo.csv
        with open(ginfo_dst) as f:
            reader = list(csv.reader(f))
        self.assertEqual(reader[0], ["genome", "completeness", "contamination"])
        self.assertEqual(reader[1], ["bin.1.fa", "98.5", "1.2"])
        self.assertEqual(reader[2], ["bin.2.fa", "85.0", "3.5"])

    def test_missing_genome_info_fallback(self):
        # Staged MAG not in raw_info.csv
        fa1 = os.path.join(self.test_dir, "mag_x.fa")
        with open(fa1, "w") as f:
            f.write(">c1\nATCG\n")

        manifest = os.path.join(self.test_dir, "manifest.tsv")
        with open(manifest, "w") as f:
            f.write(f"mag_x\t{fa1}\n")

        ginfo_src = os.path.join(self.test_dir, "empty_info.csv")
        with open(ginfo_src, "w") as f:
            f.write("genome,completeness,contamination\n")

        stage_dir = os.path.join(self.test_dir, "staged")
        out_list = os.path.join(self.test_dir, "genomes_list.txt")
        ginfo_dst = os.path.join(self.test_dir, "aligned_info.csv")

        cmd = [
            sys.executable,
            self.script,
            "--manifest", manifest,
            "--stage-dir", stage_dir,
            "--out-list", out_list,
            "--genome-info-src", ginfo_src,
            "--genome-info-dst", ginfo_dst,
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("WARNING: 1/1 staged genomes were missing", res.stderr)

        with open(ginfo_dst) as f:
            reader = list(csv.reader(f))
        self.assertEqual(reader[1], ["mag_x.fa", "0", "100"])

    def test_no_genome_info_provided(self):
        fa1 = os.path.join(self.test_dir, "mag_y.fa")
        with open(fa1, "w") as f:
            f.write(">c1\nATCG\n")

        manifest = os.path.join(self.test_dir, "manifest.tsv")
        with open(manifest, "w") as f:
            f.write(f"mag_y\t{fa1}\n")

        stage_dir = os.path.join(self.test_dir, "staged")
        out_list = os.path.join(self.test_dir, "genomes_list.txt")
        ginfo_dst = os.path.join(self.test_dir, "aligned_info.csv")

        cmd = [
            sys.executable,
            self.script,
            "--manifest", manifest,
            "--stage-dir", stage_dir,
            "--out-list", out_list,
            "--genome-info-dst", ginfo_dst,
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertTrue(os.path.exists(ginfo_dst))
        self.assertEqual(os.path.getsize(ginfo_dst), 0)

    def test_large_scale_performance(self):
        fa = os.path.join(self.test_dir, "dummy.fa")
        with open(fa, "w") as f:
            f.write(">c1\nATCG\n")

        manifest = os.path.join(self.test_dir, "large_manifest.tsv")
        with open(manifest, "w") as f:
            for i in range(10000):
                f.write(f"mag_{i}\t{fa}\n")

        ginfo_src = os.path.join(self.test_dir, "large_info.csv")
        with open(ginfo_src, "w") as f:
            f.write("genome,completeness,contamination\n")
            for i in range(10000):
                f.write(f"mag_{i}.fa,90.0,2.0\n")

        stage_dir = os.path.join(self.test_dir, "large_staged")
        out_list = os.path.join(self.test_dir, "large_genomes_list.txt")
        ginfo_dst = os.path.join(self.test_dir, "large_aligned_info.csv")

        start = time.time()
        cmd = [
            sys.executable,
            self.script,
            "--manifest", manifest,
            "--stage-dir", stage_dir,
            "--out-list", out_list,
            "--genome-info-src", ginfo_src,
            "--genome-info-dst", ginfo_dst,
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = time.time() - start
        self.assertEqual(res.returncode, 0)
        self.assertLess(elapsed, 15.0, f"10,000 items staging took {elapsed:.2f}s, expected < 15s")

        with open(out_list) as f:
            lines = [l for l in f if l.strip()]
        self.assertEqual(len(lines), 10000)


if __name__ == "__main__":
    unittest.main()
