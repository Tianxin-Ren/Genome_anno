#!/usr/bin/env python3
import unittest
import tempfile
import os
import shutil
import subprocess
import sys


class TestPartitionManifest(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_partition_manifest_")
        self.script = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "scripts", "partition_manifest.py")
        )

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _create_manifest(self, n_entries):
        path = os.path.join(self.test_dir, f"manifest_{n_entries}.tsv")
        with open(path, "w", encoding="utf-8") as f:
            for i in range(n_entries):
                f.write(f"mag_{i:04d}\t/path/to/mag_{i:04d}.fa\n")
        return path

    def test_basic_partition(self):
        manifest = self._create_manifest(25)
        out_dir = os.path.join(self.test_dir, "partitions")
        out_list = os.path.join(self.test_dir, "block_list.txt")

        cmd = [
            sys.executable,
            self.script,
            "--manifest", manifest,
            "--block-size", "10",
            "--out-dir", out_dir,
            "--out-list", out_list,
            "--no-shuffle",
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Error: {res.stderr}")

        # Should produce block_000, block_001, block_002
        with open(out_list) as f:
            blocks = [line.strip() for line in f if line.strip()]
        self.assertEqual(blocks, ["block_000", "block_001", "block_002"])

        # Check line counts in each block
        for b_id, expected_count in zip(blocks, [10, 10, 5]):
            b_path = os.path.join(out_dir, f"{b_id}.tsv")
            self.assertTrue(os.path.isfile(b_path))
            with open(b_path) as f:
                lines = [l.strip() for l in f if l.strip()]
            self.assertEqual(len(lines), expected_count)

    def test_exact_multiple(self):
        manifest = self._create_manifest(20)
        out_dir = os.path.join(self.test_dir, "partitions_exact")
        out_list = os.path.join(self.test_dir, "block_list_exact.txt")

        cmd = [
            sys.executable,
            self.script,
            "--manifest", manifest,
            "--block-size", "10",
            "--out-dir", out_dir,
            "--out-list", out_list,
            "--no-shuffle",
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Error: {res.stderr}")

        with open(out_list) as f:
            blocks = [line.strip() for line in f if line.strip()]
        self.assertEqual(blocks, ["block_000", "block_001"])

    def test_single_block_under_capacity(self):
        manifest = self._create_manifest(5)
        out_dir = os.path.join(self.test_dir, "partitions_single")
        out_list = os.path.join(self.test_dir, "block_list_single.txt")

        cmd = [
            sys.executable,
            self.script,
            "--manifest", manifest,
            "--block-size", "10000",
            "--out-dir", out_dir,
            "--out-list", out_list,
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Error: {res.stderr}")

        with open(out_list) as f:
            blocks = [line.strip() for line in f if line.strip()]
        self.assertEqual(blocks, ["block_000"])

        with open(os.path.join(out_dir, "block_000.tsv")) as f:
            lines = [l.strip() for l in f if l.strip()]
        self.assertEqual(len(lines), 5)

    def test_deterministic_shuffle(self):
        manifest = self._create_manifest(100)
        out_dir1 = os.path.join(self.test_dir, "p1")
        out_list1 = os.path.join(self.test_dir, "l1.txt")
        out_dir2 = os.path.join(self.test_dir, "p2")
        out_list2 = os.path.join(self.test_dir, "l2.txt")

        # Run 1 with seed 42
        subprocess.run(
            [sys.executable, self.script, "--manifest", manifest, "--block-size", "25",
             "--out-dir", out_dir1, "--out-list", out_list1, "--seed", "42"],
            check=True
        )
        # Run 2 with same seed 42
        subprocess.run(
            [sys.executable, self.script, "--manifest", manifest, "--block-size", "25",
             "--out-dir", out_dir2, "--out-list", out_list2, "--seed", "42"],
            check=True
        )

        with open(os.path.join(out_dir1, "block_000.tsv")) as f1, \
             open(os.path.join(out_dir2, "block_000.tsv")) as f2:
            self.assertEqual(f1.read(), f2.read())


if __name__ == "__main__":
    unittest.main()
