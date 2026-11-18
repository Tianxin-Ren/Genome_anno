#!/usr/bin/env python3
import unittest
import tempfile
import os
import shutil
import csv
import subprocess
import sys


class TestChainDrepClusters(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_chain_drep_")
        self.script = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "scripts", "chain_drep_clusters.py")
        )

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_multi_block_chaining(self):
        # Setup Block 0
        b0_dir = os.path.join(self.test_dir, "round1_blocks", "block_000", "data_tables")
        os.makedirs(b0_dir, exist_ok=True)
        b0_cdb = os.path.join(b0_dir, "Cdb.csv")
        b0_widb = os.path.join(b0_dir, "Widb.csv")

        with open(b0_cdb, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["genome", "secondary_cluster", "primary_cluster"])
            w.writerow(["bin.1.fa", "1_0", "1"])
            w.writerow(["bin.2.fa", "1_0", "1"])
            w.writerow(["bin.3.fa", "2_0", "2"])

        with open(b0_widb, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["genome", "score", "cluster"])
            w.writerow(["bin.1.fa", "95.0", "1_0"])
            w.writerow(["bin.3.fa", "88.0", "2_0"])

        # Setup Block 1
        b1_dir = os.path.join(self.test_dir, "round1_blocks", "block_001", "data_tables")
        os.makedirs(b1_dir, exist_ok=True)
        b1_cdb = os.path.join(b1_dir, "Cdb.csv")
        b1_widb = os.path.join(b1_dir, "Widb.csv")

        with open(b1_cdb, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["genome", "secondary_cluster", "primary_cluster"])
            w.writerow(["bin.4.fa", "1_0", "1"])
            w.writerow(["bin.5.fa", "1_0", "1"])
            w.writerow(["bin.6.fa", "2_0", "2"])

        with open(b1_widb, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["genome", "score", "cluster"])
            w.writerow(["bin.4.fa", "94.0", "1_0"])
            w.writerow(["bin.6.fa", "85.0", "2_0"])

        # Setup Pool dRep (pools bin.1, bin.3 from b0, and bin.4, bin.6 from b1)
        pool_dir = os.path.join(self.test_dir, "pool", "data_tables")
        os.makedirs(pool_dir, exist_ok=True)
        pool_cdb = os.path.join(pool_dir, "Cdb.csv")
        pool_widb = os.path.join(pool_dir, "Widb.csv")

        # In pool: bin.1.fa and bin.4.fa merge into final cluster 'P1_0', winner is bin.1.fa
        # bin.3.fa is 'P2_0', bin.6.fa is 'P3_0'
        with open(pool_cdb, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["genome", "secondary_cluster", "primary_cluster"])
            w.writerow(["bin.1.fa", "P1_0", "P1"])
            w.writerow(["bin.4.fa", "P1_0", "P1"])
            w.writerow(["bin.3.fa", "P2_0", "P2"])
            w.writerow(["bin.6.fa", "P3_0", "P3"])

        with open(pool_widb, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["genome", "score", "cluster"])
            w.writerow(["bin.1.fa", "95.0", "P1_0"])
            w.writerow(["bin.3.fa", "88.0", "P2_0"])
            w.writerow(["bin.6.fa", "85.0", "P3_0"])

        out_cdb = os.path.join(self.test_dir, "round1", "data_tables", "Cdb.csv")
        out_widb = os.path.join(self.test_dir, "round1", "data_tables", "Widb.csv")
        out_prov = os.path.join(self.test_dir, "round1", "data_tables", "provenance.tsv")

        cmd = [
            sys.executable,
            self.script,
            "--block-cdbs", b0_cdb, b1_cdb,
            "--block-widbs", b0_widb, b1_widb,
            "--pool-cdb", pool_cdb,
            "--pool-widb", pool_widb,
            "--out-cdb", out_cdb,
            "--out-widb", out_widb,
            "--out-provenance", out_prov,
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Error: {res.stderr}")

        # Check master Cdb
        self.assertTrue(os.path.isfile(out_cdb))
        with open(out_cdb, "r") as f:
            reader = list(csv.DictReader(f))

        cdb_map = {r["genome"]: r["secondary_cluster"] for r in reader}
        self.assertEqual(len(cdb_map), 6)
        # All of bin.1, bin.2, bin.4, bin.5 should now belong to P1_0!
        self.assertEqual(cdb_map["bin.1.fa"], "P1_0")
        self.assertEqual(cdb_map["bin.2.fa"], "P1_0")
        self.assertEqual(cdb_map["bin.4.fa"], "P1_0")
        self.assertEqual(cdb_map["bin.5.fa"], "P1_0")
        self.assertEqual(cdb_map["bin.3.fa"], "P2_0")
        self.assertEqual(cdb_map["bin.6.fa"], "P3_0")

        # Check provenance table
        self.assertTrue(os.path.isfile(out_prov))
        with open(out_prov, "r") as f:
            prov_reader = list(csv.DictReader(f, delimiter="\t"))
        prov_map = {r["genome"]: r for r in prov_reader}

        # Check bin.5.fa provenance: from block_001, winner was bin.4.fa, final winner is bin.1.fa
        b5 = prov_map["bin.5.fa"]
        self.assertEqual(b5["block_id"], "block_001")
        self.assertEqual(b5["block_winner"], "bin.4.fa")
        self.assertEqual(b5["final_secondary_cluster"], "P1_0")
        self.assertEqual(b5["final_winner"], "bin.1.fa")

    def test_single_block_passthrough(self):
        b0_dir = os.path.join(self.test_dir, "round1_blocks", "block_000", "data_tables")
        os.makedirs(b0_dir, exist_ok=True)
        b0_cdb = os.path.join(b0_dir, "Cdb.csv")
        b0_widb = os.path.join(b0_dir, "Widb.csv")

        with open(b0_cdb, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["genome", "secondary_cluster", "primary_cluster"])
            w.writerow(["bin.1.fa", "1_0", "1"])

        with open(b0_widb, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["genome", "score", "cluster"])
            w.writerow(["bin.1.fa", "95.0", "1_0"])

        out_cdb = os.path.join(self.test_dir, "round1", "Cdb.csv")
        out_widb = os.path.join(self.test_dir, "round1", "Widb.csv")

        cmd = [
            sys.executable,
            self.script,
            "--block-cdbs", b0_cdb,
            "--block-widbs", b0_widb,
            "--out-cdb", out_cdb,
            "--out-widb", out_widb,
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Error: {res.stderr}")
        self.assertTrue(os.path.isfile(out_cdb))
        self.assertTrue(os.path.isfile(out_widb))


if __name__ == "__main__":
    unittest.main()
