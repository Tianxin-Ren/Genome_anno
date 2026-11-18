import os
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class BatchPathsTest(unittest.TestCase):
    def test_batch_reference_and_outputs(self):
        binning = (ROOT / "modules/binning/binning.smk").read_text()
        abundance = (ROOT / "modules/abundance/abundance.smk").read_text()
        for batch in ("Mpf_2026", "other", None):
            config = {"binning": {}}
            if batch:
                config["batch_name"] = batch
            ns = {"os": os, "config": config, "BN_RES": "/results/binning", "BN_LOGS": "/logs/binning",
                  "AB_RES": "/results/abundance", "AB_LOGS": "/logs/abundance",
                  "AB_C": {"coverm_use_drep": True}, "ab_full_path": lambda p: "/" + p}
            for source, names in ((binning, ("BN_BATCH_RES", "BN_BATCH_LOGS")),
                                  (abundance, ("AB_COVERM_RES", "AB_COVERM_LOGS", "AB_COVERM_MERGED", "AB_BIN_RES", "AB_BIN_BATCH_RES", "AB_COVERM_REF"))):
                for name in names:
                    exec(next(line for line in source.splitlines() if line.startswith(name + " = ")), ns)
            suffix = "/batches/" + batch if batch else ""
            self.assertEqual(ns["BN_BATCH_RES"], "/results/binning" + suffix)
            self.assertEqual(ns["AB_COVERM_REF"], os.path.join(ns["BN_BATCH_RES"], "drep/dereplicated_genomes"))
            self.assertEqual(ns["AB_COVERM_RES"], "/results/abundance" + suffix)
        self.assertIn('genomes = os.path.join(BN_BATCH_RES, "drep/dereplicated_genomes")', binning)
        self.assertIn('directory(os.path.join(BN_BATCH_RES, "gtdbtk"))', binning)


if __name__ == "__main__":
    unittest.main()
