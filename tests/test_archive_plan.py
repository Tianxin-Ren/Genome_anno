import csv
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from projects.config_read.scripts import archive_plan


class ArchivePlanTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source_root = self.root / "dssg" / "home"
        self.dest_root = self.root / "union" / "home"
        self.source_patch = patch.object(archive_plan, "SOURCE_ROOT", self.source_root)
        self.dest_patch = patch.object(archive_plan, "DEST_ROOT", self.dest_root)
        self.source_patch.start()
        self.dest_patch.start()
        self.addCleanup(self.source_patch.stop)
        self.addCleanup(self.dest_patch.stop)

        self.store = self.source_root / "acct" / "output_all" / "S1"
        self.qc = self.store / "00.QC"
        self.qc.mkdir(parents=True)
        for name in ("S1_R1_001.fastq.gz", "S1_R2_001.fastq.gz"):
            (self.store / name).write_bytes(name.encode())
        for mate in (1, 2):
            (self.qc / f"S1.clean.{mate}.fq.gz").write_bytes(f"clean{mate}".encode())
        self.intermediate = self.qc / "S1.merged.fq.gz"
        self.intermediate.write_bytes(b"intermediate")
        self.union_intermediate = self.dest_root / "acct" / "output_all" / "S1" / "00.QC" / "S1.fastp_p.json"
        self.union_intermediate.parent.mkdir(parents=True)
        self.union_intermediate.write_text("{}")
        self.marker = self.qc / "done"
        self.marker.write_text("done")
        self.results = self.source_root / "acct" / "results" / "MetaG-reads" / "sylph"
        (self.results / "S1").mkdir(parents=True)
        (self.results / "S1" / "S1.sylph.tsv").write_text("profile")
        (self.results / "S1" / "gtdb_r220_S1.sylphmpa").write_text("taxonomy")
        (self.results / "merged_sylph_profiling.tsv").write_text("merged")
        self.tsv = self.root / "samples.tsv"
        self.tsv.write_text(f"RUN_id\tStore_address\nS1\t{self.store}\n")
        self.plan = self.root / "archive_plan.tsv"
        with self.plan.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=archive_plan.FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(archive_plan.inventory(self.tsv, self.results))
        self.report = self.root / "archive_apply_report.tsv"

    def test_migrates_links_and_cleans_only_intermediates(self):
        archive_plan.apply_plan(self.plan, self.report)
        retained = (
            self.store / "S1_R1_001.fastq.gz",
            self.store / "S1_R2_001.fastq.gz",
            self.qc / "S1.clean.1.fq.gz",
            self.qc / "S1.clean.2.fq.gz",
            self.results / "S1" / "S1.sylph.tsv",
            self.results / "S1" / "gtdb_r220_S1.sylphmpa",
            self.results / "merged_sylph_profiling.tsv",
        )
        for path in retained:
            self.assertTrue(path.is_symlink(), path)
            self.assertTrue(path.is_file(), path)
            self.assertTrue(archive_plan.within(path.resolve(), archive_plan.actual_root(self.dest_root)))
        self.assertFalse(self.intermediate.exists())
        self.assertFalse(self.union_intermediate.exists())
        self.assertTrue(self.marker.is_file())
        archive_plan.apply_plan(self.plan, self.report)

    def test_conflicting_union_file_stops_before_any_change(self):
        dest = archive_plan.union_destination(Path(os.path.realpath(self.store / "S1_R1_001.fastq.gz")))
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"different")
        with self.assertRaisesRegex(RuntimeError, "Conflicting destination"):
            archive_plan.apply_plan(self.plan, self.report)
        self.assertFalse((self.store / "S1_R1_001.fastq.gz").is_symlink())
        self.assertTrue(self.intermediate.is_file())

    def test_changed_source_stops_before_any_change(self):
        (self.store / "S1_R1_001.fastq.gz").write_bytes(b"changed")
        with self.assertRaisesRegex(RuntimeError, "File changed since inventory"):
            archive_plan.apply_plan(self.plan, self.report)
        self.assertFalse((self.store / "S1_R2_001.fastq.gz").is_symlink())
        self.assertTrue(self.intermediate.is_file())

    def test_sample_directory_symlink_migrates_its_real_files(self):
        alias = self.source_root / "acct" / "linked_output_all" / "S1"
        alias.parent.mkdir(parents=True)
        alias.symlink_to(self.store, target_is_directory=True)
        self.tsv.write_text(f"RUN_id\tStore_address\nS1\t{alias}\n")
        with self.plan.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=archive_plan.FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(archive_plan.inventory(self.tsv, self.results))
        archive_plan.apply_plan(self.plan, self.report)
        self.assertTrue((self.store / "S1_R1_001.fastq.gz").is_symlink())
        self.assertTrue((alias / "S1_R1_001.fastq.gz").is_file())

    def test_missing_raw_mate_keeps_intermediates(self):
        (self.store / "S1_R2_001.fastq.gz").unlink()
        with self.plan.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=archive_plan.FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(archive_plan.inventory(self.tsv, self.results))
        archive_plan.apply_plan(self.plan, self.report)
        self.assertTrue(self.intermediate.is_file())


if __name__ == "__main__":
    unittest.main()
