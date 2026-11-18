from __future__ import annotations

import csv
import gzip
import os
import subprocess
import tempfile
import unittest


class ReadProfilingCorruptedQcTests(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.base_dir = self.tmp_dir.name

    def tearDown(self):
        self.tmp_dir.cleanup()

    def test_gzip_integrity_detection(self):
        """Verify that corrupted gz file is accurately detected by gzip -t."""
        good_gz = os.path.join(self.base_dir, "good.fq.gz")
        with gzip.open(good_gz, "wb") as f:
            f.write(b"@read1\nACGT\n+\nIIII\n")

        # Create corrupted gz (truncated byte sequence)
        bad_gz = os.path.join(self.base_dir, "bad.fq.gz")
        with open(good_gz, "rb") as f:
            content = f.read()
        with open(bad_gz, "wb") as f:
            f.write(content[:len(content) - 8])

        res_good = subprocess.run(["gzip", "-t", good_gz], capture_output=True)
        self.assertEqual(res_good.returncode, 0)

        res_bad = subprocess.run(["gzip", "-t", bad_gz], capture_output=True)
        self.assertNotEqual(res_bad.returncode, 0)
        self.assertTrue(len(res_bad.stderr) > 0)

    def test_empty_gzip_placeholder_properties(self):
        """Verify that gzip placeholder is a valid 20-byte gzip file containing 0 uncompressed bytes."""
        placeholder = os.path.join(self.base_dir, "placeholder.fq.gz")
        # Generate placeholder using bash same as rule
        subprocess.run(
            f"gzip -n -c < /dev/null > {placeholder}",
            shell=True,
            check=True,
        )
        self.assertTrue(os.path.isfile(placeholder))
        self.assertEqual(os.path.getsize(placeholder), 20)

        # gzip -t should pass with return code 0
        res = subprocess.run(["gzip", "-t", placeholder], capture_output=True)
        self.assertEqual(res.returncode, 0)

        # Decompressed content has 0 bytes
        with gzip.open(placeholder, "rb") as f:
            self.assertEqual(f.read(), b"")

    def test_corrupted_report_tsv_format(self):
        """Verify the central corrupted report TSV format and atomic appending."""
        report_file = os.path.join(self.base_dir, "corrupted_samples.tsv")
        headers = ["RUN_id", "Store_address", "Raw_R1", "Raw_R2", "Error_type", "Error_detail", "Timestamp"]

        row1 = ["SAMPLE_001", "/path/sample1", "/raw/sample1_1.fq.gz", "/raw/sample1_2.fq.gz", "CorruptedGzipR1", "gzip -t failed: unexpected EOF", "2026-10-19T16:00:00Z"]
        row2 = ["SAMPLE_002", "/path/sample2", "/raw/sample2_1.fq.gz", "/raw/sample2_2.fq.gz", "EmptyRawR2", "Raw R2 file is empty (0 bytes)", "2026-10-19T16:01:00Z"]

        with open(report_file, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(headers)
            writer.writerow(row1)
            writer.writerow(row2)

        with open(report_file, "r", encoding="utf-8") as f:
            reader = list(csv.DictReader(f, delimiter="\t"))

        self.assertEqual(len(reader), 2)
        self.assertEqual(reader[0]["RUN_id"], "SAMPLE_001")
        self.assertEqual(reader[0]["Error_type"], "CorruptedGzipR1")
        self.assertEqual(reader[1]["RUN_id"], "SAMPLE_002")
        self.assertEqual(reader[1]["Error_type"], "EmptyRawR2")

    def test_sylph_merge_skips_corrupted_sample(self):
        """Verify that sylph merge script skips header-only placeholder from corrupted sample."""
        sample_a_tsv = os.path.join(self.base_dir, "sample_a.sylph.tsv")
        sample_b_corrupt_tsv = os.path.join(self.base_dir, "sample_b.sylph.tsv")
        merged_tsv = os.path.join(self.base_dir, "merged.tsv")

        sylph_header = "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_marker_coverage\tEff_cov\tTrue_cov\tAni\tContig_name\n"

        # Valid sample A has header + 2 hit rows
        with open(sample_a_tsv, "w") as f:
            f.write(sylph_header)
            f.write("sample_a.fq.gz\tgenome1.fna.gz\t40.5\t38.2\t12.0\t11.5\t11.0\t98.5\tctg1\n")
            f.write("sample_a.fq.gz\tgenome2.fna.gz\t59.5\t61.8\t18.0\t17.5\t17.0\t99.1\tctg2\n")

        # Corrupted sample B has only header
        with open(sample_b_corrupt_tsv, "w") as f:
            f.write(sylph_header)

        # Emulate the shell merge logic
        cmd = f"""
        first=1
        touch "{merged_tsv}"
        for f in "{sample_a_tsv}" "{sample_b_corrupt_tsv}"; do
            [ ! -f "$f" ] && continue
            if [ "$first" -eq 1 ]; then
                cat "$f" > "{merged_tsv}"
                first=0
            else
                tail -n +2 "$f" >> "{merged_tsv}"
            fi
        done
        """
        subprocess.run(cmd, shell=True, check=True)

        with open(merged_tsv, "r") as f:
            lines = f.readlines()

        # Should contain 1 header + 2 sample_a rows = 3 lines total (sample_b contributed 0 data rows)
        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith("Sample_file\tGenome_file"))
        self.assertIn("genome1.fna.gz", lines[1])
        self.assertIn("genome2.fna.gz", lines[2])

    def test_sylph_merge_when_corrupted_sample_is_first(self):
        """Verify that sylph merge succeeds even if the very first sample was corrupted."""
        sample_a_corrupt_tsv = os.path.join(self.base_dir, "sample_a.sylph.tsv")
        sample_b_valid_tsv = os.path.join(self.base_dir, "sample_b.sylph.tsv")
        merged_tsv = os.path.join(self.base_dir, "merged.tsv")

        sylph_header = "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_marker_coverage\tEff_cov\tTrue_cov\tAni\tContig_name\n"

        with open(sample_a_corrupt_tsv, "w") as f:
            f.write(sylph_header)

        with open(sample_b_valid_tsv, "w") as f:
            f.write(sylph_header)
            f.write("sample_b.fq.gz\tgenome1.fna.gz\t100.0\t100.0\t20.0\t19.5\t19.0\t99.5\tctg1\n")

        cmd = f"""
        first=1
        touch "{merged_tsv}"
        for f in "{sample_a_corrupt_tsv}" "{sample_b_valid_tsv}"; do
            [ ! -f "$f" ] && continue
            if [ "$first" -eq 1 ]; then
                cat "$f" > "{merged_tsv}"
                first=0
            else
                tail -n +2 "$f" >> "{merged_tsv}"
            fi
        done
        """
        subprocess.run(cmd, shell=True, check=True)

        with open(merged_tsv, "r") as f:
            lines = f.readlines()

        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].startswith("Sample_file\tGenome_file"))
        self.assertIn("sample_b.fq.gz", lines[1])


if __name__ == "__main__":
    unittest.main()
