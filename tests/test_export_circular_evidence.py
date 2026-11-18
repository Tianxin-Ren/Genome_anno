import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None

ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable


@unittest.skipIf(pysam is None, "pysam is required for BAM integration tests")
class TestCircularEvidenceExport(unittest.TestCase):
    def test_long_read_small_indel_is_not_a_structural_gap(self):
        from scripts.export_circular_evidence import span_junction
        read = pysam.AlignedSegment()
        read.reference_start = 0
        read.cigarstring = "100M2D100M"
        self.assertTrue(span_junction(read, 101, 100, max_indel=50))
        self.assertFalse(span_junction(read, 101, 100))
        read.cigarstring = "100M60D100M"
        self.assertFalse(span_junction(read, 130, 100, max_indel=50))
        read.cigarstring = "100M2N100M"
        self.assertFalse(span_junction(read, 101, 100, max_indel=50))

    def write_bam(self, path, ref_name, ref_len, records):
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": ref_name, "LN": ref_len}]}
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            for spec in records:
                rec = pysam.AlignedSegment()
                rec.query_name = spec["qname"]
                rec.flag = spec.get("flag", 0)
                rec.reference_id = 0
                rec.reference_start = spec.get("start", 0)
                rec.mapping_quality = spec.get("mapq", 60)
                rec.cigarstring = spec["cigar"]
                rec.query_sequence = "A" * rec.infer_query_length()
                rec.query_qualities = pysam.qualitystring_to_array("I" * len(rec.query_sequence))
                rec.next_reference_id = spec.get("mate_ref", -1)
                rec.next_reference_start = spec.get("mate_start", -1)
                rec.template_length = spec.get("tlen", 0)
                bam.write(rec)
        pysam.index(str(path))

    def make_inputs(self, td):
        td = Path(td)
        contig = td / "circular.fa"
        contig.write_text(">ctg1\n" + "A" * 101 + "\n")
        coords = td / "rotated_junction_coords.tsv"
        coords.write_text("contig\tseq_id\tlength\trotation_offset\trotated_junction_coord\tis_circular\nctg1\tctg1\t101\t5\t96\tYES\n")
        long_bam = td / "long.bam"
        self.write_bam(long_bam, "ctg1", 101, [
            {"qname": "primary", "start": 0, "mapq": 60, "cigar": "10M"},
            {"qname": "secondary", "start": 10, "mapq": 60, "cigar": "5M", "flag": 256},
            {"qname": "lowmapq", "start": 15, "mapq": 0, "cigar": "5M"},
            {"qname": "clip_left", "start": 30, "mapq": 0, "cigar": "2S5M"},
            {"qname": "clip_right", "start": 90, "mapq": 0, "cigar": "5M2S"},
        ])
        long_junction = td / "long_junction.bam"
        self.write_bam(long_junction, "ctg1_junction", 100, [
            {"qname": "span_insert", "start": 45, "mapq": 60, "cigar": "5M1I5M"},
            {"qname": "span_deletion", "start": 45, "mapq": 60, "cigar": "5M2D5M"},
            {"qname": "span_skip", "start": 45, "mapq": 60, "cigar": "5M2N5M"},
        ])
        short_bam = td / "short.bam"
        self.write_bam(short_bam, "ctg1", 101, [{"qname": "short", "start": 0, "mapq": 60, "cigar": "8M"}])
        rotated = td / "rotated.bam"
        self.write_bam(rotated, "ctg1", 101, [
            {"qname": "pair_ok", "flag": 99, "start": 40, "mapq": 60, "cigar": "56M", "mate_ref": 0, "mate_start": 96, "tlen": 61},
            {"qname": "pair_bad_mate", "flag": 99, "start": 40, "mapq": 60, "cigar": "56M", "mate_ref": 0, "mate_start": 96, "tlen": 61},
            {"qname": "pair_ok", "flag": 147, "start": 96, "mapq": 60, "cigar": "5M", "mate_ref": 0, "mate_start": 40, "tlen": -61},
            {"qname": "pair_bad_mate", "flag": 147, "start": 96, "mapq": 0, "cigar": "5M", "mate_ref": 0, "mate_start": 40, "tlen": -61},
        ])
        return contig, coords, long_bam, long_junction, short_bam, rotated

    def run_export(self, td, setup=True):
        out = Path(td) / "evidence"
        if setup:
            contig, coords, long_bam, long_junction, short_bam, rotated = self.make_inputs(td)
        else:
            contig = Path(td) / "circular.fa"
            coords = Path(td) / "rotated_junction_coords.tsv"
            long_bam, long_junction = Path(td) / "long.bam", Path(td) / "long_junction.bam"
            short_bam, rotated = Path(td) / "short.bam", Path(td) / "rotated.bam"
        cmd = [PYTHON, str(ROOT / "scripts/export_circular_evidence.py"), "--contig", "ctg1", "--contig-fa", str(contig), "--long-bam", str(long_bam), "--long-junction-bam", str(long_junction), "--short-bam", str(short_bam), "--short-rotated-bam", str(rotated), "--rotation-coords", str(coords), "--out-dir", str(out), "--bin-size", "10", "--junction-window", "5", "--min-mapq", "20", "--long-anchor", "2", "--short-anchor", "2"]
        result = subprocess.run(cmd, env=os.environ.copy(), text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return out

    @staticmethod
    def rows(path):
        with gzip.open(path, "rt", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_real_measurements_and_paper_style_policies(self):
        with tempfile.TemporaryDirectory() as td:
            out = self.run_export(td)
            coverage = self.rows(out / "genome_coverage.tsv.gz")
            self.assertEqual({r["track"] for r in coverage}, {"long", "short"})
            first = next(r for r in coverage if r["track"] == "long" and r["start_0based"] == "10")
            self.assertGreater(float(first["all_mean_depth"]), float(first["mean_depth"]))
            reads = self.rows(out / "junction_reads.tsv.gz")
            self.assertEqual(next(r for r in reads if r["qname"] == "span_insert")["spans_junction"].lower(), "true")
            self.assertEqual(next(r for r in reads if r["qname"] == "span_deletion")["spans_junction"].lower(), "false")
            self.assertEqual(next(r for r in reads if r["qname"] == "span_skip")["spans_junction"].lower(), "false")
            long_jdepth = next(r for r in self.rows(out / "junction_depth.tsv.gz") if r["track"] == "long" and r["relative_pos"] == "0")
            self.assertEqual(int(long_jdepth["reference_pos_0based"]), 50)
            self.assertEqual(int(long_jdepth["original_pos_0based"]), 0)
            self.assertEqual(int(long_jdepth["depth"]), 1)
            short_jdepth = next(r for r in self.rows(out / "junction_depth.tsv.gz") if r["track"] == "short" and r["relative_pos"] == "0")
            self.assertEqual(int(short_jdepth["reference_pos_0based"]), 96)
            self.assertEqual(int(short_jdepth["original_pos_0based"]), 0)
            with (out / "evidence_summary.tsv").open() as handle:
                summary = list(csv.DictReader(handle, delimiter="\t"))
            short = next(r for r in summary if r["track"] == "short")
            long = next(r for r in summary if r["track"] == "long")
            self.assertIn("junction_window_mean_depth", short)
            self.assertEqual(int(short["bridge_unique_fragments"]), 1)
            self.assertAlmostEqual(float(short["genome_mean_depth"]), 8 / 101, places=6)
            self.assertAlmostEqual(float(long["genome_mean_depth"]), 10 / 101, places=6)
            bad = [r for r in self.rows(out / "junction_reads.tsv.gz") if r["qname"] == "pair_bad_mate"]
            self.assertTrue(bad)
            self.assertTrue(all(str(r["bridges_junction"]).lower() == "false" for r in bad))
            manifest = json.loads((out / "manifest.json").read_text())
            method = manifest.get("method", {})
            self.assertIn("10.1038/s41587-025-02971-8", method.get("reference", ""))
            self.assertTrue((out / "zero_coverage_regions.tsv").exists())
            self.assertTrue((out / "long_clipping.tsv").exists())
            with (out / "long_clipping.tsv").open() as handle:
                clipping = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({r["side"] for r in clipping}, {"left", "right"})
            self.assertTrue(all(float(r["clipping_ratio"]) > 0 for r in clipping))

    def test_missing_input_fails_without_manifest(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "evidence"
            args = ["--contig", "ctg1", "--contig-fa", str(Path(td) / "missing.fa"), "--long-bam", str(Path(td) / "a.bam"), "--long-junction-bam", str(Path(td) / "b.bam"), "--short-bam", str(Path(td) / "c.bam"), "--short-rotated-bam", str(Path(td) / "d.bam"), "--rotation-coords", str(Path(td) / "coords.tsv"), "--out-dir", str(out)]
            out.mkdir(parents=True)
            (out / "manifest.json").write_text("stale")
            result = subprocess.run([PYTHON, str(ROOT / "scripts/export_circular_evidence.py")] + args, env=os.environ.copy(), capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((out / "manifest.json").exists())

    def test_empty_bams_are_measured_as_zero(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            contig, coords, long_bam, long_junction, short_bam, rotated = self.make_inputs(td)
            for path, name, length in ((long_bam, "ctg1", 101), (short_bam, "ctg1", 101), (rotated, "ctg1", 101), (long_junction, "ctg1_junction", 100)):
                self.write_bam(path, name, length, [])
            out = self.run_export(td, setup=False)
            with (out / "evidence_summary.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(all(float(r["genome_mean_depth"]) == 0 for r in rows))
            self.assertEqual(self.rows(out / "junction_reads.tsv.gz"), [])

    def test_bad_reference_length_and_rotation_fail(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            contig, coords, long_bam, long_junction, short_bam, rotated = self.make_inputs(td)
            bad = td / "bad.bam"
            self.write_bam(bad, "ctg1", 99, [])
            out = td / "evidence"
            cmd = [PYTHON, str(ROOT / "scripts/export_circular_evidence.py"), "--contig", "ctg1", "--contig-fa", str(contig), "--long-bam", str(bad), "--long-junction-bam", str(long_junction), "--short-bam", str(short_bam), "--short-rotated-bam", str(rotated), "--rotation-coords", str(coords), "--out-dir", str(out)]
            self.assertNotEqual(subprocess.run(cmd, env=os.environ.copy(), capture_output=True).returncode, 0)
