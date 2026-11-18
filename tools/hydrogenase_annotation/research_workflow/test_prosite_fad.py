import hashlib
import sys
import tempfile
import unittest
from unittest.mock import patch
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prosite_fad import _parse_output, _read_fasta, scan

ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / "references/raw/round65_prosite/officialprofile/PS51384.txt"
SCANNER = ROOT / ".env/bin/pfsearch2.3"
SEQ = "MMLPKEIMMPNDNPYALHRVKVLKVYSLTETEKLFLFRFEDPELAEKWTFKPGQFVQLTIPGVGEVPISICSSPMRKGFFELCIRKAGRVTTVVHRLKPGDTVLVRGPYGNGFPVDEWEGMDLLLIAAGLGTAPLRSVFLYAMDNRWKYGNITFINTARYGKDLLFYKELEAMKDLAEAENVKIIQSVTRDPNWPGLKGRPQQFIVEANTNPKNTAVAICGPPRMYKSVFEALINYGYRPENIFVTLERRMKCGIGKCGHCNVGTSTSWKYICKDGPVFTYFDIVSTPGLLD"


class PrositeFadTests(unittest.TestCase):
    def model(self):
        return {"prosite_fad": {"profile": str(PROFILE.relative_to(ROOT)), "scanner": str(SCANNER.relative_to(ROOT)),
            "normalization_mode": 1, "min_normalized_score": 8.5, "min_raw_score": 1658},
            "sources": {str(PROFILE.relative_to(ROOT)): hashlib.sha256(PROFILE.read_bytes()).hexdigest(),
                        str(SCANNER.relative_to(ROOT)): hashlib.sha256(SCANNER.read_bytes()).hexdigest()}}

    def test_reference_and_long_id_are_mapped_and_have_non_one_coordinates(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); annotation = td / "annotation"; annotation.mkdir()
            (annotation / "proteins.faa").write_text(">this_is_a_deliberately_long_original_protein_identifier_123456789\n" + SEQ + "\n")
            rows, bound = scan(annotation, td / "out", self.model(), ROOT)
            self.assertEqual(len(rows), 1); self.assertEqual(rows[0]["protein_id"], "this_is_a_deliberately_long_original_protein_identifier_123456789")
            self.assertGreater(rows[0]["query_start"], 1); self.assertEqual(rows[0]["level"], 0); self.assertTrue(bound)
            self.assertTrue(all(Path(p).is_file() and hashlib.sha256(Path(p).read_bytes()).hexdigest() == digest for p, digest in bound.items()))

    def test_duplicate_ids_fail_closed(self):
        with tempfile.TemporaryDirectory() as td:
            annotation = Path(td) / "annotation"; annotation.mkdir()
            (annotation / "proteins.faa").write_text(">x\nMKK\n>x\nMKK\n")
            with self.assertRaises(ValueError): scan(annotation, Path(td) / "out", self.model(), ROOT)

    def test_unknown_id_coordinates_cutoff_and_motif_fail_closed(self):
        seqs = {"p00000001": SEQ}
        sha = hashlib.sha256(PROFILE.read_bytes()).hexdigest(); scanner = hashlib.sha256(SCANNER.read_bytes()).hexdigest()
        base = ">p00000001/15-115 motif=PS51384|FAD_FR norm_score=16.5 raw_score=3000 level=0 seq_end=-178 motif_start=1 motif_end=-1"
        for text in [base.replace("p00000001", "unknown"), base.replace("15-115", "0-115"), base.replace("norm_score=16.5", "norm_score=8.4"), base.replace("PS51384|FAD_FR", "WRONG"), base + " trailer", base.replace("motif_start=1", "motif_start=0"), base.replace("seq_end=-178", "seq_end=-177")]:
            with self.assertRaises(ValueError): _parse_output(text, {"p00000001": "q"}, {"q": SEQ}, sha, scanner, (1658, 8.5, 102))

    def test_iupac_unknown_symbol_is_preserved_for_local_scanner(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); annotation = td / "annotation"; annotation.mkdir()
            seq = SEQ[:40] + 'XUBZJO' + SEQ[46:]
            (annotation / "proteins.faa").write_text(">with_unknowns\n" + seq + "\n")
            rows, _ = scan(annotation, td / 'out', self.model(), ROOT)
            self.assertEqual(_read_fasta(annotation / "proteins.faa")["with_unknowns"], seq)
            self.assertEqual(rows[0]['query_length'], len(seq))

    def test_source_change_during_scan_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); annotation = td / 'annotation'; annotation.mkdir()
            source = annotation / 'proteins.faa'; source.write_text('>q\n' + SEQ + '\n')
            real_run = subprocess.run
            def mutate(*args, **kwargs):
                result = real_run(*args, **kwargs)
                source.write_text('>changed\nMKK\n')
                return result
            with patch('prosite_fad.subprocess.run', side_effect=mutate), self.assertRaisesRegex(ValueError, 'changed during scan'):
                scan(annotation, td / 'out', self.model(), ROOT)

    def test_long_protein_is_not_truncated_by_fortran_line_buffer(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); annotation = td / 'annotation'; annotation.mkdir()
            seq = 'M' * 600 + SEQ + 'M' * 300
            (annotation / 'proteins.faa').write_text('>long\n' + seq + '\n')
            rows, _ = scan(annotation, td / 'out', self.model(), ROOT)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['query_start'], 615)
            self.assertEqual(rows[0]['query_end'], 715)
            self.assertEqual(rows[0]['query_length'], len(seq))


if __name__ == "__main__":
    unittest.main()
