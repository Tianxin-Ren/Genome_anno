import unittest
from hox_partner_audit import classify, _blast_loader

class HoxAuditTests(unittest.TestCase):
    def test_statuses_and_incomplete_guard(self):
        self.assertEqual(classify([], [], [], True), "not_detected")
        self.assertEqual(classify([], [{"identity":40}], [], True), "fragment_only")
        self.assertEqual(classify([{"bitscore":100}], [], [], True), "supported")
        self.assertEqual(classify([{"bitscore":100}], [], [{"bitscore":95}], True), "competing_family")
        self.assertEqual(classify([{"bitscore":100},{"bitscore":99}], [], [], True), "pairing_ambiguous")
        self.assertEqual(classify([{"bitscore":100}], [], [], False), "unlinked")
        self.assertEqual(classify([], [], [], True, complete=False), "search_incomplete")

    def test_real_blast_fields_are_named_and_evalue_gated(self):
        fields, load = _blast_loader()
        self.assertEqual(fields[8:12], ["evalue", "bitscore", "qlen", "slen"])
        import tempfile, pathlib
        line="q\tr\t40\t100\t1\t100\t1\t100\t1e-20\t200\t100\t100\tAAAA\tAAAA\n"
        with tempfile.TemporaryDirectory() as d:
            p=pathlib.Path(d)/"b.tsv"; p.write_text(line)
            h=load(p,{"blast_evalue":1e-10,"min_identity":30,"min_query_coverage":.65,"min_reference_coverage":.65},[])[0]
        self.assertTrue(h["passes_thresholds"])
        self.assertEqual(h["qlen"],100); self.assertEqual(h["slen"],100)

    def test_competition_is_same_query_only(self):
        self.assertEqual(classify([{"query":"q1","bitscore":100}], [], [{"query":"q2","bitscore":100}], True), "supported")

if __name__ == "__main__": unittest.main()
