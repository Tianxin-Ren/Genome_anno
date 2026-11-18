import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/query_gca_strains.sh"


class TestQueryGcaStrains(unittest.TestCase):
    def test_cached_json_fields_and_headerless_input(self):
        reports = [
            {"accession": "GCA_000000001.1", "organism": {
                "infraspecificNames": {"strain": "DSM 1234"}}},
            {"accession": "GCA_000000002.1", "assemblyInfo": {"biosample": {
                "description": {"comment": "variant of ATCC:17978"}}}},
            {"accession": "GCA_000000003.1", "assemblyInfo": {"biosample": {
                "attributes": [{"name": "culture_collection", "value": "CCOS:999"}]}}},
            {"accession": "GCA_000000004.1", "assemblyInfo": {"biosample": {
                "sampleName": "only-a-sample-name", "collectedBy": "ATCC"}}},
            {"accession": "GCA_000000005.1", "organism": {
                "infraspecificNames": {"strain": "NCIMB 123; NRRL YB-2076"}}},
            {"accession": "GCA_000000006.1", "organism": {
                "infraspecificNames": {"strain": "DSM 1234"}}},
        ]
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "input.tsv"
            source.write_text("\ufeff" + "\r\n".join(
                " GB_" + report["accession"] + " " for report in reports
            ) + "\r\nGCA_000000007.1\r\n", encoding="utf-8")
            summary = root / "cached.jsonl"
            summary.write_text("".join(json.dumps(report) + "\n" for report in reports))
            env_dir = root / "env"
            (env_dir / "bin").mkdir(parents=True)
            # Simulate a TSV projection omitting the fields carrying collection IDs.
            dataformat = env_dir / "bin/dataformat"
            dataformat.write_text(
                f"#!{sys.executable}\nimport json, sys\n"
                "print('Assembly Accession\\tOrganism Name')\n"
                "for line in sys.stdin:\n"
                "    print(json.loads(line)['accession'] + '\\tExample species')\n"
            )
            dataformat.chmod(0o755)
            output = root / "output"
            process_env = os.environ.copy()
            process_env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + process_env["PATH"]
            result = subprocess.run([
                "sh", str(SCRIPT), "-i", str(source), "-o", str(output),
                "-e", str(env_dir), "--reuse-summary", str(summary),
            ], capture_output=True, text=True, env=process_env)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Input GCA/GCF: 7", result.stdout)
            self.assertIn("Returned metadata GCA: 6", result.stdout)
            self.assertIn("GCA without returned metadata: 1", result.stdout)
            self.assertIn("Returned GCA without strain/isolate ID: 3", result.stdout)
            self.assertIn("GCA with a culture-collection-like ID: 5", result.stdout)
            self.assertIn("GCA with declared culture_collection: 1", result.stdout)
            self.assertEqual((output / "gca_statistics.log").read_text().splitlines()[0],
                             "Input GCA/GCF: 7")
            with (output / "gca_strain_presence.tsv").open() as handle:
                rows = {row["accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["GCA_000000001.1"]["culture_collection_ids"], "DSM:1234")
            self.assertEqual(rows["GCA_000000003.1"]["culture_collection_ids"], "CCOS:999")
            self.assertEqual(rows["GCA_000000004.1"]["strain_or_isolate_ids"], "NA")
            self.assertEqual(rows["GCA_000000004.1"]["culture_collection_ids"], "NA")
            self.assertIn("NCIMB:123", rows["GCA_000000005.1"]["culture_collection_ids"])
            self.assertIn("NRRL:YB-2076", rows["GCA_000000005.1"]["culture_collection_ids"])
            with (output / "gca_collection_evidence.tsv").open() as handle:
                evidence = list(csv.DictReader(handle, delimiter="\t"))
            variant = next(row for row in evidence if row["accession"] == "GCA_000000002.1")
            self.assertEqual(variant["evidence_type"], "text_candidate")
            self.assertEqual(variant["source_text"], "variant of ATCC:17978")
            self.assertIn("comment", variant["source_field"])
            with (output / "duplicated_strain_ids.tsv").open() as handle:
                duplicate = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(duplicate[0]["gca_count"], "2")

    def test_missing_argument_value(self):
        result = subprocess.run(["bash", str(SCRIPT), "-i"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("Missing value", result.stderr)


if __name__ == "__main__":
    unittest.main()
