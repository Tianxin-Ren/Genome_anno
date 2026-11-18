from __future__ import annotations

import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

from common.scripts.runtime import apply_runtime_config
from scripts.unified_annotation_integrate import (
    parse_eggnog,
    parse_interpro,
    parse_kofam,
    parse_mmseqs,
)


def unified_config(steps):
    return {
        "modules": {"unified_annotation": {"enable": True}},
        "unified_annotation": {"steps": steps},
    }


class UnifiedConfigTests(unittest.TestCase):
    def test_disabled_evidence_channels_are_valid_with_consensus(self):
        config = unified_config({
            "cds_prediction": True,
            "protein_similarity": False,
            "orthology": False,
            "ko_profiles": False,
            "protein_domains": False,
            "candidate_integration": True,
        })
        apply_runtime_config(config)
        self.assertEqual(config["unified_annotation"]["logs_dir"], "logs/unified_annotation")
        self.assertEqual(config["unified_annotation"]["tmp_dir"], "work/unified_annotation")

    def test_unimplemented_enabled_step_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "not implemented"):
            apply_runtime_config(unified_config({"crispr_prediction": True}))

    def test_protein_evidence_requires_cds(self):
        with self.assertRaisesRegex(ValueError, "cds_prediction must be true"):
            apply_runtime_config(unified_config({
                "cds_prediction": False,
                "candidate_integration": True,
            }))

    def test_disabled_placeholder_files_add_no_evidence(self):
        evidence = defaultdict(list)
        with tempfile.TemporaryDirectory() as directory:
            placeholder = Path(directory) / "disabled.tsv"
            placeholder.write_text("# evidence channel disabled by config\n")
            for parser in (parse_mmseqs, parse_eggnog, parse_kofam, parse_interpro):
                parser(str(placeholder), evidence)
        self.assertEqual(dict(evidence), {})


if __name__ == "__main__":
    unittest.main()
