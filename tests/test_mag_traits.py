from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

from scripts.merge_mag_traits import (
    parse_drep_cdb,
    parse_drep_clusters,
    parse_all_dbcan_sources,
    parse_all_camper_sources,
    parse_all_grodon2_sources,
    parse_all_genomespot_sources,
)


class MagTraitsTests(unittest.TestCase):
    def test_parse_two_round_drep_clusters(self):
        with tempfile.TemporaryDirectory() as td:
            r1_cdb = Path(td) / "r1_cdb.csv"
            r2_cdb = Path(td) / "r2_cdb.csv"

            # Round 1: MAG1 and MAG2 in cluster 1_0 (winner MAG1), MAG3 in 1_1 (winner MAG3)
            r1_cdb.write_text(
                "genome,secondary_cluster\n"
                "MAG1.fa,1_0\n"
                "MAG2.fa,1_0\n"
                "MAG3.fa,1_1\n"
            )

            # Round 2: MAG1 and MAG3 clustered together into species cluster 2_0
            r2_cdb.write_text(
                "genome,secondary_cluster\n"
                "MAG1.fa,2_0\n"
                "MAG3.fa,2_0\n"
            )

            mapping = parse_drep_clusters(str(r1_cdb), str(r2_cdb), None)
            self.assertIn("MAG1", mapping)
            self.assertIn("MAG3", mapping)
            self.assertEqual(mapping["MAG1"]["cluster"], "R2_2_0")
            self.assertEqual(mapping["MAG1"]["size"], 3)  # MAG1, MAG2, MAG3
            self.assertCountEqual(mapping["MAG1"]["members"], ["MAG1", "MAG2", "MAG3"])

    def test_parse_batch_sources(self):
        with tempfile.TemporaryDirectory() as td:
            results_dir = Path(td)
            batch1 = results_dir / "batches" / "batch_0001"
            (batch1 / "dbcan").mkdir(parents=True)
            (batch1 / "camper").mkdir(parents=True)
            (batch1 / "grodon2").mkdir(parents=True)
            (batch1 / "genomespot").mkdir(parents=True)

            # dbcan overview
            (batch1 / "dbcan" / "overview.txt").write_text(
                "Gene ID\tHMMER\tDiamond\tdbCAN_sub\t#ofTools\n"
                "MAG001___gene1\tGH13(1-100)\tGH13(1-100)\t-\t2\n"
                "MAG001___gene2\tGT2(1-50)\tGT2(1-50)\t-\t2\n"
            )

            # camper hits
            (batch1 / "camper" / "camper_hits.tsv").write_text(
                "MAG\tGene_ID\tTarget_HMM\tE_value\tScore\n"
                "MAG001\tgene3\tLaccase\t1e-10\t55.2\n"
            )

            # grodon2
            (batch1 / "grodon2" / "grodon2_growth.tsv").write_text(
                "MAG\tDoublingTime_hours\tMaxGrowthRate_per_hour\tCUB_HEG\tGrowthCategory\tMethod\tMode\tTemperature\tStatus\n"
                "MAG001\t1.5\t0.46\t0.62\tFast\tgRodon2\tpartial\t37\tSuccess\n"
            )

            # genomespot
            (batch1 / "genomespot" / "genomespot_summary.tsv").write_text(
                "MAG\tOGT\tTmin\tTmax\tpH_opt\tpH_min\tpH_max\tSalinity_opt\tSalinity_min\tSalinity_max\tOxygen_tolerance\tOxygen_prob\tStatus\n"
                "MAG001\t37.2\t10.5\t45.0\t7.2\t5.5\t8.5\t1.2\t0.0\t3.5\tAerobic\t0.98\tOK\n"
            )

            dbcan_data = parse_all_dbcan_sources(str(results_dir))
            self.assertIn("MAG001", dbcan_data)
            tot_caz, caz_classes, top_fams = dbcan_data["MAG001"]
            self.assertEqual(tot_caz, 2)
            self.assertEqual(caz_classes.get("GH"), 1)
            self.assertEqual(caz_classes.get("GT"), 1)

            camper_data = parse_all_camper_sources(str(results_dir))
            self.assertIn("MAG001", camper_data)
            self.assertEqual(camper_data["MAG001"][0], 1)
            self.assertIn("Laccase", camper_data["MAG001"][1])

            grodon_data = parse_all_grodon2_sources(str(results_dir))
            self.assertIn("MAG001", grodon_data)
            self.assertEqual(grodon_data["MAG001"]["DoublingTime_hours"], "1.5")
            self.assertEqual(grodon_data["MAG001"]["GrowthCategory"], "Fast")

            gspot_data = parse_all_genomespot_sources(str(results_dir))
            self.assertIn("MAG001", gspot_data)
            self.assertEqual(gspot_data["MAG001"]["OGT"], "37.2")
            self.assertEqual(gspot_data["MAG001"]["pH_opt"], "7.2")
            self.assertEqual(gspot_data["MAG001"]["Oxygen_tolerance"], "Aerobic")
            self.assertEqual(gspot_data["MAG001"]["Status"], "OK")

    def test_full_merge_cli(self):
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            results_dir = Path(td)
            batch1 = results_dir / "batches" / "batch_0001"
            (batch1 / "dbcan").mkdir(parents=True)
            (batch1 / "camper").mkdir(parents=True)
            (batch1 / "grodon2").mkdir(parents=True)
            (batch1 / "genomespot").mkdir(parents=True)
            (batch1 / "mag_list.tsv").write_text("MAG001\t/path/to/MAG001.fa\n")

            (batch1 / "dbcan" / "overview.txt").write_text(
                "Gene ID\tHMMER\tDiamond\tdbCAN_sub\t#ofTools\n"
                "MAG001___g1\tGH13(1-100)\tGH13(1-100)\t-\t2\n"
            )
            (batch1 / "camper" / "camper_hits.tsv").write_text(
                "MAG\tGene_ID\tTarget_HMM\tE_value\tScore\n"
                "MAG001\tg2\tLaccase\t1e-10\t50.0\n"
            )
            (batch1 / "grodon2" / "grodon2_growth.tsv").write_text(
                "MAG\tDoublingTime_hours\tMaxGrowthRate_per_hour\tCUB_HEG\tGrowthCategory\tMethod\tMode\tTemperature\tStatus\n"
                "MAG001\t2.0\t0.35\t0.58\tFast\tgRodon2\tpartial\t37\tSuccess\n"
            )
            (batch1 / "genomespot" / "genomespot_summary.tsv").write_text(
                "MAG\tOGT\tTmin\tTmax\tpH_opt\tpH_min\tpH_max\tSalinity_opt\tSalinity_min\tSalinity_max\tOxygen_tolerance\tOxygen_prob\tStatus\n"
                "MAG001\t35.0\t15.0\t42.0\t7.0\t6.0\t8.0\t0.5\t0.0\t2.0\tAerobic\t0.95\tOK\n"
            )

            out_tsv = results_dir / "summary.tsv"
            cmd = [
                "python3", "scripts/merge_mag_traits.py",
                "--results-dir", str(results_dir),
                "--output", str(out_tsv)
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, msg=res.stderr)
            self.assertTrue(out_tsv.exists())
            lines = out_tsv.read_text().splitlines()
            self.assertGreaterEqual(len(lines), 2)
            self.assertIn("MAG001", lines[1])
            self.assertIn("Fast", lines[1])
            self.assertIn("Aerobic", lines[1])
            # Check exactly 29 columns when dRep is not supplied
            header_cols = lines[0].split("\t")
            self.assertEqual(len(header_cols), 29)
            data_cols = lines[1].split("\t")
            self.assertEqual(len(data_cols), 29)

    def test_corrupted_grodon_error_message(self):
        """Ensure multi-space or quoted error messages with quotes/newlines do not shift columns."""
        with tempfile.TemporaryDirectory() as td:
            results_dir = Path(td)
            batch1 = results_dir / "batches" / "batch_0001"
            (batch1 / "grodon2").mkdir(parents=True)
            (batch1 / "mag_list.tsv").write_text("MAG_ERR\t/path/to/MAG_ERR.fa\n")

            cdb_file = results_dir / "Cdb.csv"
            cdb_file.write_text("genome,secondary_cluster\nMAG_ERR.fa,1_0\n")

            # Emulate the exact error line the user had:
            (batch1 / "grodon2" / "grodon2_growth.tsv").write_text(
                "MAG\tDoublingTime_hours\tMaxGrowthRate_per_hour\tCUB_HEG\tGrowthCategory\tMethod\tMode\tTemperature\tStatus\n"
                'MAG_ERR\tNA\tNA\tNA\tUnknown\tgRodon\tError: Invalid mode. Please pick an available prediction mode ("full", "partial", "metagenome_v1")\n'
            )

            out_tsv = results_dir / "summary.tsv"
            import subprocess
            cmd = [
                "python3", "scripts/merge_mag_traits.py",
                "--results-dir", str(results_dir),
                "--output", str(out_tsv),
                "--cdb", str(cdb_file)
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, msg=res.stderr)
            lines = out_tsv.read_text().splitlines()
            self.assertEqual(len(lines), 2)
            row = lines[1].split("\t")
            self.assertEqual(len(row), 32)
            self.assertEqual(row[0], "MAG_ERR")
            self.assertEqual(row[1], "1_0")  # dRep_Cluster
            self.assertEqual(row[30], "gRodon")  # Method
            self.assertIn("Invalid mode", row[31])  # Status

    def test_run_genomespot_runner_and_parsing(self):
        """Test GenomeSPOT runner CLI and predictions parser."""
        from scripts.run_genomespot import parse_predictions_tsv
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            pred_tsv = td_path / "test.predictions.tsv"
            pred_tsv.write_text(
                "target\tvalue\terror\tunits\tis_novel\twarning\n"
                "temperature_optimum\t37.5\t3.2\tdegrees_c\tFalse\t\n"
                "ph_optimum\t6.8\t0.5\tpH\tFalse\t\n"
                "salinity_optimum\t1.5\t0.3\tpercent\tFalse\t\n"
                "oxygen\tFacultative\t0.88\tprobability\tFalse\t\n"
            )
            parsed = parse_predictions_tsv(str(pred_tsv))
            self.assertEqual(parsed["OGT"], "37.5")
            self.assertEqual(parsed["pH_opt"], "6.8")
            self.assertEqual(parsed["Salinity_opt"], "1.5")
            self.assertEqual(parsed["Oxygen_tolerance"], "Facultative")
            self.assertEqual(parsed["Oxygen_prob"], "0.88")
            self.assertEqual(parsed["Status"], "OK")

            # Test CLI without models: should report Missing_Models_Dir
            mag_list = td_path / "mag_list.tsv"
            mag_list.write_text("MAG_MISSING\t/nonexistent/path/MAG_MISSING.fna\n")
            faa_dir = td_path / "faa"
            faa_dir.mkdir()
            out_dir = td_path / "out_gspot"

            import subprocess
            cmd = [
                "python3", "scripts/run_genomespot.py",
                "--mag_list", str(mag_list),
                "--faa_dir", str(faa_dir),
                "--out_dir", str(out_dir),
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, msg=res.stderr)
            summary_file = out_dir / "genomespot_summary.tsv"
            self.assertTrue(summary_file.exists())
            summary_lines = summary_file.read_text().splitlines()
            self.assertEqual(len(summary_lines), 2)
            self.assertIn("MAG_MISSING", summary_lines[1])
            self.assertIn("Failed(Missing_Models_Dir)", summary_lines[1])

            # Test CLI with dummy models directory: should proceed to check inputs and report Missing_FAA
            dummy_models = td_path / "models"
            dummy_models.mkdir()
            (dummy_models / "instructions.json").write_text("{}")
            out_dir2 = td_path / "out_gspot2"

            cmd2 = [
                "python3", "scripts/run_genomespot.py",
                "--mag_list", str(mag_list),
                "--faa_dir", str(faa_dir),
                "--out_dir", str(out_dir2),
                "--models", str(dummy_models),
            ]
            res2 = subprocess.run(cmd2, capture_output=True, text=True)
            self.assertEqual(res2.returncode, 0, msg=res2.stderr)
            summary_file2 = out_dir2 / "genomespot_summary.tsv"
            self.assertTrue(summary_file2.exists())
            summary_lines2 = summary_file2.read_text().splitlines()
            self.assertEqual(len(summary_lines2), 2)
            self.assertIn("MAG_MISSING", summary_lines2[1])
            self.assertIn("Missing_FAA", summary_lines2[1])

    def test_sample_namespacing_and_config(self):
        """Verify sample configuration resolution logic matches repo conventions."""
        # Simulated config with samples dictionary
        cfg1 = {
            "samples": {
                "gses_all_rawMAGs_fa_bac": {
                    "mag_dir": "data/gses_all_rawMAGs_fa_bac",
                    "genome_info": "data/gses_all_rawMAGs_fa_V0.genomeInfo.csv"
                },
                "gses_all_rawMAGs_fa_ar": {
                    "mag_dir": "data/gses_all_rawMAGs_fa_ar"
                }
            }
        }

        # 1. Active samples from dictionary
        def resolve_active(cfg):
            cli_val = cfg.get("sample") or cfg.get("sample_name")
            if isinstance(cli_val, str) and cli_val.strip() and cli_val.strip().lower() not in ["all", "*"]:
                return [cli_val.strip()]
            for cand in [cfg.get("samples"), cfg.get("sample")]:
                if isinstance(cand, dict) and cand:
                    return list(cand.keys())
                elif isinstance(cand, list) and cand:
                    return [str(x).strip() for x in cand if str(x).strip()]
                elif isinstance(cand, str) and cand.strip() and cand.strip().lower() not in ["all", "*"]:
                    return [cand.strip()]
            m_dir = cfg.get("mag_dir", "")
            if m_dir:
                base = os.path.basename(os.path.normpath(str(m_dir).strip()))
                if base:
                    return [base]
            return ["default"]

        self.assertEqual(
            resolve_active(cfg1),
            ["gses_all_rawMAGs_fa_bac", "gses_all_rawMAGs_fa_ar"]
        )

        # 2. CLI override via sample: "gses_all_rawMAGs_fa_bac"
        cfg2 = dict(cfg1)
        cfg2["sample"] = "gses_all_rawMAGs_fa_bac"
        self.assertEqual(resolve_active(cfg2), ["gses_all_rawMAGs_fa_bac"])

        # 3. Fallback to mag_dir basename when samples dict is absent
        cfg3 = {"mag_dir": "data/gses_all_rawMAGs_fa_bac"}
        self.assertEqual(resolve_active(cfg3), ["gses_all_rawMAGs_fa_bac"])

        # 4. Fallback to default
        cfg4 = {}
        self.assertEqual(resolve_active(cfg4), ["default"])

    def test_mag_traits_drep_utils(self):
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            base_dir = Path(td)
            
            # 1. Test pool-manifest
            b0 = base_dir / "block_000" / "dereplicated_genomes"
            b1 = base_dir / "block_001" / "dereplicated_genomes"
            b0.mkdir(parents=True)
            b1.mkdir(parents=True)
            (b0 / "genome_A.fa").write_text(">A\nATGC\n")
            (b0 / "genome_B.fasta").write_text(">B\nATGC\n")
            (b1 / "genome_C.fna").write_text(">C\nATGC\n")

            manifest_out = base_dir / "pooled_manifest.tsv"
            cmd_pool = [
                "python3", "scripts/mag_traits_drep_utils.py", "pool-manifest",
                "--out-manifest", str(manifest_out),
                "--block-dirs", str(b0), str(b1)
            ]
            res = subprocess.run(cmd_pool, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, msg=res.stderr)
            self.assertTrue(manifest_out.exists())
            lines = manifest_out.read_text().strip().splitlines()
            self.assertEqual(len(lines), 3)
            ids = [l.split("\t")[0] for l in lines]
            self.assertIn("genome_A", ids)
            self.assertIn("genome_B", ids)
            self.assertIn("genome_C", ids)

            # 2. Test merge-genome-info
            csv1 = base_dir / "ginfo1.csv"
            csv2 = base_dir / "ginfo2.csv"
            csv1.write_text("genome,completeness,contamination\ngenome_A.fa,95.0,1.2\n")
            csv2.write_text("genome,completeness,contamination\ngenome_C.fa,88.5,0.5\n")
            merged_csv = base_dir / "merged_ginfo.csv"

            cmd_merge = [
                "python3", "scripts/mag_traits_drep_utils.py", "merge-genome-info",
                "--out", str(merged_csv),
                "--inputs", str(csv1), str(csv2)
            ]
            res = subprocess.run(cmd_merge, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, msg=res.stderr)
            self.assertTrue(merged_csv.exists())
            m_lines = merged_csv.read_text().strip().splitlines()
            self.assertEqual(len(m_lines), 3)  # 1 header + 2 data rows
            self.assertEqual(m_lines[0], "genome,completeness,contamination")

            # 3. Test filter-genome-info
            derep_dir = base_dir / "round1_derep"
            derep_dir.mkdir()
            (derep_dir / "genome_A.fa").write_text(">A\n")

            filtered_csv = base_dir / "round2_ginfo.csv"
            cmd_filter = [
                "python3", "scripts/mag_traits_drep_utils.py", "filter-genome-info",
                "--derep-dir", str(derep_dir),
                "--src", str(merged_csv),
                "--dst", str(filtered_csv)
            ]
            res = subprocess.run(cmd_filter, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, msg=res.stderr)
            self.assertTrue(filtered_csv.exists())
            f_lines = filtered_csv.read_text().strip().splitlines()
            self.assertEqual(len(f_lines), 2)  # header + genome_A
            self.assertIn("genome_A.fa,95.0,1.2", f_lines[1])

    def test_partition_manifest_taxonomy(self):
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            base_dir = Path(td)
            manifest = base_dir / "manifest.tsv"
            taxonomy = base_dir / "taxonomy.tsv"
            out_dir = base_dir / "partitions"
            out_list = base_dir / "block_list.txt"

            manifest.write_text(
                "MAG1\t/path/MAG1.fa\n"
                "MAG2\t/path/MAG2.fa\n"
                "MAG3\t/path/MAG3.fa\n"
                "MAG4\t/path/MAG4.fa\n"
                "MAG5\t/path/MAG5.fa\n"
            )

            taxonomy.write_text(
                "bin\tcompleteness\tcontamination\tGC\tN50\tsize\tclassification\n"
                "MAG1\t90.0\t1.0\t0.5\t5000\t2000000\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales\n"
                "MAG2\t85.0\t2.0\t0.5\t6000\t2100000\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales\n"
                "MAG3\t92.0\t0.5\t0.6\t7000\t3000000\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales\n"
                "MAG4\t80.0\t3.0\t0.4\t4000\t1800000\td__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhizobiales\n"
                "MAG5\t75.0\t4.0\t0.4\t3000\t1500000\tUnclassified\n"
            )

            cmd = [
                "python3", "scripts/partition_manifest.py",
                "--manifest", str(manifest),
                "--taxonomy", str(taxonomy),
                "--block-size", "1",
                "--out-dir", str(out_dir),
                "--out-list", str(out_list)
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, msg=res.stderr)
            self.assertTrue(out_list.exists())
            blocks = out_list.read_text().strip().splitlines()
            for b in blocks:
                self.assertTrue(b.startswith("tax_"), f"Block {b} should start with tax_")
                self.assertTrue((out_dir / f"{b}.tsv").exists())

            out_list10 = base_dir / "block_list10.txt"
            cmd10 = [
                "python3", "scripts/partition_manifest.py",
                "--manifest", str(manifest),
                "--taxonomy", str(taxonomy),
                "--block-size", "10",
                "--min-phylum-size", "0",
                "--out-dir", str(out_dir),
                "--out-list", str(out_list10)
            ]
            res10 = subprocess.run(cmd10, capture_output=True, text=True)
            self.assertEqual(res10.returncode, 0, msg=res10.stderr)
            blocks10 = out_list10.read_text().strip().splitlines()
            self.assertIn("tax_p__Bacteroidota", blocks10)
            self.assertIn("tax_p__Pseudomonadota", blocks10)
            self.assertIn("tax_unclassified", blocks10)

            # Test min_phylum_size = 3 (Bacteroidota has 2, Pseudomonadota has 2 -> both aggregated into tax_minor_phyla)
            out_list_agg = base_dir / "block_list_agg.txt"
            cmd_agg = [
                "python3", "scripts/partition_manifest.py",
                "--manifest", str(manifest),
                "--taxonomy", str(taxonomy),
                "--block-size", "10",
                "--min-phylum-size", "3",
                "--out-dir", str(out_dir),
                "--out-list", str(out_list_agg)
            ]
            res_agg = subprocess.run(cmd_agg, capture_output=True, text=True)
            self.assertEqual(res_agg.returncode, 0, msg=res_agg.stderr)
            blocks_agg = out_list_agg.read_text().strip().splitlines()
            self.assertIn("tax_minor_phyla", blocks_agg)
            self.assertIn("tax_unclassified", blocks_agg)
            self.assertNotIn("tax_p__Bacteroidota", blocks_agg)

    def test_union_blocks(self):
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            base_dir = Path(td)
            b_list = base_dir / "block_list.txt"
            b_list.write_text("tax_p__Bacteroidota\ntax_p__Pseudomonadota\n")

            # Setup block 1
            b1_dir = base_dir / "tax_p__Bacteroidota"
            (b1_dir / "dereplicated_genomes").mkdir(parents=True)
            (b1_dir / "data_tables").mkdir(parents=True)
            (b1_dir / "dereplicated_genomes" / "MAG1.fa").write_text(">MAG1\nATGC\n")
            (b1_dir / "data_tables" / "Cdb.csv").write_text(
                "genome,secondary_cluster,primary_cluster\n"
                "MAG1.fa,1_0,1\n"
                "MAG2.fa,1_0,1\n"
            )
            (b1_dir / "data_tables" / "Widb.csv").write_text(
                "genome,score,cluster\n"
                "MAG1.fa,98.5,1_0\n"
            )
            (b1_dir / "genomeInfo.csv").write_text(
                "genome,completeness,contamination\n"
                "MAG1.fa,90.0,1.0\n"
                "MAG2.fa,85.0,2.0\n"
            )

            # Setup block 2
            b2_dir = base_dir / "tax_p__Pseudomonadota"
            (b2_dir / "dereplicated_genomes").mkdir(parents=True)
            (b2_dir / "data_tables").mkdir(parents=True)
            (b2_dir / "dereplicated_genomes" / "MAG3.fa").write_text(">MAG3\nGGCC\n")
            (b2_dir / "data_tables" / "Cdb.csv").write_text(
                "genome,secondary_cluster,primary_cluster\n"
                "MAG3.fa,1_0,1\n"
                "MAG4.fa,1_0,1\n"
            )
            (b2_dir / "data_tables" / "Widb.csv").write_text(
                "genome,score,cluster\n"
                "MAG3.fa,95.0,1_0\n"
            )
            (b2_dir / "genomeInfo.csv").write_text(
                "genome,completeness,contamination\n"
                "MAG3.fa,92.0,0.5\n"
                "MAG4.fa,80.0,3.0\n"
            )

            out_derep = base_dir / "final_derep"
            out_cdb = base_dir / "final_cdb.csv"
            out_widb = base_dir / "final_widb.csv"
            out_ginfo = base_dir / "final_ginfo.csv"

            cmd = [
                "python3", "scripts/mag_traits_drep_utils.py", "union-blocks",
                "--block-list", str(b_list),
                "--blocks-dir", str(base_dir),
                "--out-derep", str(out_derep),
                "--out-cdb", str(out_cdb),
                "--out-widb", str(out_widb),
                "--out-ginfo", str(out_ginfo)
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, msg=res.stderr)

            # Check dereplicated genomes unioned
            self.assertTrue((out_derep / "MAG1.fa").exists())
            self.assertTrue((out_derep / "MAG3.fa").exists())

            # Check Cdb namespacing
            cdb_lines = out_cdb.read_text().strip().splitlines()
            self.assertEqual(len(cdb_lines), 5)  # 1 header + 4 rows
            self.assertIn("MAG1.fa,tax_p__Bacteroidota_1_0,tax_p__Bacteroidota_1", cdb_lines)
            self.assertIn("MAG3.fa,tax_p__Pseudomonadota_1_0,tax_p__Pseudomonadota_1", cdb_lines)

            # Check Widb namespacing
            widb_lines = out_widb.read_text().strip().splitlines()
            self.assertEqual(len(widb_lines), 3)  # 1 header + 2 winners
            self.assertIn("MAG1.fa,98.5,tax_p__Bacteroidota_1_0", widb_lines)
            self.assertIn("MAG3.fa,95.0,tax_p__Pseudomonadota_1_0", widb_lines)

            # Check genomeInfo
            ginfo_lines = out_ginfo.read_text().strip().splitlines()
            self.assertEqual(len(ginfo_lines), 5)  # 1 header + 4 genomes


if __name__ == "__main__":
    unittest.main()


