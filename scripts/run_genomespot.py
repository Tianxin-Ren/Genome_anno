#!/usr/bin/env python3
"""
Batch runner for GenomeSPOT phenotypic trait prediction on MAG batches.

Usage:
    python3 run_genomespot.py \
        --mag_list <mag_list.tsv> \
        --faa_dir <faa_dir> \
        --out_dir <out_dir> \
        [--models <models_dir>] \
        [--threads <threads>] \
        [--save_features]
"""

import os
import sys
import argparse
import subprocess
import tempfile
import gzip
import shutil
import re
from pathlib import Path


def find_models_dir(user_models_dir=None):
    """
    Discovers GenomeSPOT trained models directory using multiple fallbacks:
    1. User-supplied CLI argument
    2. Environment variable GENOMESPOT_MODELS
    3. Conda environment paths ($CONDA_PREFIX/share/genomespot/models, etc.)
    4. Python package internal models (genome_spot package dir)
    5. Standard HPC cluster database paths
    """
    candidates = []

    if user_models_dir:
        u_path = os.path.abspath(os.path.expanduser(user_models_dir))
        candidates.extend([
            u_path,
            os.path.join(u_path, "models"),
            os.path.join(u_path, "GenomeSPOT", "models"),
        ])

    env_var = os.environ.get("GENOMESPOT_MODELS")
    if env_var:
        e_path = os.path.abspath(os.path.expanduser(env_var))
        candidates.extend([
            e_path,
            os.path.join(e_path, "models"),
        ])

    conda_prefix = os.environ.get("CONDA_PREFIX", sys.prefix)
    if conda_prefix:
        candidates.extend([
            os.path.join(conda_prefix, "share", "genomespot", "models"),
            os.path.join(conda_prefix, "share", "GenomeSPOT", "models"),
            os.path.join(conda_prefix, "models"),
        ])

    # Try python package directory
    try:
        import genome_spot
        pkg_dir = os.path.dirname(genome_spot.__file__)
        candidates.extend([
            os.path.join(pkg_dir, "models"),
            os.path.join(pkg_dir, "..", "models"),
        ])
    except ImportError:
        pass

    # Standard cluster paths
    candidates.extend([
        "/dssg/home/acct-schhui/gses0/Database/GenomeSPOT/models",
        "/dssg/home/acct-schhui/gses0/Database/GenomeSPOT/GenomeSPOT/models",
        "/dssg/home/acct-schhui/gses0/Database/GenomeSPOT",
        "/dssg/home/acct-clswfp/clswfp-icdli1/Database/GenomeSPOT/models",
        "/lustre/home/acct-clswfp/clswfp-icdli1/DataBase/GenomeSPOT/models",
    ])

    for p in candidates:
        if not p or not os.path.isdir(p):
            continue
        # Direct check for instructions.json
        if os.path.isfile(os.path.join(p, "instructions.json")):
            return os.path.abspath(p)
        if os.path.isfile(os.path.join(p, "models", "instructions.json")):
            return os.path.abspath(os.path.join(p, "models"))
        # Fallback check for any .joblib or instruction files
        try:
            files = os.listdir(p)
            if any("instruction" in f.lower() or f.endswith(".joblib") or f.endswith(".pkl") for f in files):
                return os.path.abspath(p)
        except OSError:
            pass

    return None


def parse_predictions_tsv(pred_tsv):
    """
    Parses GenomeSPOT <prefix>.predictions.tsv:
    Columns: target \t value \t error \t units \t is_novel \t warning
    """
    data = {
        "OGT": "NA",
        "Tmin": "NA",
        "Tmax": "NA",
        "pH_opt": "NA",
        "pH_min": "NA",
        "pH_max": "NA",
        "Salinity_opt": "NA",
        "Salinity_min": "NA",
        "Salinity_max": "NA",
        "Oxygen_tolerance": "NA",
        "Oxygen_prob": "NA",
        "Status": "OK"
    }

    target_map = {
        "temperature_optimum": "OGT",
        "temperature_min": "Tmin",
        "temperature_max": "Tmax",
        "ph_optimum": "pH_opt",
        "ph_min": "pH_min",
        "ph_max": "pH_max",
        "salinity_optimum": "Salinity_opt",
        "salinity_min": "Salinity_min",
        "salinity_max": "Salinity_max",
    }

    if not os.path.exists(pred_tsv) or os.path.getsize(pred_tsv) == 0:
        data["Status"] = "Missing_Output"
        return data

    try:
        with open(pred_tsv, "r", encoding="utf-8", errors="ignore") as f:
            header = f.readline()
            for line in f:
                parts = [p.strip() for p in line.split("\t")]
                if len(parts) < 2:
                    continue
                target = parts[0].lower()
                val = parts[1]

                if target in target_map:
                    data[target_map[target]] = val
                elif "oxygen" in target:
                    # Oxygen row: value might be classification (e.g. Aerobic), error/units might hold probability
                    data["Oxygen_tolerance"] = val
                    if len(parts) >= 3 and parts[2] != "":
                        data["Oxygen_prob"] = parts[2]
    except Exception as e:
        data["Status"] = f"Parse_Error({e})"

    return data


def run_single_mag(mag, fna_path, faa_path, out_prefix, models_dir, save_features=False):
    """
    Executes GenomeSPOT for a single MAG with error isolation.
    """
    # Verify input FAA
    if not os.path.exists(faa_path) or os.path.getsize(faa_path) == 0:
        return {
            "OGT": "NA", "Tmin": "NA", "Tmax": "NA",
            "pH_opt": "NA", "pH_min": "NA", "pH_max": "NA",
            "Salinity_opt": "NA", "Salinity_min": "NA", "Salinity_max": "NA",
            "Oxygen_tolerance": "NA", "Oxygen_prob": "NA",
            "Status": "Missing_FAA"
        }

    # Handle .gz FNA if necessary
    temp_fna = None
    input_fna = fna_path
    if fna_path.endswith(".gz"):
        fd, temp_fna = tempfile.mkstemp(suffix=".fna")
        os.close(fd)
        with gzip.open(fna_path, "rb") as gz_in, open(temp_fna, "wb") as f_out:
            shutil.copyfileobj(gz_in, f_out)
        input_fna = temp_fna

    if not os.path.exists(input_fna) or os.path.getsize(input_fna) == 0:
        if temp_fna and os.path.exists(temp_fna):
            os.remove(temp_fna)
        return {
            "OGT": "NA", "Tmin": "NA", "Tmax": "NA",
            "pH_opt": "NA", "pH_min": "NA", "pH_max": "NA",
            "Salinity_opt": "NA", "Salinity_min": "NA", "Salinity_max": "NA",
            "Oxygen_tolerance": "NA", "Oxygen_prob": "NA",
            "Status": "Missing_FNA"
        }

    cmd = [
        sys.executable, "-m", "genome_spot.genome_spot",
        "--contigs", input_fna,
        "--proteins", faa_path,
        "--output", out_prefix
    ]

    if models_dir:
        cmd.extend(["--models", models_dir])

    # Execute
    res_status = "OK"
    try:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False
        )
        if proc.returncode != 0:
            err_msg = proc.stderr.strip().splitlines()
            if err_msg:
                print(f"ERROR: [GenomeSPOT] MAG {mag} failed with:\n" + "\n".join(err_msg[-10:]), file=sys.stderr)
            last_err = err_msg[-1] if err_msg else f"ExitCode_{proc.returncode}"
            # Clean up error message for status column
            if "FileNotFoundError" in last_err and "instructions.json" in "\n".join(err_msg):
                clean_err = "Missing_instructions.json"
            elif "FileNotFoundError" in last_err:
                match = re.search(r"No such file or directory:\s*['\"]([^'\"]+)['\"]", last_err)
                if match:
                    clean_err = f"FileNotFound({os.path.basename(match.group(1))})"
                else:
                    clean_err = "FileNotFound"
            else:
                clean_err = re.sub(r'[\r\n\t,]+', ' ', last_err)[:50]
            res_status = f"Failed({clean_err})"
    except Exception as e:
        res_status = f"Execution_Error({e})"
    finally:
        if temp_fna and os.path.exists(temp_fna):
            try:
                os.remove(temp_fna)
            except OSError:
                pass

    pred_tsv = f"{out_prefix}.predictions.tsv"
    feat_json = f"{out_prefix}.features.json"

    if not save_features and os.path.exists(feat_json):
        try:
            os.remove(feat_json)
        except OSError:
            pass

    if os.path.exists(pred_tsv) and os.path.getsize(pred_tsv) > 0:
        parsed = parse_predictions_tsv(pred_tsv)
        if res_status != "OK" and parsed["Status"] == "OK":
            parsed["Status"] = res_status
        return parsed
    else:
        # Create an empty placeholder file so downstream checks can see it was attempted
        try:
            with open(pred_tsv, "w", encoding="utf-8") as pf:
                pf.write("target\tvalue\terror\tunits\tis_novel\twarning\n")
        except OSError:
            pass
        return {
            "OGT": "NA",
            "Tmin": "NA",
            "Tmax": "NA",
            "pH_opt": "NA",
            "pH_min": "NA",
            "pH_max": "NA",
            "Salinity_opt": "NA",
            "Salinity_min": "NA",
            "Salinity_max": "NA",
            "Oxygen_tolerance": "NA",
            "Oxygen_prob": "NA",
            "Status": res_status if res_status != "OK" else "No_Predictions"
        }


def main():
    parser = argparse.ArgumentParser(description="Run GenomeSPOT batch trait prediction.")
    parser.add_argument("--mag_list", required=True, help="Path to batch mag_list.tsv (mag\\tfna_path)")
    parser.add_argument("--faa_dir", required=True, help="Path to batch faa directory")
    parser.add_argument("--out_dir", required=True, help="Output directory for GenomeSPOT results")
    parser.add_argument("--models", default=None, help="Optional models directory path")
    parser.add_argument("--threads", type=int, default=2, help="Number of parallel worker threads")
    parser.add_argument("--save_features", action="store_true", help="Save intermediate .features.json")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    models_dir = find_models_dir(args.models)
    if models_dir:
        print(f"INFO: [GenomeSPOT] Using models directory: {models_dir}")
    else:
        err_msg = (
            "================================================================================\n"
            "ERROR: [GenomeSPOT] Pre-trained models directory NOT FOUND!\n"
            "GenomeSPOT requires pre-trained models (specifically 'instructions.json' and associated '*.joblib' files).\n"
            f"Searched user CLI path: {args.models}\n"
            "Please download/copy the models on the cluster via:\n"
            "    git clone https://github.com/cultivarium/GenomeSPOT.git\n"
            "    mkdir -p /dssg/home/acct-schhui/gses0/Database/GenomeSPOT\n"
            "    cp -r GenomeSPOT/models /dssg/home/acct-schhui/gses0/Database/GenomeSPOT/models\n"
            "================================================================================"
        )
        print(err_msg, file=sys.stderr)
        print(err_msg)

    # Read MAG list
    mag_entries = []
    with open(args.mag_list, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                mag_entries.append((parts[0].strip(), parts[1].strip()))
            elif len(parts) == 1:
                mag_entries.append((parts[0].strip(), ""))

    summary_rows = []
    headers = [
        "MAG",
        "OGT",
        "Tmin",
        "Tmax",
        "pH_opt",
        "pH_min",
        "pH_max",
        "Salinity_opt",
        "Salinity_min",
        "Salinity_max",
        "Oxygen_tolerance",
        "Oxygen_prob",
        "Status"
    ]

    print(f"INFO: [GenomeSPOT] Processing {len(mag_entries)} MAGs in batch...")

    for mag, fna_path in mag_entries:
        faa_path = os.path.join(args.faa_dir, f"{mag}.faa")
        out_prefix = os.path.join(args.out_dir, mag)

        if not models_dir:
            pred_res = {
                "OGT": "NA", "Tmin": "NA", "Tmax": "NA",
                "pH_opt": "NA", "pH_min": "NA", "pH_max": "NA",
                "Salinity_opt": "NA", "Salinity_min": "NA", "Salinity_max": "NA",
                "Oxygen_tolerance": "NA", "Oxygen_prob": "NA",
                "Status": "Failed(Missing_Models_Dir)"
            }
        else:
            pred_res = run_single_mag(
                mag=mag,
                fna_path=fna_path,
                faa_path=faa_path,
                out_prefix=out_prefix,
                models_dir=models_dir,
                save_features=args.save_features
            )

        row = [
            mag,
            str(pred_res.get("OGT", "NA")),
            str(pred_res.get("Tmin", "NA")),
            str(pred_res.get("Tmax", "NA")),
            str(pred_res.get("pH_opt", "NA")),
            str(pred_res.get("pH_min", "NA")),
            str(pred_res.get("pH_max", "NA")),
            str(pred_res.get("Salinity_opt", "NA")),
            str(pred_res.get("Salinity_min", "NA")),
            str(pred_res.get("Salinity_max", "NA")),
            str(pred_res.get("Oxygen_tolerance", "NA")),
            str(pred_res.get("Oxygen_prob", "NA")),
            str(pred_res.get("Status", "OK"))
        ]
        summary_rows.append(row)

    summary_file = os.path.join(args.out_dir, "genomespot_summary.tsv")
    with open(summary_file, "w", encoding="utf-8") as sf:
        sf.write("\t".join(headers) + "\n")
        for row in summary_rows:
            sf.write("\t".join(row) + "\n")

    print(f"INFO: [GenomeSPOT] Completed batch. Summary written to {summary_file}")


if __name__ == "__main__":
    main()
