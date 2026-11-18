#!/usr/bin/env python3
"""
serialize_kofam.py
Serializes KOfam HMM profiles + ko_list into database.pkl and generates kofam_manifest.json.
Stores complete provenance:
- release_date
- profiles_sha256
- ko_list_sha256
- database_pkl_sha256
- pykofamsearch_version
- pyhmmer_version
- python_version
- serialization_timestamp
- profile_count & curated_profile_count
"""

import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys

def compute_sha256(filepath):
    if not os.path.isfile(filepath):
        return "file_not_found"
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def count_ko_list(ko_list_path):
    total = 0
    with_threshold = 0
    if os.path.isfile(ko_list_path):
        import gzip
        open_func = gzip.open if ko_list_path.endswith(".gz") else open
        with open_func(ko_list_path, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                total += 1
                if len(parts) >= 2 and parts[1] != "-" and parts[1] != "":
                    with_threshold += 1
    return total, with_threshold

def main():
    parser = argparse.ArgumentParser(description="Serialize KOfam models and generate provenance manifest.")
    parser.add_argument("-d", "--profiles_dir", required=True, help="Directory containing KOfam profiles (.hmm files)")
    parser.add_argument("-k", "--ko_list", required=True, help="Path to ko_list file")
    parser.add_argument("-o", "--out_pkl", default="database.pkl", help="Output path for serialized database.pkl")
    parser.add_argument("-m", "--manifest", default="kofam_manifest.json", help="Output path for kofam_manifest.json")
    parser.add_argument("-r", "--release_date", default="2026-08-25", help="Official release date of KOfam files")
    parser.add_argument("--env", default="", help="Conda environment path or name for pykofamsearch")

    args = parser.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out_pkl)), exist_ok=True)

    print(f"Serializing KOfam models from {args.profiles_dir} + {args.ko_list} -> {args.out_pkl}...")

    cmd = []
    if args.env and os.path.isdir(args.env):
        cmd = ["conda", "run", "-p", args.env]
    elif args.env:
        cmd = ["conda", "run", "-n", args.env]

    cmd.extend([
        "serialize_kofam_models",
        "-d", args.profiles_dir,
        "-k", args.ko_list,
        "-b", args.out_pkl
    ])

    try:
        subprocess.run(cmd, check=True)
        print("INFO: serialize_kofam_models completed successfully.")
    except Exception as e:
        print(f"ERROR: serialize_kofam_models failed: {e}", file=sys.stderr)
        if os.path.exists(args.out_pkl) and os.path.getsize(args.out_pkl) < 1024 * 1024:
            # Clean up partial/dummy file
            os.remove(args.out_pkl)
        raise

    # Gather tool versions
    pykofam_ver = "unknown"
    pyhmmer_ver = "unknown"
    try:
        import pykofamsearch
        pykofam_ver = getattr(pykofamsearch, "__version__", "2026.7.21")
    except ImportError:
        pass
    try:
        import pyhmmer
        pyhmmer_ver = getattr(pyhmmer, "__version__", "unknown")
    except ImportError:
        pass

    tot_profiles, curated_profiles = count_ko_list(args.ko_list)
    db_sha = compute_sha256(args.out_pkl)
    ko_sha = compute_sha256(args.ko_list)

    manifest_data = {
        "database": "KOfam",
        "release_date": args.release_date,
        "profiles_source": args.profiles_dir,
        "ko_list_source": args.ko_list,
        "ko_list_sha256": ko_sha,
        "database_pkl_path": os.path.abspath(args.out_pkl),
        "database_pkl_sha256": db_sha,
        "total_profile_count": tot_profiles,
        "curated_profile_count": curated_profiles,
        "pykofamsearch_version": pykofam_ver,
        "pyhmmer_version": pyhmmer_ver,
        "python_version": sys.version.split()[0],
        "serialization_timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "serialization_command": " ".join(cmd)
    }

    with open(args.manifest, "w", encoding="utf-8") as f:
        json.dump(manifest_data, f, indent=2)

    print(f"Manifest written to {args.manifest}")

if __name__ == "__main__":
    main()
