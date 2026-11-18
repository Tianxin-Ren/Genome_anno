import os
import glob
import re
import yaml

# Load config to get paths
conf_dir = "projects/config_phy"
samples_yaml = os.path.join(conf_dir, "samples.yaml")
phy_yaml = os.path.join(conf_dir, "phylogeny.yaml")

with open(samples_yaml) as f:
    samples_cfg = yaml.safe_load(f)
with open(phy_yaml) as f:
    phy_cfg = yaml.safe_load(f)

project_root = samples_cfg.get("project_root", "")
gtdbtk_tasks = samples_cfg.get("gtdbtk_tasks", {})
task_info = gtdbtk_tasks.get("gses_all_mags_V0", {})

mag_dir = task_info.get("mag_dir", "")
ext = task_info.get("extension", "fa")

# Resolve absolute path
if not os.path.isabs(mag_dir) and project_root:
    mag_dir = os.path.normpath(os.path.join(project_root, mag_dir))

print(f"MAG directory: {mag_dir}")
mag_files = sorted(glob.glob(os.path.join(mag_dir, f"*.{ext}")))
print(f"Total MAG files on disk: {len(mag_files)}")

genome_to_path = {}
for f in mag_files:
    gid = os.path.splitext(os.path.basename(f))[0]
    genome_to_path[gid] = os.path.abspath(f)

# Scan batch runs
results_dir = os.path.join(project_root, "results/phylogeny")
batch_runs_dir = os.path.join(results_dir, "gtdbtk", "gses_all_mags_V0", "batch_runs")
print(f"Batch runs directory: {batch_runs_dir}")

summary_files = glob.glob(os.path.join(batch_runs_dir, "batch_*", "gtdbtk.summary.tsv"))
print(f"Found {len(summary_files)} summary files")

matches = 0
mismatches = 0
sample_mismatches = []
sample_matches = []

for sf in sorted(summary_files, key=lambda x: int(re.search(r'batch_(\d+)', x).group(1))):
    bid = re.search(r'batch_(\d+)', sf).group(1)
    if os.path.exists(sf) and os.path.getsize(sf) > 0:
        file_matches = 0
        file_mismatches = 0
        with open(sf, "r") as f:
            header = f.readline()
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split("\t")
                    if parts:
                        gid = parts[0]
                        if gid in genome_to_path:
                            file_matches += 1
                            matches += 1
                            if len(sample_matches) < 5:
                                sample_matches.append(gid)
                        else:
                            file_mismatches += 1
                            mismatches += 1
                            if len(sample_mismatches) < 10:
                                sample_mismatches.append((bid, gid))
        # print(f"Batch {bid}: matches={file_matches}, mismatches={file_mismatches}")

print("\n=== DEBUG SUMMARY ===")
print(f"Total genome matches: {matches}")
print(f"Total genome mismatches: {mismatches}")
if sample_matches:
    print(f"Sample matches: {sample_matches}")
if sample_mismatches:
    print(f"Sample mismatches (Batch, GID): {sample_mismatches[:10]}")
