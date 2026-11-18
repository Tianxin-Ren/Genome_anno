import os
import glob
import re
import yaml

conf_dir = "projects/config_phy"
samples_yaml = os.path.join(conf_dir, "samples.yaml")

with open(samples_yaml) as f:
    samples_cfg = yaml.safe_load(f)

project_root = samples_cfg.get("project_root", "")
gtdbtk_tasks = samples_cfg.get("gtdbtk_tasks", {})
task_info = gtdbtk_tasks.get("gses_all_mags_V0", {})

mag_dir = task_info.get("mag_dir", "")
ext = task_info.get("extension", "fa")

if not os.path.isabs(mag_dir) and project_root:
    mag_dir = os.path.normpath(os.path.join(project_root, mag_dir))

mag_files = sorted(glob.glob(os.path.join(mag_dir, f"*.{ext}")))
total_mag_count = len(mag_files)

genome_to_path = {}
for f in mag_files:
    gid = os.path.splitext(os.path.basename(f))[0]
    genome_to_path[gid] = os.path.abspath(f)

results_dir = os.path.join(project_root, "results/phylogeny")
batch_runs_dir = os.path.join(results_dir, "gtdbtk", "gses_all_mags_V0", "batch_runs")

# Scan both bacteria and archaea summaries inside gtdbtk_out folders
raw_summaries = glob.glob(os.path.join(batch_runs_dir, "batch_*", "gtdbtk_out", "gtdbtk.*.summary.tsv"))

unique_completed_genomes = set()
total_rows_scanned = 0

for sf in raw_summaries:
    if os.path.exists(sf) and os.path.getsize(sf) > 0:
        with open(sf, "r") as f:
            header = f.readline()
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split("\t")
                    if parts:
                        gid = parts[0]
                        if gid in genome_to_path:
                            unique_completed_genomes.add(gid)
                            total_rows_scanned += 1

unique_completed_count = len(unique_completed_genomes)
remaining_count = total_mag_count - unique_completed_count

print("\n=== TRUE RUN COMPARISON (Bacteria + Archaea) ===")
print(f"Total MAG files on disk: {total_mag_count}")
print(f"Unique genomes successfully run (with raw results): {unique_completed_count}")
print(f"Remaining genomes to be run: {remaining_count}")
if unique_completed_count > 0:
    duplicate_count = total_rows_scanned - unique_completed_count
    print(f"Estimated duplicate runs: {duplicate_count}")
