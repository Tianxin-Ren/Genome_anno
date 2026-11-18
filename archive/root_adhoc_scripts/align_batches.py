import os
import glob
import re
import yaml
import shutil

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

# Resolve paths
results_dir = os.path.join(project_root, "results/phylogeny")
batch_runs_dir = os.path.join(results_dir, "gtdbtk", "gses_all_mags_V0", "batch_runs")
batches_dir = os.path.join(results_dir, "gtdbtk", "gses_all_mags_V0", "batches")

print(f"MAG Directory: {mag_dir}")
mag_files = sorted(glob.glob(os.path.join(mag_dir, f"*.{ext}")))
print(f"Total MAG files on disk: {len(mag_files)}")

genome_to_path = {}
for f in mag_files:
    gid = os.path.splitext(os.path.basename(f))[0]
    genome_to_path[gid] = os.path.abspath(f)

# 1. Scan completed batch runs
completed_batch_folders = []
if os.path.isdir(batch_runs_dir):
    for bd in glob.glob(os.path.join(batch_runs_dir, "batch_*")):
        if os.path.isdir(bd):
            m = re.search(r'batch_(\d+)$', bd)
            if m:
                completed_batch_folders.append((int(m.group(1)), bd))

completed_batch_folders.sort()
print(f"Found {len(completed_batch_folders)} completed batch folders in batch_runs.")

# 2. Re-create batch_*.txt files with Deduplication
recreated_count = 0
total_lines_written = 0
active_bids = set()
seen_gids = set() # Store seen genomes to ensure uniqueness across all batches

for bid, bd in completed_batch_folders:
    raw_summaries = glob.glob(os.path.join(bd, "gtdbtk_out", "gtdbtk.*.summary.tsv"))
    genomes_in_summaries = []
    
    for sf in raw_summaries:
        if os.path.exists(sf) and os.path.getsize(sf) > 0:
            try:
                with open(sf, "r") as f:
                    header = f.readline()
                    for line in f:
                        line = line.strip()
                        if line:
                            parts = line.split("\t")
                            if parts:
                                gid = parts[0]
                                # Skip duplicate genomes we have already assigned to an earlier batch
                                if gid in seen_gids:
                                    continue
                                path = genome_to_path.get(gid)
                                if path:
                                    genomes_in_summaries.append((path, gid))
                                    seen_gids.add(gid)
            except Exception as e:
                print(f"Error reading summary {sf}: {e}")
                
    # Even if genomes_in_summaries is empty (all were duplicate/QC failed),
    # we still keep the batch.txt as empty or delete it.
    # To keep Snakemake completely happy and matched, we always write the file.
    active_bids.add(bid)
    batch_txt = os.path.join(batches_dir, f"batch_{bid}.txt")
    # Ensure parent dir exists
    os.makedirs(os.path.dirname(batch_txt), exist_ok=True)
    
    # Write clean genomes list (no duplicates)
    with open(batch_txt, "w") as out_f:
        for path, gid in genomes_in_summaries:
            out_f.write(f"{path}\t{gid}\n")
            total_lines_written += 1
    recreated_count += 1

print(f"Re-created {recreated_count} batch txt files under {batches_dir} with strict deduplication.")
print(f"Total unique genomes resolved and written: {len(seen_gids)}")
print(f"Total lines written across all batch txt files: {total_lines_written}")

# 3. Clean up/delete redundant batch_*.txt files that don't match completed runs
deleted_count = 0
for txt_path in glob.glob(os.path.join(batches_dir, "batch_*.txt")):
    filename = os.path.basename(txt_path)
    m = re.match(r'^batch_(\d+)\.txt$', filename)
    if m:
        bid = int(m.group(1))
        if bid not in active_bids:
            try:
                os.remove(txt_path)
                deleted_count += 1
            except Exception as e:
                print(f"Error deleting redundant file {filename}: {e}")

print(f"Deleted {deleted_count} redundant batch txt files from batches directory.")
print("\nDeduplication alignment finished successfully!")
