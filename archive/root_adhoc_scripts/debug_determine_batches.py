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

def py_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    return os.path.normpath(os.path.join(project_root, rel_path))

results_dir = os.path.join(project_root, "results/phylogeny")

print(f"MAG directory: {mag_dir}")
mag_files = sorted(glob.glob(os.path.join(mag_dir, f"*.{ext}")))
print(f"Total MAG files on disk: {len(mag_files)}")

genome_to_path = {}
for f in mag_files:
    gid = os.path.splitext(os.path.basename(f))[0]
    genome_to_path[gid] = os.path.abspath(f)

batch_runs_dir = os.path.join(results_dir, "gtdbtk", "gses_all_mags_V0", "batch_runs")
batches_dir = os.path.join(results_dir, "gtdbtk", "gses_all_mags_V0", "batches")

classified_paths = set()
completed_batches = {}
completed_batch_ids = set()

# Scan completed batches
if os.path.isdir(batch_runs_dir):
    batch_folders = glob.glob(os.path.join(batch_runs_dir, "batch_*"))
    print(f"Found {len(batch_folders)} batch directories on disk.")
    
    for bd in sorted(batch_folders, key=lambda x: int(re.search(r'batch_(\d+)', x).group(1)) if re.search(r'batch_(\d+)', x) else 999999):
        if not os.path.isdir(bd):
            continue
        m = re.search(r'batch_(\d+)$', bd)
        if m:
            bid = int(m.group(1))
            bid_str = str(bid)
            
            raw_summaries = glob.glob(os.path.join(bd, "gtdbtk_out", "gtdbtk.*.summary.tsv"))
            is_valid = raw_summaries and any(os.path.exists(sf) and os.path.getsize(sf) > 0 for sf in raw_summaries)
            
            if is_valid:
                paths_in_batch = completed_batches.setdefault(bid_str, [])
                completed_batch_ids.add(bid)
                
                batch_txt = os.path.join(batches_dir, f"batch_{bid}.txt")
                if os.path.exists(batch_txt) and os.path.getsize(batch_txt) > 0:
                    try:
                        with open(batch_txt, "r") as f:
                            for line in f:
                                line = line.strip()
                                if line:
                                    parts = line.split("\t")
                                    if len(parts) >= 2:
                                        gid = parts[1]
                                        path = genome_to_path.get(gid)
                                        if path and path not in paths_in_batch:
                                            paths_in_batch.append(path)
                                            classified_paths.add(path)
                    except Exception as e:
                        print(f"Error reading {batch_txt}: {e}")
                
                # Always read raw summaries as well to handle out-of-sync batch_txt files
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
                                            path = genome_to_path.get(gid)
                                            if path and path not in paths_in_batch:
                                                paths_in_batch.append(path)
                                                classified_paths.add(path)
                        except Exception as e:
                            print(f"Error reading summary {sf}: {e}")
                
                if not paths_in_batch:
                    completed_batches.pop(bid_str, None)
                    completed_batch_ids.discard(bid)
            else:
                print(f"Batch {bid} is NOT valid (no non-empty raw summaries found in gtdbtk_out)")

uncompleted_files = []
for f in mag_files:
    abs_f = os.path.abspath(f)
    if abs_f not in classified_paths:
        uncompleted_files.append(f)

print(f"\nCompleted batch IDs: {sorted(list(completed_batch_ids))}")
print(f"Total classified paths resolved: {len(classified_paths)}")
print(f"Total uncompleted files: {len(uncompleted_files)}")

if uncompleted_files:
    print("\nFirst 10 uncompleted files:")
    for f in uncompleted_files[:10]:
        gid = os.path.splitext(os.path.basename(f))[0]
        # Check if this gid is in genome_to_path
        in_map = gid in genome_to_path
        print(f"  - Path: {f} (GID: {gid}, in_map: {in_map})")
        
        # Check if this gid appears in any batch_*.txt
        found_in_txt = []
        for txt in glob.glob(os.path.join(batches_dir, "batch_*.txt")):
            try:
                with open(txt, "r") as tf:
                    for line in tf:
                        if gid in line:
                            found_in_txt.append(os.path.basename(txt))
            except Exception:
                pass
        print(f"    Found in batch text files: {found_in_txt}")
