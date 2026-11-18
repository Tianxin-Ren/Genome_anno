import os
import glob
import re
import shutil
import yaml
import sys

# Set DRY_RUN = True by default to prevent accidental modifications.
# Pass --run to execute the modifications.
DRY_RUN = True
if len(sys.argv) > 1 and sys.argv[1].lower() in ["--run", "run", "false"]:
    DRY_RUN = False

print("=============================================================")
print(f"GTDB-Tk Batch Cleanup Script (batches + batch_runs)")
print(f"Mode: {'DRY RUN (preview only)' if DRY_RUN else 'EXECUTE (making changes)'}")
if DRY_RUN:
    print("To execute changes, run: python3 cleanup_duplicate_batches.py --run")
print("=============================================================\n")

# Load config to get project_root
conf_dir = "projects/config_phy"
samples_yaml = os.path.join(conf_dir, "samples.yaml")

if not os.path.exists(samples_yaml):
    print(f"ERROR: samples config file not found at {samples_yaml}")
    exit(1)

with open(samples_yaml) as f:
    samples_cfg = yaml.safe_load(f)

project_root = samples_cfg.get("project_root", "")
results_dir = os.path.join(project_root, "results/phylogeny")
gtdbtk_dir = os.path.join(results_dir, "gtdbtk", "gses_all_mags_V0")

batch_runs_dir = os.path.join(gtdbtk_dir, "batch_runs")
batches_dir = os.path.join(gtdbtk_dir, "batches")

duplicates_runs_dir = os.path.join(gtdbtk_dir, "batch_runs_duplicates")
duplicates_batches_dir = os.path.join(gtdbtk_dir, "batches_duplicates")

print(f"Project Root: {project_root}")
print(f"Batch Runs Dir: {batch_runs_dir}")
print(f"Batches Dir: {batches_dir}")
print(f"Duplicates Runs Temp Destination: {duplicates_runs_dir}")
print(f"Duplicates Batches Temp Destination: {duplicates_batches_dir}\n")

if not os.path.isdir(batch_runs_dir) and not os.path.isdir(batches_dir):
    print("ERROR: Neither batch_runs nor batches directory exists.")
    exit(1)

if not DRY_RUN:
    os.makedirs(duplicates_runs_dir, exist_ok=True)
    os.makedirs(duplicates_batches_dir, exist_ok=True)

# 1. Scan both batches/ and batch_runs/ for all batch IDs
batch_ids_set = set()

# Scan batches/ for batch_*.txt
if os.path.isdir(batches_dir):
    for f in glob.glob(os.path.join(batches_dir, "batch_*.txt")):
        filename = os.path.basename(f)
        m = re.match(r'^batch_(\d+)\.txt$', filename)
        if m:
            batch_ids_set.add(int(m.group(1)))

# Scan batch_runs/ for batch_*
if os.path.isdir(batch_runs_dir):
    for p in glob.glob(os.path.join(batch_runs_dir, "batch_*")):
        if os.path.isdir(p):
            dirname = os.path.basename(p)
            m = re.match(r'^batch_(\d+)$', dirname)
            if m:
                batch_ids_set.add(int(m.group(1)))

sorted_batch_ids = sorted(list(batch_ids_set))
print(f"Discovered {len(sorted_batch_ids)} unique batch IDs from both folders.\n")

seen_genomes = set()
kept_batches = []       # list of (bid, genomes_set)
duplicate_batches = []  # list of (bid, genomes_set, reason)

# 2. Analyze each batch ID to determine if it is a duplicate
for bid in sorted_batch_ids:
    genomes_in_batch = set()
    source = ""
    
    # Option A: Read classified genomes from raw summaries in batch_runs/batch_{bid}/gtdbtk_out/
    bd = os.path.join(batch_runs_dir, f"batch_{bid}")
    raw_summaries = glob.glob(os.path.join(bd, "gtdbtk_out", "gtdbtk.*.summary.tsv"))
    
    for sf in raw_summaries:
        if os.path.exists(sf) and os.path.getsize(sf) > 0:
            with open(sf, "r") as f:
                header = f.readline()
                for line in f:
                    line = line.strip()
                    if line:
                        parts = line.split("\t")
                        if parts:
                            genomes_in_batch.add(parts[0])
            source = "gtdbtk_out summaries"
            
    # Option B: Read planned genomes from batches/batch_{bid}.txt (if no runs/summaries found)
    batch_txt = os.path.join(batches_dir, f"batch_{bid}.txt")
    if not genomes_in_batch and os.path.exists(batch_txt):
        with open(batch_txt, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        genomes_in_batch.add(parts[1])
                    elif len(parts) == 1:
                        # Fallback to file name without extension
                        gid = os.path.splitext(os.path.basename(parts[0]))[0]
                        genomes_in_batch.add(gid)
        source = "batch txt file"
        
    # Check if this batch is empty
    if not genomes_in_batch:
        print(f"Batch {bid}: empty (no raw summaries and no batch txt, or empty). Marked for removal.")
        duplicate_batches.append((bid, genomes_in_batch, "empty"))
        continue
        
    # Deduplication logic: a batch is kept if it has at least one genome not seen in previous kept batches
    new_genomes = genomes_in_batch - seen_genomes
    if not new_genomes:
        print(f"Batch {bid}: duplicate/fully redundant (all {len(genomes_in_batch)} genomes already covered). Marked for removal.")
        duplicate_batches.append((bid, genomes_in_batch, "duplicate"))
    else:
        print(f"Batch {bid}: kept (source: {source}, contains {len(genomes_in_batch)} genomes, {len(new_genomes)} are new).")
        kept_batches.append((bid, genomes_in_batch))
        seen_genomes.update(genomes_in_batch)

print("\n=============================================")
print("Analysis Summary:")
print(f"  Total Batch IDs: {len(sorted_batch_ids)}")
print(f"  Batches to Keep: {len(kept_batches)}")
print(f"  Batches to Remove: {len(duplicate_batches)}")
print(f"  Total Unique Genomes Covered: {len(seen_genomes)}")
print("=============================================\n")

# 3. Move duplicate batches (both folders and txt files)
print("Moving duplicate batches to backup directories:")
for bid, _, reason in duplicate_batches:
    # A. Move run folder
    bd = os.path.join(batch_runs_dir, f"batch_{bid}")
    if os.path.exists(bd):
        dst_run = os.path.join(duplicates_runs_dir, f"batch_{bid}")
        print(f"  [MOVE] Run folder: {bd} -> {dst_run}")
        if not DRY_RUN:
            if os.path.exists(dst_run):
                if os.path.isdir(dst_run):
                    shutil.rmtree(dst_run)
                else:
                    os.remove(dst_run)
            shutil.move(bd, dst_run)
            
    # B. Move batch file
    batch_txt = os.path.join(batches_dir, f"batch_{bid}.txt")
    if os.path.exists(batch_txt):
        dst_batch = os.path.join(duplicates_batches_dir, f"batch_{bid}.txt")
        print(f"  [MOVE] Batch file: {batch_txt} -> {dst_batch}")
        if not DRY_RUN:
            if os.path.exists(dst_batch):
                os.remove(dst_batch)
            shutil.move(batch_txt, dst_batch)

# 4. Rename kept batches (Step 1: Temp prefix to prevent name collisions)
print("\nRenaming kept batches (Step 1: Move to temporary name):")
for i, (bid, _) in enumerate(kept_batches):
    target_bid = i + 1
    
    temp_run_name = f"temp_rename_batch_{target_bid}"
    temp_batch_name = f"temp_rename_batch_{target_bid}.txt"
    
    # A. Rename run folder
    bd = os.path.join(batch_runs_dir, f"batch_{bid}")
    if os.path.exists(bd):
        dst_temp_run = os.path.join(batch_runs_dir, temp_run_name)
        print(f"  [RENAME-1] Run folder: {bd} -> {dst_temp_run}")
        if not DRY_RUN:
            if os.path.exists(dst_temp_run):
                if os.path.isdir(dst_temp_run):
                    shutil.rmtree(dst_temp_run)
                else:
                    os.remove(dst_temp_run)
            shutil.move(bd, dst_temp_run)
            
    # B. Rename batch file
    batch_txt = os.path.join(batches_dir, f"batch_{bid}.txt")
    if os.path.exists(batch_txt):
        dst_temp_batch = os.path.join(batches_dir, temp_batch_name)
        print(f"  [RENAME-1] Batch file: {batch_txt} -> {dst_temp_batch}")
        if not DRY_RUN:
            if os.path.exists(dst_temp_batch):
                os.remove(dst_temp_batch)
            shutil.move(batch_txt, dst_temp_batch)

# 5. Rename kept batches (Step 2: Temp prefix to final sequential name)
print("\nRenaming kept batches (Step 2: Move to final sequential name):")
for i, (bid, _) in enumerate(kept_batches):
    target_bid = i + 1
    
    temp_run_name = f"temp_rename_batch_{target_bid}"
    temp_batch_name = f"temp_rename_batch_{target_bid}.txt"
    
    final_run_name = f"batch_{target_bid}"
    final_batch_name = f"batch_{target_bid}.txt"
    
    # A. Rename run folder from temp to final
    temp_run_path = os.path.join(batch_runs_dir, temp_run_name)
    final_run_path = os.path.join(batch_runs_dir, final_run_name)
    
    if DRY_RUN:
        bd = os.path.join(batch_runs_dir, f"batch_{bid}")
        if os.path.exists(bd):
            print(f"  [RENAME-2] Run folder: {bd} (temp) -> {final_run_path}")
            print(f"    [REMOVE] Force rebuild: {os.path.join(final_run_path, 'gtdbtk.summary.tsv')}")
    else:
        if os.path.exists(temp_run_path):
            print(f"  [RENAME-2] Run folder: {temp_run_path} -> {final_run_path}")
            if os.path.exists(final_run_path):
                if os.path.isdir(final_run_path):
                    shutil.rmtree(final_run_path)
                else:
                    os.remove(final_run_path)
            shutil.move(temp_run_path, final_run_path)
            
            # Delete the old merged summary file to force Snakemake to rebuild it with BOTH bacteria and archaea
            summary_file = os.path.join(final_run_path, "gtdbtk.summary.tsv")
            if os.path.exists(summary_file):
                print(f"    [REMOVE] Force rebuild: {summary_file}")
                os.remove(summary_file)
                
    # B. Rename batch file from temp to final
    temp_batch_path = os.path.join(batches_dir, temp_batch_name)
    final_batch_path = os.path.join(batches_dir, final_batch_name)
    
    if DRY_RUN:
        batch_txt = os.path.join(batches_dir, f"batch_{bid}.txt")
        if os.path.exists(batch_txt):
            print(f"  [RENAME-2] Batch file: {batch_txt} (temp) -> {final_batch_path}")
    else:
        if os.path.exists(temp_batch_path):
            print(f"  [RENAME-2] Batch file: {temp_batch_path} -> {final_batch_path}")
            if os.path.exists(final_batch_path):
                os.remove(final_batch_path)
            shutil.move(temp_batch_path, final_batch_path)

print(f"\nBatch cleanup {'preview (dry-run)' if DRY_RUN else 'execution'} finished successfully!")
