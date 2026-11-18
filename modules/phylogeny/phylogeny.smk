# ==============================================================================
# Phylogeny Module: GToTree and IQ-TREE Workflow (Multi-task Edition)
# ==============================================================================
import os
import glob
import re

localrules: gtdbtk_prepare_batch_files, merge_gtdbtk_batch_summary

# 1. Global Path Logic (Namespaced)
# ------------------------------------------------------------------------------
# Resolve project_root globally. All relative paths will be prefixed with this.
PY_ROOT = config.get("project_root", "")
PY_ENV_DIR = config.get("conda_env_dir", "")
PY_C = config.get("phylogeny", {})

# Load the dictionary of phylogeny tasks from samples.yaml
PHY_TASKS = config.get("phylogeny_tasks", {})

# Load the dictionary of GTDB-Tk tasks from samples.yaml
GTDBTK_TASKS = config.get("gtdbtk_tasks", {})
GTDBTK_C = PY_C.get("gtdbtk", {})
GTDBTK_ENABLE = GTDBTK_C.get("enable", False)

# Load the dictionary of subtyping tasks from samples.yaml
SUBTYPING_TASKS = config.get("subtyping_tasks", {})
SUBTYPING_C = PY_C.get("subtyping", {})
SUBTYPING_ENABLE = SUBTYPING_C.get("enable", False)

def py_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    # Prefix relative paths with the global project_root
    return os.path.normpath(os.path.join(PY_ROOT, rel_path))

def py_get_env(name):
    env_name = config.get("environments", {}).get(name)
    if not env_name and name == "gtdbtk":
        env_name = GTDBTK_C.get("conda_envs", {}).get("gtdbtk")
    if not env_name:
        env_name = PY_C.get("conda_envs", {}).get(name, name)
    if not PY_ENV_DIR or os.path.isabs(env_name): return env_name
    return os.path.join(PY_ENV_DIR, env_name)


# Results directories
PY_RES = py_full_path(PY_C.get("results_dir", "results/phylogeny"))
PY_LOGS = py_full_path(PY_C.get("logs_dir", os.path.join("logs", "phylogeny")))

# Dynamic batch calculation and incremental loading function for GTDB-Tk
def determine_gtdbtk_batches(sample_name, task_info, global_c, results_dir):
    mag_dir = py_full_path(task_info.get("mag_dir", ""))
    ext = task_info.get("extension", "fa")
    batch_size = int(global_c.get("batch_size", 50))
    
    if not os.path.isdir(mag_dir):
        return {}, []
        
    # Find all MAG files
    mag_files = sorted(glob.glob(os.path.join(mag_dir, f"*.{ext}")))
    if not mag_files:
        return {}, []
        
    # Build mapping from genome ID to path
    genome_to_path = {}
    for f in mag_files:
        gid = os.path.splitext(os.path.basename(f))[0]
        genome_to_path[gid] = os.path.abspath(f)
        
    batch_runs_dir = os.path.join(results_dir, "gtdbtk", sample_name, "batch_runs")
    batches_dir = os.path.join(results_dir, "gtdbtk", sample_name, "batches")
    
    classified_paths = set()
    completed_batches = {}
    completed_batch_ids = set()
    
    # 1. Scan completed batches
    if os.path.isdir(batch_runs_dir):
        # We look for folders matching batch_*
        batch_folders = glob.glob(os.path.join(batch_runs_dir, "batch_*"))
        for bd in batch_folders:
            if not os.path.isdir(bd):
                continue
            m = re.search(r'batch_(\d+)$', bd)
            if m:
                bid = int(m.group(1))
                bid_str = str(bid)
                
                # A batch is completed if raw summaries exist inside gtdbtk_out and are not empty
                raw_summaries = glob.glob(os.path.join(bd, "gtdbtk_out", "gtdbtk.*.summary.tsv"))
                if raw_summaries and any(os.path.exists(sf) and os.path.getsize(sf) > 0 for sf in raw_summaries):
                    paths_in_batch = completed_batches.setdefault(bid_str, [])
                    completed_batch_ids.add(bid)
                    
                    # 1. Try to read the original batch definition file (to get failures)
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
                        except Exception:
                            pass
                    
                    # 2. Always also read raw summaries (to handle cases where batch_txt is out of sync or missing)
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
                            except Exception:
                                pass
                                    
                    # If somehow no paths could be resolved for this batch, clean it up from completed dict
                    if not paths_in_batch:
                        completed_batches.pop(bid_str, None)
                        completed_batch_ids.discard(bid)

    # 2. Identify remaining (uncompleted) MAG files
    uncompleted_files = []
    for f in mag_files:
        abs_f = os.path.abspath(f)
        if abs_f not in classified_paths:
            uncompleted_files.append(f)
            
    # 3. Group uncompleted MAG files into chunks and assign unused batch IDs
    new_batches = {}
    if uncompleted_files:
        next_bid = 1
        for i in range(0, len(uncompleted_files), batch_size):
            # Find the next unused batch ID
            while next_bid in completed_batch_ids:
                next_bid += 1
            
            bid_str = str(next_bid)
            new_batches[bid_str] = uncompleted_files[i : i + batch_size]
            next_bid += 1
            
    # Combine completed and new batches
    all_batches = {}
    all_batches.update(completed_batches)
    all_batches.update(new_batches)
    
    batch_ids = sorted(list(all_batches.keys()), key=int)
    return all_batches, batch_ids

# Pre-calculate batches for all active GTDB-Tk tasks
GTDBTK_BATCH_DATA = {} # Maps task_name -> (batches_dict, batch_ids_list)
if GTDBTK_ENABLE and GTDBTK_TASKS:
    for task_name, task_info in GTDBTK_TASKS.items():
        batches, batch_ids = determine_gtdbtk_batches(
            task_name, task_info, GTDBTK_C, PY_RES
        )
        GTDBTK_BATCH_DATA[task_name] = (batches, batch_ids)


def py_get_genome_id(path):
    """Extract GToTree genome identifier (filename stripped of compression and fasta extension)."""
    base = os.path.basename(path)
    if base.endswith(".gz"):
        base = base[:-3]
    base = os.path.splitext(base)[0]
    return base

# Helper function to discover input genomes for a specific task configuration
def find_genomes_for_task(task_info):
    genome_paths = []
    
    # 1. Individual files or glob patterns
    files = task_info.get("genome_files", [])
    if files:
        if isinstance(files, str):
            files = [files]
        for f in files:
            p = py_full_path(f)
            if os.path.isfile(p):
                genome_paths.append(p)
            else:
                for gp in glob.glob(p):
                    if os.path.isfile(gp):
                        genome_paths.append(gp)
                        
    # 2. Directories to scan for fasta files
    dirs = task_info.get("genomes_dirs", [])
    if dirs:
        if isinstance(dirs, str):
            dirs = [dirs]
        for d in dirs:
            p = py_full_path(d)
            if os.path.isdir(p):
                # Scan common fasta extensions (and gzipped variants)
                for ext in ["*.fna", "*.fasta", "*.fa", "*.fna.gz", "*.fasta.gz", "*.fa.gz"]:
                    for f_path in glob.glob(os.path.join(p, ext)):
                        if os.path.isfile(f_path):
                            genome_paths.append(f_path)
                            
    # 3. Path list file (e.g. gtotree_genome_path.tsv)
    path_list_file = task_info.get("genome_path_tsv") or task_info.get("genome_path_list")
    if path_list_file:
        p = py_full_path(path_list_file)
        if os.path.isfile(p):
            with open(p, "r") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        # Extract first field if it's a TSV
                        gp = py_full_path(line.split("\t")[0])
                        if os.path.isfile(gp):
                            genome_paths.append(gp)
                        elif glob.glob(gp):
                            for matches in glob.glob(gp):
                                if os.path.isfile(matches):
                                    genome_paths.append(matches)
                            
    # Deduplicate while preserving order:
    # 1. By absolute file path
    # 2. By GToTree genome identifier (basename without fasta/gz extension)
    # GToTree requires every genome to have a distinct base identifier.
    # If two files have the same identifier (e.g. u1202559.fa in two dirs, or u1202559.fa and u1202559.fa.gz),
    # GToTree's gtt-reorder-fasta crashes with BioPython "ValueError: Duplicate key '<id>'".
    seen_paths = set()
    seen_ids = set()
    dedup_paths = []
    for gp in genome_paths:
        abs_p = os.path.abspath(gp)
        if abs_p in seen_paths:
            continue
        seen_paths.add(abs_p)
        
        gid = py_get_genome_id(abs_p)
        if gid in seen_ids:
            print(f"[WARNING] Duplicate genome identifier '{gid}' detected for '{abs_p}'. Skipping duplicate to prevent GToTree crash.")
            continue
        seen_ids.add(gid)
        dedup_paths.append(abs_p)
            
    return dedup_paths


# 2. Target Registration
# ------------------------------------------------------------------------------
PY_TARGETS = []

if PY_C.get("enable", False) and PHY_TASKS:
    # Iterate through all configured phylogeny tasks and register their final output .treefiles
    # Output structure: results/phylogeny/iqtree/{sample}/{sample}.iqtree.treefile
    for task_name in PHY_TASKS:
        PY_TARGETS.append(os.path.join(PY_RES, f"iqtree/{task_name}/{task_name}.iqtree.treefile"))

# Add GTDB-Tk targets if enabled and active
if GTDBTK_ENABLE and GTDBTK_TASKS:
    for task_name in GTDBTK_TASKS:
        # Register the final merged classification summary for each task
        PY_TARGETS.append(os.path.join(PY_RES, f"gtdbtk/{task_name}/merged_gtdbtk_summary.tsv"))

# Add subtyping targets if enabled and active
if SUBTYPING_ENABLE and SUBTYPING_TASKS:
    for task_name in SUBTYPING_TASKS:
        # Register the final per_query.tsv classification output
        PY_TARGETS.append(os.path.join(PY_RES, f"subtyping/{task_name}/subtype_classification_per_query.tsv"))

TARGETS.extend(PY_TARGETS)



# 3. Rules
# ------------------------------------------------------------------------------

rule phylogeny_prepare_genomes:
    output:
        # CRITICAL FIX: Write the paths file OUTSIDE of GToTree's output directory
        # GToTree's -F flag deletes the entire output directory (-o) upon starting, 
        # so keeping the input list inside the output folder causes GToTree to delete its own input.
        genome_list = os.path.join(PY_RES, "gtotree/genome_lists/{sample}_paths.txt")
    log:
        os.path.join(PY_LOGS, "gtotree/prepare_genomes_{sample}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    run:
        task_name = wildcards.sample
        task_info = PHY_TASKS.get(task_name)
        if not task_info:
            raise ValueError(f"ERROR: Phylogeny task '{task_name}' is not defined in samples.yaml")
        
        with open(log[0], "w") as log_f:
            log_f.write(f"=== PHYLOGENY GENOME INPUT SCAN FOR TASK: {task_name} ===\n")
            log_f.write(f"Global Project Root: {PY_ROOT}\n")
            log_f.write(f"HMM Profile: {task_info.get('hmm_profile', 'Bacteria')}\n")
            log_f.write("Scanning genomes configured under task...\n")
            
            paths = find_genomes_for_task(task_info)
            if not paths:
                msg = (
                    f"ERROR: No valid genomes found for phylogeny task '{task_name}'! "
                    "Please configure one of the following under 'phylogeny_tasks' in samples.yaml:\n"
                    "  - 'genome_path_tsv' with a valid list file\n"
                    "  - 'genomes_dirs' with directories containing .fa/.fasta/.fna files\n"
                    "  - 'genome_files' with direct file paths\n"
                )
                log_f.write(msg)
                raise ValueError(msg)
            
            log_f.write(f"Successfully discovered {len(paths)} genomes:\n")
            for p in paths:
                log_f.write(f"  - {p}\n")
            
            # Ensure the output directory exists
            os.makedirs(os.path.dirname(output.genome_list), exist_ok=True)
            with open(output.genome_list, "w") as out_f:
                for p in paths:
                    out_f.write(f"{p}\n")
            log_f.write(f"Genome list written to: {output.genome_list}\n")


rule run_gtotree:
    input:
        genome_list = os.path.join(PY_RES, "gtotree/genome_lists/{sample}_paths.txt")
    output:
        msa = os.path.join(PY_RES, "gtotree/{sample}/Aligned_SCGs.faa")
    log:
        os.path.join(PY_LOGS, "gtotree/run_{sample}.log")
    threads: get_threads(64)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = py_get_env("gtotree"),
        hmm_profile = lambda w: PHY_TASKS.get(w.sample, {}).get("hmm_profile", "Bacteria"),
        out_dir = lambda w: os.path.join(PY_RES, "gtotree", w.sample)
    shell:
        """
        set -euo pipefail
        
        # Ensure the parent directory of GToTree's output directory exists
        mkdir -p "$(dirname "{params.out_dir}")"
        
        # Execute GToTree using designated Conda environment.
        # -F flag will automatically force overwrite and clean up the {params.out_dir} safely,
        # without affecting the input {input.genome_list} which is safely isolated outside.
        conda run -p {params.env} GToTree \
            -j {threads} \
            -f {input.genome_list} \
            -H {params.hmm_profile} \
            -o {params.out_dir} \
            -F 2>&1 | tee {log}
        """


rule run_iqtree:
    input:
        msa = os.path.join(PY_RES, "gtotree/{sample}/Aligned_SCGs.faa")
    output:
        treefile = os.path.join(PY_RES, "iqtree/{sample}/{sample}.iqtree.treefile")
    log:
        os.path.join(PY_LOGS, "iqtree/{sample}/{sample}.log")
    threads: get_threads(64)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = py_get_env("iqtree"),
        gtotree_out_dir = lambda w: os.path.join(PY_RES, "gtotree", w.sample),
        out_prefix = lambda w: os.path.join(PY_RES, "iqtree", w.sample, f"{w.sample}.iqtree"),
        mem_g = lambda wildcards, threads: max(1, int(threads * 8)),
        executable = PY_C.get("iqtree_executable", "iqtree3"),
        bootstrap = PY_C.get("iqtree_bootstrap", "-B 1000"),  # -B 1000 is Ultrafast Bootstrap in IQ-TREE 2.x/3.x
        alrt = PY_C.get("iqtree_alrt", "-alrt 1000"),
        model = PY_C.get("iqtree_model", "MFP"),
        extra = PY_C.get("iqtree_extra", "")
    shell:
        """
        set -euo pipefail
        
        mkdir -p "$(dirname "{output.treefile}")"
        
        # 1. GToTree v1.8.12 stores its partition file at run_files/Partitions.txt
        PART_FILE=""
        if [ -f "{params.gtotree_out_dir}/run_files/Partitions.txt" ]; then
            PART_FILE="-p {params.gtotree_out_dir}/run_files/Partitions.txt"
            echo "INFO: Enabling partition model with {params.gtotree_out_dir}/run_files/Partitions.txt"
        fi
        
        # 2. Defensive check for partition-resampling parameters
        # If no partition file was detected/used, passing partition-based options
        # like --sampling GENESITE will cause IQ-TREE to crash.
        EXTRA_ARGS="{params.extra}"
        if [ -z "$PART_FILE" ]; then
            echo "WARNING: run_files/Partitions.txt not found. Stripping partition-specific sampling arguments to prevent IQ-TREE crash."
            # Strip --sampling <value> or -bsam <value> using robust regex patterns
            EXTRA_ARGS=$(echo "$EXTRA_ARGS" | sed -E 's/--sampling[ =][^ ]+//g' | sed -E 's/-bsam[ =][^ ]+//g')
        fi
        
        # 3. Dynamic Memory Option Integration
        # IQ-TREE 3.1.2 does not support the --mem option when running partition models (-p).
        # We dynamically omit the memory cap parameter if a partition file is active.
        MEM_ARG=""
        if [ -z "$PART_FILE" ]; then
            MEM_ARG="--mem {params.mem_g}G"
        else
            echo "INFO: Partition model is active. Omiting --mem parameter for IQ-TREE 3.1.2 compatibility."
        fi
        
        # 4. Handle checkpoint mismatch during transition to partition models
        # If transitioning from a non-partition model to a partition model, 
        # IQ-TREE's previous ModelFinder checkpoint file (.model.gz) will be incompatible and cause a crash.
        # We automatically detect this mismatch and delete the stale .model.gz file to allow fresh partition ModelFinder.
        if [ -n "$PART_FILE" ] && [ -f "{params.out_prefix}.model.gz" ]; then
            echo "INFO: Transitioning to partition model. Deleting incompatible .model.gz checkpoint file to prevent crash."
            rm -f "{params.out_prefix}.model.gz"
        fi
        
        # Execute IQ-TREE using designated Conda environment
        conda run -p {params.env} {params.executable} \
            $MEM_ARG \
            -T {threads} \
            -m {params.model} \
            -s {input.msa} \
            $PART_FILE \
            {params.alrt} \
            {params.bootstrap} \
            $EXTRA_ARGS \
            --prefix {params.out_prefix} 2>&1 | tee {log}
        """


# ==============================================================================
# GTDB-Tk Batch Classification and Merging Rules
# ==============================================================================

# Helper to get all expected batch summaries for a given sample
def get_gtdbtk_batch_summaries(wildcards):
    task_name = wildcards.sample
    _, batch_ids = GTDBTK_BATCH_DATA.get(task_name, ({}, []))
    return [os.path.join(PY_RES, f"gtdbtk/{task_name}/batch_runs/batch_{bid}/gtdbtk.summary.tsv") for bid in batch_ids]

rule gtdbtk_prepare_batch_files:
    output:
        batch_file = os.path.join(PY_RES, "gtdbtk/{sample}/batches/batch_{batch_id}.txt")
    log:
        os.path.join(PY_LOGS, "gtdbtk/prepare_{sample}_batch_{batch_id}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        mag_paths = lambda wildcards: GTDBTK_BATCH_DATA.get(wildcards.sample, ({}, []))[0].get(wildcards.batch_id, [])
    run:
        task_name = wildcards.sample
        batch_id = wildcards.batch_id
        mag_paths = params.mag_paths
        
        if not mag_paths:
            raise ValueError(f"ERROR: Batch ID '{batch_id}' for task '{task_name}' is empty or not found in params.")

            
        with open(log[0], "w") as log_f:
            log_f.write(f"Preparing batch file for task: {task_name}, batch: {batch_id}\n")
            log_f.write(f"Contains {len(mag_paths)} MAGs\n")
            
            # Ensure parent directory exists
            os.makedirs(os.path.dirname(output.batch_file), exist_ok=True)
            
            # Generate the proposed content
            lines_to_write = []
            for path in mag_paths:
                genome_id = os.path.splitext(os.path.basename(path))[0]
                lines_to_write.append(f"{path}\t{genome_id}\n")
                log_f.write(f"  - {genome_id} -> {path}\n")
            content = "".join(lines_to_write)
            
            # Write only if content differs to preserve timestamp and prevent triggering Snakemake downstreams
            write_needed = True
            if os.path.exists(output.batch_file):
                try:
                    with open(output.batch_file, "r") as f:
                        if f.read() == content:
                            write_needed = False
                            log_f.write("Batch file content is identical to existing file. Skipping write to preserve timestamp.\n")
                except Exception:
                    pass
            
            if write_needed:
                with open(output.batch_file, "w") as out_f:
                    out_f.write(content)

rule run_gtdbtk_batch:
    input:
        batch_file = os.path.join(PY_RES, "gtdbtk/{sample}/batches/batch_{batch_id}.txt")
    output:
        outdir = directory(os.path.join(PY_RES, "gtdbtk/{sample}/batch_runs/batch_{batch_id}/gtdbtk_out"))
    log:
        os.path.join(PY_LOGS, "gtdbtk/{sample}_batch_{batch_id}.log")
    threads: get_threads(64)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = py_get_env("gtdbtk"),
        db = GTDBTK_C.get("db", "")
    shell:
        """
        set -euo pipefail
        
        # Defensive check: if raw summary files already exist in output directory, do not run GTDB-Tk again.
        if ls "{output.outdir}"/gtdbtk.*.summary.tsv >/dev/null 2>&1; then
            # Verify that at least one is non-empty
            non_empty=0
            for sf in "{output.outdir}"/gtdbtk.*.summary.tsv; do
                if [ -s "$sf" ]; then
                    non_empty=1
                    break
                fi
            done
            if [ "$non_empty" -eq 1 ]; then
                echo "INFO: GTDB-Tk raw summary files already exist in {output.outdir}. Skipping classify_wf execution." | tee -a "{log}"
                exit 0
            fi
        fi
        
        mkdir -p "{output.outdir}"
        
        if [ -n "{params.db}" ]; then
            export GTDBTK_DATA_PATH="{params.db}"
        fi
        
        echo "=== RUNNING GTDB-TK BATCH {wildcards.batch_id} ===" > "{log}"
        echo "Conda Env: {params.env}" >> "{log}"
        echo "Threads: {threads}" >> "{log}"
        
        # Run GTDB-Tk classify_wf using batch file
        conda run -p {params.env} gtdbtk classify_wf \
            --batchfile "{input.batch_file}" \
            --out_dir "{output.outdir}" \
            --cpus {threads} 2>&1 | tee -a "{log}"
            
        echo "=== BATCH {wildcards.batch_id} FINISHED ===" >> "{log}"
        """

rule merge_gtdbtk_batch_summary:
    input:
        outdir = os.path.join(PY_RES, "gtdbtk/{sample}/batch_runs/batch_{batch_id}/gtdbtk_out")
    output:
        summary = os.path.join(PY_RES, "gtdbtk/{sample}/batch_runs/batch_{batch_id}/gtdbtk.summary.tsv")
    log:
        os.path.join(PY_LOGS, "gtdbtk/{sample}_merge_batch_{batch_id}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        """
        set -euo pipefail
        
        summary_out="{output.summary}"
        
        # Dynamically discover summary files to support different database versions (ar122 vs ar53, bac120)
        bac_sum=$(ls "{input.outdir}"/gtdbtk.bac*.summary.tsv 2>/dev/null | head -n 1 || true)
        ar_sum=$(ls "{input.outdir}"/gtdbtk.ar*.summary.tsv 2>/dev/null | head -n 1 || true)
        
        if [ -n "$bac_sum" ] && [ -f "$bac_sum" ] && [ -n "$ar_sum" ] && [ -f "$ar_sum" ]; then
            echo "INFO: Found both bacteria ($bac_sum) and archaea ($ar_sum) summaries. Merging..." > "{log}"
            head -n 1 "$bac_sum" > "$summary_out"
            tail -n +2 "$bac_sum" >> "$summary_out"
            tail -n +2 "$ar_sum" >> "$summary_out"
        elif [ -n "$bac_sum" ] && [ -f "$bac_sum" ]; then
            echo "INFO: Found bacteria summary ($bac_sum) only. Copying..." > "{log}"
            cp "$bac_sum" "$summary_out"
        elif [ -n "$ar_sum" ] && [ -f "$ar_sum" ]; then
            echo "INFO: Found archaea summary ($ar_sum) only. Copying..." > "{log}"
            cp "$ar_sum" "$summary_out"
        else
            echo "WARNING: No GTDB-Tk summary files found for batch {wildcards.batch_id}." > "{log}"
            # Write a dummy header to avoid breaking the merge step
            echo -e "user_genome\tclassification\tfastani_reference\tfastani_ani\tfastani_align_fraction\tfastani_reference_radius\tfastani_taxonomy\tred_value\tred_taxonomy\tpplacer_taxonomy\tclassification_method\tnote\tother_significantly_alike_genomes\tred_values_other_significantly_alike_genomes\twarnings" > "$summary_out"
        fi
        """

rule gtdbtk_merge_results:
    input:
        get_gtdbtk_batch_summaries
    output:
        merged = os.path.join(PY_RES, "gtdbtk/{sample}/merged_gtdbtk_summary.tsv")
    log:
        os.path.join(PY_LOGS, "gtdbtk/{sample}_merge.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    run:
        with open(log[0], "w") as log_f:
            log_f.write(f"Merging {len(input)} batch summaries for task {wildcards.sample}...\n")
            
            # Ensure parent directory of merged output exists
            os.makedirs(os.path.dirname(output.merged), exist_ok=True)
            
            seen_genomes = set()
            header_written = False
            with open(output.merged, "w") as out_f:
                for path in input:
                    log_f.write(f"  Reading: {path}\n")
                    if not os.path.exists(path) or os.path.getsize(path) == 0:
                        log_f.write(f"    Warning: Empty or missing summary {path}\n")
                        continue
                    
                    with open(path, "r") as in_f:
                        lines = in_f.readlines()
                        if not lines:
                            continue
                        
                        header = lines[0]
                        data_lines = lines[1:]
                        
                        if not header_written:
                            out_f.write(header)
                            header_written = True
                        
                        for dl in data_lines:
                            parts = dl.split("\t")
                            if parts:
                                gid = parts[0]
                                if gid not in seen_genomes:
                                    out_f.write(dl)
                                    seen_genomes.add(gid)
                                else:
                                    log_f.write(f"    Skipping duplicate genome: {gid}\n")
            
            log_f.write(f"Merged output written to: {output.merged}\n")


rule run_protein_subtyping:
    input:
        ref_fasta = lambda w: py_full_path(SUBTYPING_TASKS[w.sample].get("ref_fasta", "")),
        query_fasta = lambda w: py_full_path(SUBTYPING_TASKS[w.sample].get("query_fasta", ""))
    output:
        per_query = os.path.join(PY_RES, "subtyping/{sample}/subtype_classification_per_query.tsv"),
        profile = os.path.join(PY_RES, "subtyping/{sample}/subtype_classification_profile.tsv")
    log:
        os.path.join(PY_LOGS, "subtyping/{sample}.log")
    threads: lambda w: SUBTYPING_C.get("threads", 8)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        out_dir = lambda w: os.path.join(PY_RES, "subtyping", w.sample),
        metadata = lambda w: py_full_path(SUBTYPING_TASKS[w.sample].get("metadata", "")),
        tree_method = SUBTYPING_C.get("tree_method", "iqtree"),
        script = "scripts/classify_subtypes_placement.py",
        env = py_get_env("subtyping")
    shell:
        """
        set -euo pipefail
        
        # Build optional metadata argument
        META_ARG=""
        if [ -n "{params.metadata}" ] && [ -f "{params.metadata}" ]; then
            META_ARG="--metadata {params.metadata}"
        fi
        
        # Determine conda run parameter flag (name vs path)
        ENV_FLAG="-n"
        if [[ "{params.env}" == /* ]] || [[ "{params.env}" == ./* ]]; then
            ENV_FLAG="-p"
        fi
        
        # Run the Python subtyping script inside the conda env
        conda run $ENV_FLAG "{params.env}" python3 {params.script} \
            --ref-fasta {input.ref_fasta} \
            --query-fasta {input.query_fasta} \
            --out-dir {params.out_dir} \
            --threads {threads} \
            --tree-method {params.tree_method} \
            $META_ARG 2>&1 | tee {log}
        """
