# ==============================================================================
# Protein Phylogeny Module: MAFFT and IQ-TREE Workflow (Multi-task Edition)
# ==============================================================================
import os
import glob
import re

# 1. Global Path Logic (Namespaced)
# ------------------------------------------------------------------------------
PROT_ROOT = config.get("project_root", "")
PROT_ENV_DIR = config.get("conda_env_dir", "")
PROT_C = config.get("protein_phylogeny", {})

# Load the dictionary of protein phylogeny tasks from samples.yaml
PROT_TASKS = config.get("protein_phylogeny_tasks", {})

def prot_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    # Prefix relative paths with the global project_root
    return os.path.normpath(os.path.join(PROT_ROOT, rel_path))

def prot_get_env(name):
    env_name = config.get("environments", {}).get(name) or PROT_C.get("conda_envs", {}).get(name, name)
    if not PROT_ENV_DIR: return env_name
    return os.path.join(PROT_ENV_DIR, env_name)

# Results directories
PROT_RES = prot_full_path(PROT_C.get("results_dir", "results/protein_phylogeny"))
PROT_LOGS = prot_full_path(PROT_C.get("logs_dir", os.path.join("logs", "protein_phylogeny")))

# Helper function to discover input protein sequences for a specific task
def find_proteins_for_task(task_info):
    fasta_paths = []
    
    # 1. Individual files or glob patterns
    files = task_info.get("fasta_files", [])
    if files:
        if isinstance(files, str):
            files = [files]
        for f in files:
            p = prot_full_path(f)
            if os.path.isfile(p):
                fasta_paths.append(p)
            else:
                for gp in glob.glob(p):
                    if os.path.isfile(gp):
                        fasta_paths.append(gp)
                        
    # 2. Directories to scan for fasta files
    dirs = task_info.get("fasta_dirs", [])
    if dirs:
        if isinstance(dirs, str):
            dirs = [dirs]
        for d in dirs:
            p = prot_full_path(d)
            if os.path.isdir(p):
                # Scan common fasta extensions (and gzipped variants)
                for ext in ["*.faa", "*.fasta", "*.fa", "*.faa.gz", "*.fasta.gz", "*.fa.gz"]:
                    for f_path in glob.glob(os.path.join(p, ext)):
                        if os.path.isfile(f_path):
                            fasta_paths.append(f_path)
                            
    # Deduplicate while preserving order
    seen = set()
    dedup_paths = []
    for gp in fasta_paths:
        abs_p = os.path.abspath(gp)
        if abs_p not in seen:
            seen.add(abs_p)
            dedup_paths.append(abs_p)
            
    return dedup_paths


# 2. Target Registration
# ------------------------------------------------------------------------------
PROT_TARGETS = []

if PROT_C.get("enable", False) and PROT_TASKS:
    # Iterate through configured tasks and register final tree files
    # Output structure: results/protein_phylogeny/iqtree/{sample}/{sample}.iqtree.treefile
    for task_name in PROT_TASKS:
        PROT_TARGETS.append(os.path.join(PROT_RES, f"iqtree/{task_name}/{task_name}.iqtree.treefile"))

TARGETS.extend(PROT_TARGETS)


# 3. Rules
# ------------------------------------------------------------------------------

rule protein_prepare_sequences:
    output:
        raw_faa = os.path.join(PROT_RES, "raw/{sample}_raw.faa")
    log:
        os.path.join(PROT_LOGS, "prepare_{sample}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    run:
        task_name = wildcards.sample
        task_info = PROT_TASKS.get(task_name)
        if not task_info:
            raise ValueError(f"ERROR: Protein phylogeny task '{task_name}' is not defined in samples.yaml")
        
        with open(log[0], "w") as log_f:
            log_f.write(f"=== PROTEIN SEQUENCE INPUT SCAN FOR TASK: {task_name} ===\n")
            log_f.write(f"Global Project Root: {PROT_ROOT}\n")
            log_f.write("Scanning input files and directories...\n")
            
            paths = find_proteins_for_task(task_info)
            if not paths:
                msg = (
                    f"ERROR: No valid FASTA files found for protein phylogeny task '{task_name}'! "
                    "Please configure one of the following under 'protein_phylogeny_tasks' in samples.yaml:\n"
                    "  - 'fasta_files' with direct paths/globs\n"
                    "  - 'fasta_dirs' with directories containing .fa/.fasta/.faa files\n"
                )
                log_f.write(msg)
                raise ValueError(msg)
            
            log_f.write(f"Successfully discovered {len(paths)} FASTA files:\n")
            for p in paths:
                log_f.write(f"  - {p}\n")
            
            # Defensive FASTA reading and header cleaning
            # Sanitizes header fields: removes white space, replaces invalid symbols for trees/MAFFT with underscores,
            # and appends numeric suffixes to resolve duplicates to ensure total name uniqueness.
            log_f.write("\nProcessing and sanitizing sequences...\n")
            seen_headers = set()
            seq_counter = 0
            
            # Ensure the output directory exists
            os.makedirs(os.path.dirname(output.raw_faa), exist_ok=True)
            
            import gzip
            
            with open(output.raw_faa, "w") as out_f:
                for path in paths:
                    # Open with gzip if compressed
                    is_gz = path.endswith(".gz")
                    open_func = gzip.open if is_gz else open
                    mode = "rt" if is_gz else "r"
                    
                    with open_func(path, mode) as in_f:
                        curr_header = ""
                        curr_seq = []
                        for line in in_f:
                            line = line.strip()
                            if not line:
                                continue
                            if line.startswith(">"):
                                if curr_header:
                                    # Process and write previous sequence
                                    seq_counter += 1
                                    write_sanitized_seq(curr_header, curr_seq, out_f, seen_headers, log_f)
                                curr_header = line[1:]
                                curr_seq = []
                            else:
                                curr_seq.append(line)
                        # Don't forget the last sequence
                        if curr_header:
                            seq_counter += 1
                            write_sanitized_seq(curr_header, curr_seq, out_f, seen_headers, log_f)
                            
            log_f.write(f"\nSequence preparation finished. Processed {seq_counter} sequences.\n")
            log_f.write(f"Clean FASTA written to: {output.raw_faa}\n")

# Python helper to sanitize headers and write sequence in FASTA format
def write_sanitized_seq(header, seq_lines, out_f, seen_headers, log_f):
    # Extract only the sequence ID (everything before first whitespace)
    parts = header.split(None, 1)
    seq_id = parts[0]
    
    # 1. Strip and clean invalid characters for Newick tree structures and IQ-TREE parsing:
    # IQ-TREE and Newick format reserve the following characters: ':', ',', ';', '(', ')', '[', ']'
    # We replace them, along with whitespaces and quotes, with underscores.
    clean_id = re.sub(r'[^a-zA-Z0-9_\-\.]', '_', seq_id)
    
    # 2. Prevent duplicate names by appending a suffix if needed
    final_id = clean_id
    counter = 1
    while final_id in seen_headers:
        final_id = f"{clean_id}_{counter}"
        counter += 1
    
    seen_headers.add(final_id)
    
    # Log any naming alterations
    if final_id != seq_id:
        log_f.write(f"  [Sanitize] Renamed '{seq_id}' -> '{final_id}'\n")
        
    out_f.write(f">{final_id}\n")
    out_f.write("".join(seq_lines) + "\n")


rule run_mafft:
    input:
        raw_faa = os.path.join(PROT_RES, "raw/{sample}_raw.faa")
    output:
        msa = os.path.join(PROT_RES, "mafft/{sample}/{sample}.aligned.faa")
    log:
        os.path.join(PROT_LOGS, "mafft/run_{sample}.log")
    threads: get_threads(64)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = prot_get_env("mafft"),
        executable = PROT_C.get("mafft_executable", "mafft"),
        mode = lambda w: PROT_TASKS.get(w.sample, {}).get("mafft_mode") or PROT_C.get("mafft_mode", "--auto"),
        extra = lambda w: PROT_TASKS.get(w.sample, {}).get("mafft_extra") or PROT_C.get("mafft_extra", "")
    shell:
        """
        set -euo pipefail
        
        mkdir -p "$(dirname "{output.msa}")"
        
        echo "=== MAFFT ALIGNMENT START ===" > "{log}"
        echo "Conda Env: {params.env}" >> "{log}"
        echo "Threads: {threads}" >> "{log}"
        echo "Mode: {params.mode}" >> "{log}"
        
        # Execute MAFFT within designated Conda environment.
        # Captures stderr in log and redirects stdout to the aligned output FASTA file.
        conda run -p {params.env} {params.executable} \
            --thread {threads} \
            {params.mode} \
            {params.extra} \
            "{input.raw_faa}" > "{output.msa}" 2>> "{log}"
            
        echo "=== MAFFT ALIGNMENT COMPLETE ===" >> "{log}"
        """


rule run_protein_iqtree:
    input:
        msa = os.path.join(PROT_RES, "mafft/{sample}/{sample}.aligned.faa")
    output:
        treefile = os.path.join(PROT_RES, "iqtree/{sample}/{sample}.iqtree.treefile")
    log:
        os.path.join(PROT_LOGS, "iqtree/{sample}/{sample}.log")
    threads: get_threads(64)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = prot_get_env("iqtree"),
        out_prefix = lambda w: os.path.join(PROT_RES, "iqtree", w.sample, f"{w.sample}.iqtree"),
        mem_g = lambda wildcards, threads: max(1, int(threads * 8)),
        executable = PROT_C.get("iqtree_executable", "iqtree3"),
        bootstrap = lambda w: PROT_TASKS.get(w.sample, {}).get("iqtree_bootstrap") or PROT_C.get("iqtree_bootstrap", "-B 1000"),
        alrt = lambda w: PROT_TASKS.get(w.sample, {}).get("iqtree_alrt") or PROT_C.get("iqtree_alrt", "-alrt 1000"),
        model = lambda w: PROT_TASKS.get(w.sample, {}).get("iqtree_model") or PROT_C.get("iqtree_model", "MFP"),
        extra = lambda w: PROT_TASKS.get(w.sample, {}).get("iqtree_extra") or PROT_C.get("iqtree_extra", "")
    shell:
        """
        set -euo pipefail
        
        mkdir -p "$(dirname "{output.treefile}")"
        
        # Execute IQ-TREE using designated Conda environment.
        # Since this is a single protein family ML tree reconstruction (without partition models),
        # the --mem parameter is fully compatible and safe to use.
        conda run -p {params.env} {params.executable} \
            --mem {params.mem_g}G \
            -T {threads} \
            -m {params.model} \
            -s {input.msa} \
            {params.alrt} \
            {params.bootstrap} \
            {params.extra} \
            --prefix {params.out_prefix} 2>&1 | tee {log}
        """
