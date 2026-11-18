#!/usr/bin/env python3
"""
Batch Prodigal + HMM Annotation Pipeline for large-scale MAG annotation.

Output structure:
  {output_dir}/{batch_name}/{mag_name}/{mag_name}.faa
  {output_dir}/{batch_name}/{mag_name}/{mag_name}.gff3
  {output_dir}/{batch_name}/{mag_name}/advanced_hmm/unified_all_raw.tbl
  {output_dir}/{batch_name}/{mag_name}/{mag_name}_advanced_hmm_resolved.tsv
  {output_dir}/{batch_name}/{mag_name}/{mag_name}_functional_summary.tsv
"""
import sys
import os
import glob
import re
import math
import shlex
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed


def run_cmd(cmd):
    """运行 shell 命令，不捕获 stdout/stderr，让输出直接流入 Snakemake 的 tee 日志。"""
    try:
        subprocess.run(cmd, shell=True, check=True)
        return True, ""
    except subprocess.CalledProcessError as e:
        return False, f"Command failed with exit code {e.returncode}: {cmd}"


def find_mags(genomes_dir):
    """Finds all DNA genomes in the specified directory/directories using high-performance os.scandir."""
    dirs = [genomes_dir] if isinstance(genomes_dir, str) else genomes_dir
    valid_exts = {".fna", ".fasta", ".fa", ".fna.gz", ".fasta.gz", ".fa.gz"}
    mag_files = {}
    for d in dirs:
        if not os.path.exists(d):
            print(f"WARNING: genomes directory does not exist: {d}", file=sys.stderr, flush=True)
            continue
        try:
            with os.scandir(d) as entries:
                for entry in entries:
                    if entry.is_file():
                        name = entry.name
                        # Skip protein files
                        if name.endswith((".faa", ".faa.gz")):
                            continue
                        # Find matching extension
                        for ext in valid_exts:
                            if name.endswith(ext):
                                # Extract sample name by removing the extension
                                stem = name[:-len(ext)]
                                f_path = entry.path
                                if stem not in mag_files:
                                    mag_files[stem] = f_path
                                break
        except Exception as e:
            print(f"WARNING: Failed to scan directory {d}: {e}", file=sys.stderr, flush=True)
    return mag_files


def run_prodigal_single(mag_name, mag_path, base_out_dir, conda_env, output_ffn=False):
    """Runs prodigal for a single genome.
    
    Args:
        base_out_dir: Base output directory, already includes batch_name prefix,
                      e.g. results/external_annotation/test
    """
    mag_out_dir = os.path.join(base_out_dir, mag_name)
    faa_out = os.path.join(mag_out_dir, f"{mag_name}.faa")
    gff_out = os.path.join(mag_out_dir, f"{mag_name}.gff3")
    ffn_out = os.path.join(mag_out_dir, f"{mag_name}.ffn")

    os.makedirs(mag_out_dir, exist_ok=True)

    # Cache check: if outputs exist and FAA has content, skip
    if os.path.exists(faa_out) and os.path.exists(gff_out) and os.path.getsize(faa_out) > 0:
        if not output_ffn or (os.path.exists(ffn_out) and os.path.getsize(ffn_out) > 0):
            return mag_name, faa_out, gff_out, True

    # Handle gzipped input
    input_fa = mag_path
    temp_fa = None
    if mag_path.endswith(".gz"):
        temp_fa = os.path.join(mag_out_dir, f"temp_{mag_name}.fa")
        run_cmd(f"gunzip -c {shlex.quote(mag_path)} > {shlex.quote(temp_fa)}")
        input_fa = temp_fa

    if output_ffn:
        cmd = f"pyrodigal -i {shlex.quote(input_fa)} -a {shlex.quote(faa_out)} -d {shlex.quote(ffn_out)} -o {shlex.quote(gff_out)} -f gff -p meta"
    else:
        cmd = f"pyrodigal -i {shlex.quote(input_fa)} -a {shlex.quote(faa_out)} -o {shlex.quote(gff_out)} -f gff -p meta"
        
    if conda_env:
        cmd = f"conda run -p {shlex.quote(conda_env)} {cmd}"

    success, err = run_cmd(cmd)

    # Clean up temp
    if temp_fa and os.path.exists(temp_fa):
        os.remove(temp_fa)

    if not success:
        print(f"ERROR: Prodigal failed on {mag_name}: {err}", file=sys.stderr, flush=True)
        return mag_name, None, None, False

    # Ensure output files exist (prodigal may predict 0 genes and not write files)
    check_paths = [faa_out, gff_out]
    if output_ffn:
        check_paths.append(ffn_out)
    for fpath in check_paths:
        if not os.path.exists(fpath):
            with open(fpath, "w"):
                pass

    return mag_name, faa_out, gff_out, False


def combine_faa(batch_genomes, base_out_dir, batch_id, combined_faa_path=None):
    """Combines FAA files of a batch, prepending the genome name to prevent header collisions."""
    if not combined_faa_path:
        combined_faa_path = os.path.join(base_out_dir, f"batch_{batch_id}_combined.faa")
    total_proteins = 0
    with open(combined_faa_path, "w") as out:
        for mag_name, faa_path in batch_genomes:
            if not faa_path or not os.path.exists(faa_path):
                continue
            with open(faa_path) as f:
                for line in f:
                    if line.startswith(">"):
                        orig_header = line[1:].strip()
                        prot_id = orig_header.split()[0]
                        unique_id = f"{mag_name}|{prot_id}"
                        rest = orig_header[len(prot_id):]
                        out.write(f">{unique_id}{rest}\n")
                        total_proteins += 1
                    else:
                        out.write(line)
    print(f"Combined {total_proteins} proteins from {len(batch_genomes)} genomes.", flush=True)
    return combined_faa_path


def split_fasta(input_path, num_chunks, output_dir, batch_id):
    """Splits a FASTA file into num_chunks roughly equal size chunks by writing directly to files."""
    if num_chunks <= 1:
        # Create a copy as chunk 1 to keep cleanup/concat logic uniform
        chunk_path = os.path.join(output_dir, f"batch_{batch_id}_chunk_1.faa")
        import shutil
        shutil.copy2(input_path, chunk_path)
        return [chunk_path]

    chunk_paths = [
        os.path.join(output_dir, f"batch_{batch_id}_chunk_{idx+1}.faa")
        for idx in range(num_chunks)
    ]
    
    # Open all chunk files
    chunk_files = [open(p, "w") for p in chunk_paths]
    
    try:
        chunk_idx = -1
        with open(input_path) as f:
            for line in f:
                if line.startswith(">"):
                    chunk_idx = (chunk_idx + 1) % num_chunks
                if chunk_idx >= 0:
                    chunk_files[chunk_idx].write(line)
    finally:
        for cf in chunk_files:
            cf.close()
            
    # Remove empty chunk files (if any)
    active_paths = []
    for p in chunk_paths:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            active_paths.append(p)
        else:
            if os.path.exists(p):
                os.remove(p)
    return active_paths


def split_hmm_results(combined_tbl, base_out_dir, mag_names, threads=16):
    """Splits the combined HMM table back into individual genome tbl files using a parallel Map-Reduce hash partitioning strategy."""
    import os
    import sys
    import math
    import time
    import shutil
    import threading
    from concurrent.futures import ThreadPoolExecutor, as_completed

    if not os.path.exists(combined_tbl):
        print(f"WARNING: Combined HMM table not found: {combined_tbl}", file=sys.stderr, flush=True)
        return

    file_size = os.path.getsize(combined_tbl)
    if file_size == 0:
        print("Combined HMM table is empty.", flush=True)
        # Create empty header tbl files for all MAGs
        for name in mag_names:
            mag_hmm_dir = os.path.join(base_out_dir, name, "advanced_hmm")
            os.makedirs(mag_hmm_dir, exist_ok=True)
            tbl_path = os.path.join(mag_hmm_dir, "unified_all_raw.tbl")
            with open(tbl_path, "wb") as out:
                out.write(b"# target name        accession   query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu ov env dom rep inc description of target\n")
                out.write(b"#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ ----- ----- --- --- --- --- --- --- --- ---\n")
        return

    print(f"Starting Map-Reduce HMM split for {len(mag_names)} MAGs (File size: {file_size / (1024*1024*1024):.2f} GB, threads={threads})", flush=True)
    start_time = time.time()

    mag_names_set = set(name.encode('utf-8') for name in mag_names)
    num_partitions = 128
    temp_dir = os.path.join(base_out_dir, "_split_temp")
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir, exist_ok=True)

    # Stage 1: Sequential Partitioning (Map) with Block Reads (Bypasses GIL bottlenecks)
    print("[Map Stage] Partitioning combined table into buckets (block-reads)...", flush=True)
    
    # Open all partition files
    part_files = {}
    for p_idx in range(num_partitions):
        p_path = os.path.join(temp_dir, f"part_{p_idx}.tbl")
        part_files[p_idx] = open(p_path, "wb", buffering=10*1024*1024)  # 10MB write buffer

    total_hits = 0
    try:
        import zlib
        with open(combined_tbl, "rb") as f:
            block_size = 16 * 1024 * 1024  # 16MB block size
            block = f.read(block_size)
            leftover = b""
            
            while block:
                if leftover:
                    block = leftover + block
                    leftover = b""
                
                last_newline = block.rfind(b'\n')
                if last_newline == -1:
                    leftover = block
                    block = f.read(block_size)
                    continue
                
                data_to_process = block[:last_newline]
                leftover = block[last_newline+1:]
                
                block = f.read(block_size)
                
                lines = data_to_process.split(b'\n')
                for line in lines:
                    if not line or line.startswith(b'#'):
                        continue
                    pipe_idx = line.find(b'|')
                    if pipe_idx == -1:
                        continue
                    mag_name_bytes = line[:pipe_idx]
                    if mag_name_bytes in mag_names_set:
                        part_idx = zlib.crc32(mag_name_bytes) % num_partitions
                        part_files[part_idx].write(line)
                        part_files[part_idx].write(b'\n')
                        total_hits += 1
                        
            if leftover:
                lines = leftover.split(b'\n')
                for line in lines:
                    if not line or line.startswith(b'#'):
                        continue
                    pipe_idx = line.find(b'|')
                    if pipe_idx == -1:
                        continue
                    mag_name_bytes = line[:pipe_idx]
                    if mag_name_bytes in mag_names_set:
                        part_idx = zlib.crc32(mag_name_bytes) % num_partitions
                        part_files[part_idx].write(line)
                        part_files[part_idx].write(b'\n')
                        total_hits += 1
    finally:
        for pf in part_files.values():
            pf.close()

    map_time = time.time()
    print(f"[Map Stage] Completed. Found {total_hits:,} total hits. Time: {map_time - start_time:.2f} seconds.", flush=True)

    # Stage 2: Parallel Write (Reduce)
    print("[Reduce Stage] Writing hits to individual MAG folders...", flush=True)
    num_reduce_threads = min(threads, 16)  # Cap at 16 to avoid excessive disk write concurrency/thrashing
    
    written_mags = set()
    written_mags_lock = threading.Lock()

    def process_partition(part_idx):
        p_path = os.path.join(temp_dir, f"part_{part_idx}.tbl")
        if not os.path.exists(p_path):
            return set()

        hits = {}
        # Read and group lines in memory
        with open(p_path, "rb") as pf:
            for line in pf:
                pipe_idx = line.find(b'|')
                if pipe_idx == -1:
                    continue
                mag_name_bytes = line[:pipe_idx]
                clean_line = line[pipe_idx+1:]
                hits.setdefault(mag_name_bytes, []).append(clean_line)

        # Write to final files
        local_written = set()
        for mag_name_bytes, lines in hits.items():
            mag_name = mag_name_bytes.decode('utf-8')
            mag_hmm_dir = os.path.join(base_out_dir, mag_name, "advanced_hmm")
            tbl_path = os.path.join(mag_hmm_dir, "unified_all_raw.tbl")

            os.makedirs(mag_hmm_dir, exist_ok=True)
            with open(tbl_path, "wb") as out:
                out.write(b"# target name        accession   query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu ov env dom rep inc description of target\n")
                out.write(b"#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ ----- ----- --- --- --- --- --- --- --- ---\n")
                out.writelines(lines)
            local_written.add(mag_name)

        # Clean up temp file for this partition
        try:
            os.remove(p_path)
        except OSError:
            pass

        return local_written

    with ThreadPoolExecutor(max_workers=num_reduce_threads) as executor:
        futures = [executor.submit(process_partition, p) for p in range(num_partitions)]
        for fut in as_completed(futures):
            try:
                local_written = fut.result()
                if local_written:
                    with written_mags_lock:
                        written_mags.update(local_written)
            except Exception as e:
                print(f"ERROR: Partition reduction failed: {e}", file=sys.stderr, flush=True)
                raise

    reduce_time = time.time()
    print(f"[Reduce Stage] Completed. Wrote hits to {len(written_mags)} MAGs. Time: {reduce_time - map_time:.2f} seconds.", flush=True)

    # Stage 3: Write empty headers for MAGs with zero hits
    print("[Final Stage] Creating empty headers for remaining MAGs...", flush=True)
    remaining_mags = [name for name in mag_names if name not in written_mags]
    
    def write_empty_mag(mag_name):
        mag_hmm_dir = os.path.join(base_out_dir, mag_name, "advanced_hmm")
        tbl_path = os.path.join(mag_hmm_dir, "unified_all_raw.tbl")
        os.makedirs(mag_hmm_dir, exist_ok=True)
        with open(tbl_path, "wb") as out:
            out.write(b"# target name        accession   query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu ov env dom rep inc description of target\n")
            out.write(b"#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ ----- ----- --- --- --- --- --- --- --- ---\n")

    # Use parallel threads to write empty files to speed it up
    if remaining_mags:
        with ThreadPoolExecutor(max_workers=num_reduce_threads) as executor:
            executor.map(write_empty_mag, remaining_mags)

    # Remove temporary directory
    try:
        shutil.rmtree(temp_dir)
    except OSError:
        pass

    end_time = time.time()
    print(f"Parallel Map-Reduce HMM split completed in {end_time - start_time:.2f} seconds. Total hits: {total_hits:,}", flush=True)


def resolve_and_aggregate_single(name, base_out_dir, ko_list_path, run_aggregation, py_exe, biopython_env_path, use_conda, in_mem_context=None, target_profiles="", target_set=None, force_overwrite=False, negative_map="", competition_fraction=0.10):
    """Resolves and aggregates annotations for a single MAG, with support for in-memory and subprocess execution."""
    import os
    import shlex
    import sys
    
    mag_dir = os.path.join(base_out_dir, name)
    unified_tbl = os.path.join(mag_dir, "advanced_hmm/unified_all_raw.tbl")
    faa_file = os.path.join(mag_dir, f"{name}.faa")
    out_resolved = os.path.join(mag_dir, f"{name}_advanced_hmm_resolved.tsv")
    gff_file = os.path.join(mag_dir, f"{name}.gff3")
    out_summary = os.path.join(mag_dir, f"{name}_functional_summary.tsv")

    # Cache check for resolution
    resolve_exists = os.path.exists(out_resolved) and os.path.getsize(out_resolved) > 0
    if force_overwrite:
        resolve_exists = False
    
    # Cache check for aggregation (if run_aggregation is enabled)
    agg_exists = True
    if run_aggregation:
        agg_exists = os.path.exists(out_summary) and os.path.getsize(out_summary) > 0
    if force_overwrite:
        agg_exists = False

    # If both resolution and aggregation results are already completed, skip instantly!
    if resolve_exists and agg_exists:
        return True

    success = True
    
    # 1. HMM Resolution
    if not resolve_exists:
        if in_mem_context and in_mem_context.get("has_pandas"):
            try:
                # In-memory execution: import resolve_single_mag
                from scripts.resolve_annotations import resolve_single_mag
                resolve_single_mag(
                    unified_path=unified_tbl,
                    pfam_path="/dev/null",
                    ko_thresholds=in_mem_context["ko_thresholds"],
                    faa_path=faa_file,
                    output_path=out_resolved,
                    evalue_cutoff=1e-5,
                    log_factor=0.2,
                    target_set=target_set if target_set is not None else set(),
                    target_negatives=(in_mem_context or {}).get("target_negatives", {}),
                    competition_fraction=(in_mem_context or {}).get("competition_fraction", 0.10)
                )
                success = True
            except Exception as e:
                print(f"ERROR: In-memory resolve failed for {name}: {e}", file=sys.stderr, flush=True)
                success = False
        else:
            # Subprocess execution fallback
            if use_conda:
                resolve_cmd = (
                    f"conda run -p {shlex.quote(biopython_env_path)} python3 scripts/resolve_annotations.py "
                    f"--unified {shlex.quote(unified_tbl)} "
                    f"--pfam /dev/null "
                    f"--ko_list {shlex.quote(ko_list_path)} "
                    f"--faa {shlex.quote(faa_file)} "
                    f"--output {shlex.quote(out_resolved)} "
                    f"--evalue_cutoff 1e-5 "
                    f"--negative_map {shlex.quote((in_mem_context or {}).get('negative_map', negative_map))} "
                    f"--competition_fraction {shlex.quote(str((in_mem_context or {}).get('competition_fraction', competition_fraction)))} "
                    f"--target_profiles {shlex.quote(target_profiles)}"
                )
            else:
                resolve_cmd = (
                    f"{shlex.quote(py_exe)} scripts/resolve_annotations.py "
                    f"--unified {shlex.quote(unified_tbl)} "
                    f"--pfam /dev/null "
                    f"--ko_list {shlex.quote(ko_list_path)} "
                    f"--faa {shlex.quote(faa_file)} "
                    f"--output {shlex.quote(out_resolved)} "
                    f"--evalue_cutoff 1e-5 "
                    f"--negative_map {shlex.quote((in_mem_context or {}).get('negative_map', negative_map))} "
                    f"--competition_fraction {shlex.quote(str((in_mem_context or {}).get('competition_fraction', competition_fraction)))} "
                    f"--target_profiles {shlex.quote(target_profiles)}"
                )
            success, err = run_cmd(resolve_cmd)
            if not success:
                print(f"ERROR: subprocess resolve failed for {name}: {err}", file=sys.stderr, flush=True)

    # 2. Aggregation
    if success and run_aggregation and not agg_exists:
        if in_mem_context and in_mem_context.get("has_pandas"):
            try:
                from scripts.aggregate_functional import aggregate_functional
                aggregate_functional(
                    faa_path=faa_file,
                    egg_path="none",
                    kofam_path="none",
                    advanced_hmm_path=out_resolved,
                    gff_path=gff_file,
                    sample_name=name,
                    alias=name,
                    output_path=out_summary
                )
                success = True
            except Exception as e:
                print(f"ERROR: In-memory aggregation failed for {name}: {e}", file=sys.stderr, flush=True)
                success = False
        else:
            # Subprocess execution fallback
            if use_conda:
                agg_cmd = (
                    f"conda run -p {shlex.quote(biopython_env_path)} python3 scripts/aggregate_functional.py "
                    f"{shlex.quote(faa_file)} none none {shlex.quote(out_resolved)} "
                    f"{shlex.quote(gff_file)} {shlex.quote(name)} {shlex.quote(name)} "
                    f"{shlex.quote(out_summary)}"
                )
            else:
                agg_cmd = (
                    f"{shlex.quote(py_exe)} scripts/aggregate_functional.py "
                    f"{shlex.quote(faa_file)} none none {shlex.quote(out_resolved)} "
                    f"{shlex.quote(gff_file)} {shlex.quote(name)} {shlex.quote(name)} "
                    f"{shlex.quote(out_summary)}"
                )
            success_agg, err_agg = run_cmd(agg_cmd)
            if not success_agg:
                print(f"ERROR: subprocess aggregation failed for {name}: {err_agg}", file=sys.stderr, flush=True)
                success = False

    return success


def main():
    p = argparse.ArgumentParser(description="Batch Prodigal & HMM Annotation for MAGs")
    p.add_argument("--genomes-dir", required=True, action="append",
                    help="Directory containing MAG fasta files (can specify multiple)")
    p.add_argument("--output-dir", required=True,
                    help="Base output directory (e.g. results/external_annotation)")
    p.add_argument("--batch-name", required=True,
                    help="Batch/sample name. Outputs go to {output-dir}/{batch-name}/{mag}/")
    p.add_argument("--hmm-db", required=True, help="HMM database path")
    p.add_argument("--ko-list", required=True, help="KEGG ko_list threshold file")
    p.add_argument("--negative-map", default="", help="Target-specific target_to_negative.tsv")
    p.add_argument("--competition-fraction", type=float, default=0.10,
                   help="Near-tie fraction; calibrate against validation data")
    p.add_argument("--conda-env-dir", help="Conda environment directory")
    p.add_argument("--bakta-env", default="bakta_1.12.0")
    p.add_argument("--biopython-env", default="biopython_1.86")
    p.add_argument("--prodigal-env", help="Conda environment for Prodigal (defaults to bakta-env if not specified)")
    p.add_argument("--hmm-split-chunks", type=int, default=20,
                   help="Number of chunks to split the combined proteins for parallel HMM search")
    p.add_argument("--only-predict-combine", action="store_true",
                   help="Only run Prodigal and combine FAAs, then exit")
    p.add_argument("--only-resolve-aggregate", action="store_true",
                   help="Only split HMM results, resolve and aggregate, then exit")
    p.add_argument("--combined-tbl", help="Combined HMM table file (required for only-resolve-aggregate)")
    p.add_argument("--protein-to-mag", help="JSON file mapping protein ID to MAG name")
    p.add_argument("--combined-faa", help="Path to write the combined FAA file")
    p.add_argument("--batch-size", type=int, default=500)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--run-aggregation", action="store_true",
                    help="Run functional annotation aggregation after HMM resolution")
    p.add_argument("--target-profiles", default="",
                    help="Comma-separated list of target profiles to find the minimum E-value for")
    p.add_argument("--force-overwrite", action="store_true",
                    help="Force overwrite existing HMM resolution and aggregation results")
    p.add_argument("--output-ffn", action="store_true",
                    help="Also output nucleotide gene sequences (.ffn) during Pyrodigal prediction")
    args = p.parse_args()

    # Parse target profiles set
    target_set = set()
    if args.target_profiles:
        items = re.split(r'[\s,]+', args.target_profiles)
        target_set = {x.strip() for x in items if x.strip()}

    # Construct the actual output base: {output_dir}/{batch_name}
    base_out_dir = os.path.join(args.output_dir, args.batch_name)
    os.makedirs(base_out_dir, exist_ok=True)

    # 1. Scan MAGs
    if not args.only_resolve_aggregate:
        print("=" * 60, flush=True)
        print(f"Batch: {args.batch_name}", flush=True)
        print("Scanning MAGs...", flush=True)
        mag_files = find_mags(args.genomes_dir)
        print(f"Found {len(mag_files)} MAG files to process.", flush=True)
        if not mag_files:
            print("No MAG files found. Exiting.", flush=True)
            sys.exit(0)

    # Resolve conda env
    bakta_env_path = args.bakta_env
    biopython_env_path = args.biopython_env
    prodigal_env_path = args.prodigal_env if args.prodigal_env else args.bakta_env
    if args.conda_env_dir:
        bakta_env_path = os.path.join(args.conda_env_dir, args.bakta_env)
        biopython_env_path = os.path.join(args.conda_env_dir, args.biopython_env)
        if args.prodigal_env:
            prodigal_env_path = os.path.join(args.conda_env_dir, args.prodigal_env)
        else:
            prodigal_env_path = bakta_env_path

    if args.only_resolve_aggregate:
        if not args.combined_tbl:
            print("ERROR: --combined-tbl is required in --only-resolve-aggregate mode.", file=sys.stderr)
            sys.exit(1)
        print("Scanning genomes directory (this might take a few seconds)...", flush=True)
        # Scan genomes-dir to get ALL MAG names, supporting genomes with 0 predicted genes
        mag_files = find_mags(args.genomes_dir)
        mag_names = sorted(mag_files.keys())
        
        print(f"\n{'=' * 60}", flush=True)
        print(f"--- Resolving and Aggregating HMM results for {len(mag_names)} MAGs ---", flush=True)
        
        # 1. Split HMM results
        print("Splitting HMM results back to individual MAGs...", flush=True)
        split_hmm_results(args.combined_tbl, base_out_dir, mag_names, threads=args.threads)
        
        # 2. Run HMM Resolution (and optional aggregation)
        print(f"Resolving HMM annotations...", flush=True)
        res_count = 0
        resolve_failed = []
        resolve_workers = max(1, args.threads)
        
        # Determine python path inside biopython_env to bypass conda run startup overhead
        py_exe = "python3"
        use_conda = False
        if biopython_env_path:
            candidate = os.path.join(biopython_env_path, "bin/python3")
            if not os.path.exists(candidate):
                candidate = os.path.join(biopython_env_path, "bin/python")
            if os.path.exists(candidate):
                py_exe = candidate
            else:
                use_conda = True

        # Check if pandas is available in the current environment to enable in-memory execution
        has_pandas = False
        in_mem_context = None
        try:
            import pandas as pd
            from scripts.resolve_annotations import load_ko_thresholds
            ko_thresholds = load_ko_thresholds(args.ko_list)
            from scripts.resolve_annotations import load_target_negatives
            in_mem_context = {
                "has_pandas": True,
                "ko_thresholds": ko_thresholds,
                "target_negatives": load_target_negatives(args.negative_map),
                "negative_map": args.negative_map,
                "competition_fraction": args.competition_fraction,
            }
            print("Pandas is importable in parent environment. Enabling high-speed In-Memory HMM resolution.", flush=True)
        except ImportError:
            print("Pandas not found in parent environment. Falling back to high-concurrency subprocesses.", flush=True)

        with ThreadPoolExecutor(max_workers=resolve_workers) as executor:
            resolve_futures = {
                executor.submit(
                    resolve_and_aggregate_single,
                    name=name,
                    base_out_dir=base_out_dir,
                    ko_list_path=args.ko_list,
                    run_aggregation=args.run_aggregation,
                    py_exe=py_exe,
                    biopython_env_path=biopython_env_path,
                    use_conda=use_conda,
                    in_mem_context=in_mem_context,
                    target_profiles=args.target_profiles,
                    target_set=target_set,
                    force_overwrite=args.force_overwrite,
                    negative_map=args.negative_map,
                    competition_fraction=args.competition_fraction
                ): name
                for name in mag_names
            }
            for fut in as_completed(resolve_futures):
                name = resolve_futures[fut]
                try:
                    if fut.result():
                        res_count += 1
                    else:
                        resolve_failed.append(name)
                except Exception as e:
                    print(f"EXCEPTION resolving {name}: {e}", file=sys.stderr, flush=True)
                    resolve_failed.append(name)

        print(f"Done: {res_count} resolved, {len(resolve_failed)} failed.", flush=True)
        if resolve_failed:
            sys.exit(1)
        sys.exit(0)

    if not args.only_resolve_aggregate:
        mag_items = list(mag_files.items())
        num_batches = math.ceil(len(mag_items) / args.batch_size)
    else:
        mag_items = []
        num_batches = 0

    failed_genomes = []
    total_processed = 0

    for b in range(num_batches):
        batch_start = b * args.batch_size
        batch_end = min((b + 1) * args.batch_size, len(mag_items))
        batch_data = mag_items[batch_start:batch_end]
        batch_id = b + 1
        batch_names = [name for name, _ in batch_data]
        print(f"\n{'=' * 60}", flush=True)
        print(f"--- Processing Batch {batch_id}/{num_batches} ({len(batch_data)} MAGs) ---", flush=True)

        # 2. Parallel Prodigal
        print(f"[Step 1/4] Running Prodigal gene prediction (threads={args.threads})...", flush=True)
        batch_genomes = []
        skipped_count = 0
        prodigal_failed = []
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = {
                executor.submit(run_prodigal_single, name, path, base_out_dir, prodigal_env_path, args.output_ffn): name
                for name, path in batch_data
            }
            for fut in as_completed(futures):
                name, faa_out, gff_out, skipped = fut.result()
                if faa_out:
                    batch_genomes.append((name, faa_out))
                    if skipped:
                        skipped_count += 1
                else:
                    prodigal_failed.append(name)

        if prodigal_failed:
            print(f"WARNING: Prodigal failed for {len(prodigal_failed)} genomes: "
                  f"{', '.join(prodigal_failed[:10])}{'...' if len(prodigal_failed) > 10 else ''}",
                  file=sys.stderr, flush=True)
            failed_genomes.extend(prodigal_failed)

        # Count non-empty FAA files
        non_empty_count = sum(1 for _, faa in batch_genomes
                             if faa and os.path.exists(faa) and os.path.getsize(faa) > 0)
        print(f"Prodigal done: {len(batch_genomes)} succeeded ({skipped_count} cached), "
              f"{non_empty_count} with predicted genes, {len(prodigal_failed)} failed.", flush=True)

        if not batch_genomes:
            print(f"No genomes with successful Prodigal in batch {batch_id}. Skipping.", flush=True)
            continue

        # 3. Combine FAA
        print("[Step 2/4] Combining predicted protein sequences...", flush=True)
        combined_faa = combine_faa(batch_genomes, base_out_dir, batch_id, args.combined_faa)

        if args.only_predict_combine:
            print(f"Batch {batch_id} predict and combine done. Skipping HMM search and resolve.", flush=True)
            # Clean up combined FAA is skipped, because we want it to stay
            continue

        # Check if combined FAA is empty
        if os.path.getsize(combined_faa) == 0:
            print(f"WARNING: Combined FAA is empty (0 predicted proteins across all genomes). "
                  f"Skipping HMM for batch {batch_id}.", flush=True)
            # Still create empty tbl files so resolve can run
            split_hmm_results("/dev/null", base_out_dir, batch_names, threads=1)
        else:
            # 4. Run HMMsearch in parallel via splitting combined FAA
            combined_tbl = os.path.join(base_out_dir, f"batch_{batch_id}_combined_hmm.tbl")
            print(f"[Step 3/4] Splitting and running hmmsearch in parallel (threads={args.threads}, chunks={args.hmm_split_chunks})...", flush=True)
            
            # Split combined_faa into chunks
            chunk_faa_paths = split_fasta(combined_faa, args.hmm_split_chunks, base_out_dir, batch_id)
            
            # Run hmmsearch on each chunk in parallel
            num_actual_chunks = len(chunk_faa_paths)
            chunk_tbl_paths = []
            
            # Calculate workers and threads per job
            hmm_workers = min(args.threads, num_actual_chunks)
            cpu_per_job = max(1, args.threads // hmm_workers)
            
            print(f"Running {num_actual_chunks} HMM search jobs (concurrency={hmm_workers}, threads per job={cpu_per_job})...", flush=True)
            
            def run_hmm_chunk(idx, chunk_faa):
                chunk_tbl = os.path.join(base_out_dir, f"batch_{batch_id}_combined_hmm_chunk_{idx+1}.tbl")
                hmm_cmd = (f"hmmsearch --cpu {cpu_per_job} "
                           f"--tblout {shlex.quote(chunk_tbl)} -E 1e-5 "
                           f"{shlex.quote(args.hmm_db)} {shlex.quote(chunk_faa)} > /dev/null")
                if bakta_env_path:
                    hmm_cmd = f"conda run -p {shlex.quote(bakta_env_path)} {hmm_cmd}"
                success, err = run_cmd(hmm_cmd)
                return chunk_tbl, success, err
            
            hmm_failed = False
            failed_err_msg = ""
            with ThreadPoolExecutor(max_workers=hmm_workers) as executor:
                futures = {
                    executor.submit(run_hmm_chunk, idx, path): idx
                    for idx, path in enumerate(chunk_faa_paths)
                }
                for fut in as_completed(futures):
                    idx = futures[fut]
                    chunk_tbl, success, err = fut.result()
                    if success:
                        chunk_tbl_paths.append(chunk_tbl)
                    else:
                        hmm_failed = True
                        failed_err_msg = err
            
            # Clean up chunk FAA files
            for p in chunk_faa_paths:
                if os.path.exists(p):
                    os.remove(p)
                    
            if hmm_failed:
                print(f"FATAL: hmmsearch failed for batch {batch_id}: {failed_err_msg}",
                      file=sys.stderr, flush=True)
                print(f"Marking all {len(batch_names)} genomes in batch as FAILED.",
                      file=sys.stderr, flush=True)
                failed_genomes.extend(batch_names)
                for p in chunk_tbl_paths:
                    if os.path.exists(p):
                        os.remove(p)
                if os.path.exists(combined_faa):
                    os.remove(combined_faa)
                continue
                
            # Concatenate chunk tbl files into combined_tbl
            chunk_tbl_paths.sort()
            with open(combined_tbl, "w") as out:
                for chunk_tbl in chunk_tbl_paths:
                    if os.path.exists(chunk_tbl):
                        with open(chunk_tbl) as f:
                            for line in f:
                                out.write(line)
                                
            # Clean up chunk tbl files
            for p in chunk_tbl_paths:
                if os.path.exists(p):
                    os.remove(p)
                    
            # 5. Split HMM results
            print("Splitting HMM results back to individual MAGs...", flush=True)
            split_hmm_results(combined_tbl, base_out_dir, batch_names, threads=args.threads)
            
            if os.path.exists(combined_tbl):
                os.remove(combined_tbl)

        # Clean up combined FAA
        if not args.only_predict_combine and os.path.exists(combined_faa):
            os.remove(combined_faa)

        # 6. Run HMM Resolution (and optional aggregation)
        print(f"[Step 4/4] Resolving HMM annotations...", flush=True)
        res_count = 0
        resolve_failed = []

        # 降低 resolve 并发度，避免 conda 锁竞争和内存暴涨
        resolve_workers = min(args.threads, 8)

        # Determine python path inside biopython_env to bypass conda run startup overhead
        py_exe = "python3"
        use_conda = False
        if biopython_env_path:
            candidate = os.path.join(biopython_env_path, "bin/python3")
            if not os.path.exists(candidate):
                candidate = os.path.join(biopython_env_path, "bin/python")
            if os.path.exists(candidate):
                py_exe = candidate
            else:
                use_conda = True

        # Check if pandas is available in the current environment to enable in-memory execution
        has_pandas = False
        in_mem_context = None
        try:
            import pandas as pd
            from scripts.resolve_annotations import load_ko_thresholds
            from scripts.resolve_annotations import load_target_negatives
            ko_thresholds = load_ko_thresholds(args.ko_list)
            in_mem_context = {
                "has_pandas": True,
                "ko_thresholds": ko_thresholds,
                "target_negatives": load_target_negatives(args.negative_map),
                "negative_map": args.negative_map,
                "competition_fraction": args.competition_fraction,
            }
            print("Pandas is importable in parent environment. Enabling high-speed In-Memory HMM resolution.", flush=True)
        except ImportError:
            print("Pandas not found in parent environment. Falling back to high-concurrency subprocesses.", flush=True)

        with ThreadPoolExecutor(max_workers=resolve_workers) as executor:
            resolve_futures = {
                executor.submit(
                    resolve_and_aggregate_single,
                    name=name,
                    base_out_dir=base_out_dir,
                    ko_list_path=args.ko_list,
                    run_aggregation=args.run_aggregation,
                    py_exe=py_exe,
                    biopython_env_path=biopython_env_path,
                    use_conda=use_conda,
                    in_mem_context=in_mem_context,
                    target_profiles=args.target_profiles,
                    target_set=target_set,
                    force_overwrite=args.force_overwrite,
                    negative_map=args.negative_map,
                    competition_fraction=args.competition_fraction
                ): name
                for name in batch_names if name not in prodigal_failed
            }
            for fut in as_completed(resolve_futures):
                name = resolve_futures[fut]
                try:
                    if fut.result():
                        res_count += 1
                    else:
                        resolve_failed.append(name)
                except Exception as e:
                    print(f"EXCEPTION resolving {name}: {e}", file=sys.stderr, flush=True)
                    resolve_failed.append(name)

        if resolve_failed:
            failed_genomes.extend(resolve_failed)

        total_processed += res_count
        print(f"Batch {batch_id} done: {res_count} resolved, {len(resolve_failed)} failed.", flush=True)

    # Final summary
    print(f"\n{'=' * 60}", flush=True)
    print(f"FINAL SUMMARY", flush=True)
    print(f"{'=' * 60}", flush=True)
    print(f"Batch name:            {args.batch_name}", flush=True)
    print(f"Total MAGs scanned:    {len(mag_items)}", flush=True)
    print(f"Successfully resolved: {total_processed}", flush=True)
    print(f"Failed:                {len(failed_genomes)}", flush=True)

    if failed_genomes:
        unique_failed = sorted(set(failed_genomes))
        print(f"\nFAILED genomes ({len(unique_failed)}):", file=sys.stderr, flush=True)
        for name in unique_failed:
            print(f"  - {name}", file=sys.stderr, flush=True)
        sys.exit(1)
    else:
        print("\nBatch Annotation completed successfully for all genomes!", flush=True)


if __name__ == "__main__":
    main()
