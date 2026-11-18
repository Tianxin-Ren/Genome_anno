#!/usr/bin/env python3
import os
import sys
import argparse
import logging
import shutil
import time

# Configure logging
logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")

def transfer_path(src, dst, mode="symlink", dry_run=False):
    """Safely transfers a file or directory from src to dst based on the mode ('symlink' or 'copy').
    Returns the absolute path of the destination if successful, None otherwise.
    """
    if not os.path.exists(src):
        logging.warning(f"Source path does not exist, skipping: {src}")
        return None
    
    src_abs = os.path.abspath(src)
    dst_abs = os.path.abspath(dst)
    
    # Check if destination exists
    if os.path.lexists(dst_abs):
        if mode == "symlink" and os.path.islink(dst_abs):
            link_target = os.readlink(dst_abs)
            if os.path.abspath(link_target) == src_abs:
                logging.info(f"Link already exists and points to target, skipping touch: {dst}")
                return None
        logging.info(f"Destination already exists, skipping migration & touch: {dst}")
        return None

    if dry_run:
        logging.info(f"[DRY-RUN] Would {mode}: {dst} -> {src}")
        return dst_abs
    
    try:
        os.makedirs(os.path.dirname(dst_abs), exist_ok=True)
        if mode == "symlink":
            os.symlink(src_abs, dst_abs)
            logging.info(f"Linked: {dst} -> {src}")
        elif mode == "copy":
            if os.path.isdir(src_abs):
                shutil.copytree(src_abs, dst_abs)
                logging.info(f"Copied directory: {dst} <- {src}")
            else:
                shutil.copy2(src_abs, dst_abs)
                logging.info(f"Copied file: {dst} <- {src}")
        return dst_abs
    except Exception as e:
        logging.error(f"Failed to {transfer_verb(mode)} {dst} from {src}: {e}")
        return None

def transfer_verb(mode):
    return "link" if mode == "symlink" else "copy"

def touch_recursive(path, mtime):
    """Recursively updates the access and modification times of path and all its children."""
    if not os.path.exists(path):
        return
    try:
        os.utime(path, (mtime, mtime))
        if os.path.isdir(path):
            for root, dirs, files in os.walk(path):
                for d in dirs:
                    try:
                        os.utime(os.path.join(root, d), (mtime, mtime))
                    except Exception:
                        pass
                for f in files:
                    try:
                        os.utime(os.path.join(root, f), (mtime, mtime))
                    except Exception:
                        pass
    except Exception as e:
        logging.warning(f"Could not touch {path}: {e}")

def migrate_sample(sample, old_sample_dir, new_results_dir, mode="symlink", dry_run=False):
    logging.info(f"==================== Migrating Sample: {sample} ({mode} mode) ====================")
    
    # Define new module result roots
    pp_dir = os.path.join(new_results_dir, "preprocess", sample)
    as_dir = os.path.join(new_results_dir, "assembly", sample)
    bn_dir = os.path.join(new_results_dir, "binning", sample)
    ab_dir = os.path.join(new_results_dir, "abundance", sample)
    
    # Track paths by stage for chronological touching
    stages = {
        "1_preprocess": [],
        "2_assembly_raw": [],
        "3_assembly_filtered": [],
        "4_binning": [],
        "5_refinement": [],
        "6_abundance": []
    }
    
    # 1. QC (preprocess & abundance prep)
    old_fastp_dir = os.path.join(old_sample_dir, "01_qc", "fastp")
    # Preprocess module mapping
    stages["1_preprocess"].append(transfer_path(os.path.join(old_fastp_dir, f"{sample}_R1.clean.fastq.gz"), os.path.join(pp_dir, "short_read", f"{sample}_R1.clean.fastq.gz"), mode, dry_run))
    stages["1_preprocess"].append(transfer_path(os.path.join(old_fastp_dir, f"{sample}_R2.clean.fastq.gz"), os.path.join(pp_dir, "short_read", f"{sample}_R2.clean.fastq.gz"), mode, dry_run))
    stages["1_preprocess"].append(transfer_path(os.path.join(old_fastp_dir, f"{sample}_fastp.json"), os.path.join(pp_dir, "short_read", f"{sample}.fastp.json"), mode, dry_run))
    stages["1_preprocess"].append(transfer_path(os.path.join(old_fastp_dir, f"{sample}_fastp.html"), os.path.join(pp_dir, "short_read", f"{sample}.fastp.html"), mode, dry_run))
    
    # Abundance module preprocess mapping (abundance expects genomics reads under its folder)
    stages["1_preprocess"].append(transfer_path(os.path.join(old_fastp_dir, f"{sample}_R1.clean.fastq.gz"), os.path.join(ab_dir, "preprocess", "genomics", f"{sample}_R1.clean.fastq.gz"), mode, dry_run))
    stages["1_preprocess"].append(transfer_path(os.path.join(old_fastp_dir, f"{sample}_R2.clean.fastq.gz"), os.path.join(ab_dir, "preprocess", "genomics", f"{sample}_R2.clean.fastq.gz"), mode, dry_run))

    # 2. Assembly (megahit)
    old_megahit_dir = os.path.join(old_sample_dir, "02_assembly", f"{sample}_megahit")
    stages["2_assembly_raw"].append(transfer_path(os.path.join(old_megahit_dir, f"{sample}.contigs.fa"), os.path.join(as_dir, "megahit", f"{sample}.contigs.fa"), mode, dry_run))
    stages["3_assembly_filtered"].append(transfer_path(os.path.join(old_megahit_dir, f"{sample}.contigs.1000.fa"), os.path.join(as_dir, "megahit", f"{sample}.contigs_filtered_1000.fa"), mode, dry_run))
    stages["2_assembly_raw"].append(transfer_path(os.path.join(old_megahit_dir, f"{sample}.log"), os.path.join(new_results_dir, "assembly", "logs", "megahit", f"{sample}.log"), mode, dry_run))

    # 3. Binning (metawrap input bins)
    old_binning_dir = os.path.join(old_sample_dir, "03_binning", f"{sample}_metawrap_bins")
    for bin_type in ["concoct_bins", "maxbin2_bins", "metabat2_bins", "work_files"]:
        stages["4_binning"].append(transfer_path(os.path.join(old_binning_dir, bin_type), os.path.join(bn_dir, "metawrap", bin_type), mode, dry_run))

    # 4. Refinement (metawrap output bins & stats)
    old_refinement_dir = os.path.join(old_sample_dir, "04_refinement", f"{sample}_refined_bins")
    
    # Map refined final bins (metawrap_50_10_bins) to metawrap_50_10_bins
    stages["5_refinement"].append(transfer_path(os.path.join(old_refinement_dir, "metawrap_50_10_bins"), os.path.join(bn_dir, "refinement", "metawrap_50_10_bins"), mode, dry_run))
    stages["5_refinement"].append(transfer_path(os.path.join(old_refinement_dir, "metawrap_50_10_bins.stats"), os.path.join(bn_dir, "refinement", "metawrap_50_10_bins.stats"), mode, dry_run))
    stages["5_refinement"].append(transfer_path(os.path.join(old_refinement_dir, "metawrap_50_10_bins.contigs"), os.path.join(bn_dir, "refinement", "metawrap_50_10_bins.contigs"), mode, dry_run))
    
    # Map intermediate refinement dirs & files
    for folder in ["concoct_bins", "maxbin2_bins", "metabat2_bins", "figures", "work_files"]:
        stages["5_refinement"].append(transfer_path(os.path.join(old_refinement_dir, folder), os.path.join(bn_dir, "refinement", folder), mode, dry_run))
    for file in ["concoct_bins.stats", "concoct_bins.contigs", "maxbin2_bins.stats", "maxbin2_bins.contigs", "metabat2_bins.stats", "metabat2_bins.contigs"]:
        stages["5_refinement"].append(transfer_path(os.path.join(old_refinement_dir, file), os.path.join(bn_dir, "refinement", file), mode, dry_run))

    # 5. Taxonomy (GTDB-Tk results)
    old_tax_dir = os.path.join(old_sample_dir, "05_taxonomy")
    stages["5_refinement"].append(transfer_path(os.path.join(old_tax_dir, "gtdbtk_out"), os.path.join(bn_dir, "taxonomy", "gtdbtk_out"), mode, dry_run))
    stages["5_refinement"].append(transfer_path(os.path.join(old_tax_dir, "gtdbtk.summary.tsv"), os.path.join(bn_dir, "taxonomy", "gtdbtk.summary.tsv"), mode, dry_run))
    
    # 6. Abundance (coverm results)
    old_abundance_dir = os.path.join(old_sample_dir, "06_abundance")
    stages["6_abundance"].append(transfer_path(os.path.join(old_abundance_dir, f"{sample}.coverm.tsv"), os.path.join(ab_dir, "genome_abundance", f"{sample}.coverm.tsv"), mode, dry_run))

    # Apply staged timestamps to prevent Snakemake from thinking inputs are newer than outputs
    if not dry_run:
        base_time = time.time() - 3600  # Start 1 hour ago to be safe
        stage_order = sorted(stages.keys())
        for idx, stage_name in enumerate(stage_order):
            stage_time = base_time + (idx * 60) # Increment by 60s per stage
            logging.info(f"Touching stage {stage_name} files with timestamp offset +{idx * 60}s...")
            for path in stages[stage_name]:
                if path:
                    touch_recursive(path, stage_time)

def main():
    parser = argparse.ArgumentParser(description="Migrate old Genome_anno sample outputs to the new modular results folder structure.")
    parser.add_argument("--old-dir", required=True, help="Path to the directory containing old sample folders (e.g. results/Mariana_metagenomics/samples)")
    parser.add_argument("--new-dir", required=True, help="Path to the new results directory (e.g. results)")
    parser.add_argument("--samples", default="6442,6447,6457,6459,6461,6494", help="Comma-separated list of sample names to migrate")
    parser.add_argument("--mode", choices=["symlink", "copy"], default="symlink", help="Migration mode: 'symlink' (symbolic links, recommended) or 'copy' (physical files copy)")
    parser.add_argument("--dry-run", action="store_true", help="Print the changes that would be made without executing")
    
    args = parser.parse_args()
    
    sample_list = [s.strip() for s in args.samples.split(",") if s.strip()]
    
    logging.info(f"Starting migration for samples: {sample_list}")
    logging.info(f"Old samples directory: {args.old_dir}")
    logging.info(f"New results directory: {args.new_dir}")
    logging.info(f"Transfer mode: {args.mode}")
    
    if args.dry_run:
        logging.info("--- DRY-RUN MODE ENABLED (No changes will be written) ---")
        
    for sample in sample_list:
        old_sample_dir = os.path.join(args.old_dir, sample)
        if not os.path.exists(old_sample_dir):
            logging.error(f"Old sample directory does not exist: {old_sample_dir}. Skipping sample.")
            continue
        
        migrate_sample(sample, old_sample_dir, args.new_dir, args.mode, args.dry_run)
        
    logging.info("Migration check completed.")

if __name__ == "__main__":
    main()
