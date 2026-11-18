#!/usr/bin/env python3
"""
combine_mag_reference.py

Constructs a combined multi-MAG reference index and annotation asset for
genomic/transcriptomic read quantification (Bowtie2/HISAT2 + featureCounts).

Integrates directly with `mag_annotation`:
1. Concatenates reference MAG genomes into a unified multi-FASTA reference (--out-fa).
   Ensures contig headers are unique across all genomes (e.g., prefixing with `{genome}___`
   if needed, or preserving existing canonical `{genome}|...` / `{genome}___...` headers).
2. Generates a synchronized combined GFF3 annotation (--out-gff) matching the exact
   contig headers in the FASTA and canonical gene IDs ({genome}___{contig}___gene00001)
   from `mag_annotation`.
3. Populates GFF attributes with `ID`, `gene_id`, `locus_tag`, and `protein_id`
   for seamless compatibility with featureCounts.
"""

import argparse
import glob
import gzip
import os
import re
import sys


def open_file(path, mode="rt"):
    """Transparently opens plain text or gzip files."""
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None, errors="replace" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None, errors="replace" if "t" in mode else None)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build combined multi-MAG reference FASTA and GFF3 synchronized with mag_annotation."
    )
    parser.add_argument(
        "--ref-dir",
        required=True,
        help="Path to directory containing MAG fasta files or path to a single fasta file.",
    )
    parser.add_argument(
        "--out-fa",
        required=True,
        help="Path to destination combined reference FASTA file.",
    )
    parser.add_argument(
        "--out-gff",
        required=True,
        help="Path to destination combined GFF3 annotation file.",
    )
    parser.add_argument(
        "--mag-ext",
        default="fa",
        help="MAG file extension (e.g. fa, fna, fasta, default: fa).",
    )
    parser.add_argument(
        "--genes-dir",
        default="",
        help="Optional path to directory containing per-MAG {genome}.gff or {genome}.gff.gz from mag_annotation.",
    )
    parser.add_argument(
        "--annotation-tsv",
        default="",
        help="Optional path to mag_annotation output (e.g. protein_KO_best.tsv or protein_KO_all_strict.tsv).",
    )
    return parser.parse_args()


def discover_genome_fastas(ref_dir, mag_ext="fa"):
    """Finds all genome fasta files in ref_dir or returns single file."""
    if os.path.isfile(ref_dir):
        base = os.path.basename(ref_dir)
        for ext in [".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"]:
            if base.endswith(ext):
                base = base[:-len(ext)]
                break
        return [(base, ref_dir)]

    if not os.path.isdir(ref_dir):
        raise FileNotFoundError(f"Reference directory or file not found: {ref_dir}")

    valid_exts = [
        f".{mag_ext}", f".{mag_ext}.gz",
        ".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz"
    ]
    seen_stems = {}
    genomes = []

    for entry in sorted(os.listdir(ref_dir)):
        matched_ext = None
        for ext in valid_exts:
            if entry.endswith(ext) and not entry.endswith((".faa", ".faa.gz", ".ffn", ".ffn.gz", ".gff", ".gff.gz")):
                matched_ext = ext
                break
        if matched_ext:
            stem = entry[:-len(matched_ext)]
            if stem not in seen_stems:
                full_path = os.path.join(ref_dir, entry)
                seen_stems[stem] = full_path
                genomes.append((stem, full_path))

    return genomes


def read_fasta_headers(fasta_path):
    """Extracts all seq IDs from a FASTA file without loading all sequences."""
    seq_ids = []
    with open_file(fasta_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                seq_ids.append(seq_id)
    return seq_ids


def find_genome_gff(genome_id, genes_dir=""):
    """
    Searches for pre-computed {genome_id}.gff[.gz] from mag_annotation.
    Checks:
    1. Direct file matches in genes_dir ({genome_id}.gff, {genome_id}.gff.gz, etc.)
    2. Subdirectories in genes_dir (e.g. {genes_dir}/*/{genome_id}.gff*, {genes_dir}/*/genes/{genome_id}.gff*)
    3. Common repository output locations: results/mag_annotation/*/genes/{genome_id}.gff*
    """
    candidates = []

    if genes_dir:
        # Direct paths
        for ext in [".gff", ".gff.gz", ".genes.gff", ".genes.gff.gz"]:
            candidates.append(os.path.join(genes_dir, f"{genome_id}{ext}"))
            candidates.append(os.path.join(genes_dir, "genes", f"{genome_id}{ext}"))

        # Subdirectories in genes_dir (e.g. if genes_dir points to results/mag_annotation or dataset folder)
        if os.path.isdir(genes_dir):
            for sub in glob.glob(os.path.join(genes_dir, "*", "genes", f"{genome_id}.gff*")):
                candidates.append(sub)
            for sub in glob.glob(os.path.join(genes_dir, "*", f"{genome_id}.gff*")):
                candidates.append(sub)

    # Standard fallback paths within the project or cluster environment
    for pattern in [
        os.path.join("results", "mag_annotation", "*", "genes", f"{genome_id}.gff*"),
        os.path.join("results", "mag_annotation", "genes", f"{genome_id}.gff*"),
        f"/dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin/results/mag_annotation/*/genes/{genome_id}.gff*",
    ]:
        for match in glob.glob(pattern):
            candidates.append(match)

    for c in candidates:
        if os.path.isfile(c):
            return c

    return None


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out_fa)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_gff)), exist_ok=True)

    print(f"[combine_mag_reference] Scanning references in: {args.ref_dir}")
    genomes = discover_genome_fastas(args.ref_dir, args.mag_ext)
    if not genomes:
        raise ValueError(f"No MAG FASTA files discovered in: {args.ref_dir}")
    print(f"[combine_mag_reference] Discovered {len(genomes)} genome(s).")

    # Step 1: Analyze contig namespaces across all genomes
    contig_counts = {}
    for genome_id, fasta_path in genomes:
        for c_id in read_fasta_headers(fasta_path):
            contig_counts[c_id] = contig_counts.get(c_id, 0) + 1

    has_duplicates = any(count > 1 for count in contig_counts.values())
    if has_duplicates:
        print("[combine_mag_reference] Duplicate contig names detected across MAGs. Contigs will be prefixed with '{genome}___'.")
    else:
        print("[combine_mag_reference] All contig headers are globally unique across MAGs.")

    # Step 2: Combine FASTAs and GFFs in lockstep
    total_contigs = 0
    total_genes = 0

    with open(args.out_fa, "w", encoding="utf-8") as out_fa_f, \
         open(args.out_gff, "w", encoding="utf-8") as out_gff_f:

        out_gff_f.write("##gff-version 3\n")

        for genome_id, fasta_path in genomes:
            contig_rename_map = {}

            # Read FASTA records
            with open_file(fasta_path, "rt") as in_fa:
                for line in in_fa:
                    if line.startswith(">"):
                        orig_id = line[1:].strip().split()[0]
                        # Decide target contig ID in combined reference
                        if has_duplicates and not (orig_id.startswith(f"{genome_id}___") or orig_id.startswith(f"{genome_id}|")):
                            target_id = f"{genome_id}___{orig_id}"
                        else:
                            target_id = orig_id
                        contig_rename_map[orig_id] = target_id
                        out_fa_f.write(f">{target_id}\n")
                        total_contigs += 1
                    else:
                        out_fa_f.write(line)

            # Resolve GFF for this genome (strictly from mag_annotation, never re-predict)
            gff_path = find_genome_gff(genome_id, args.genes_dir)
            if not gff_path:
                raise FileNotFoundError(
                    f"Pre-computed gene annotation GFF for genome '{genome_id}' was not found. "
                    f"Per strict requirement, re-predicting genes is forbidden (mag_annotation results must be used). "
                    f"Searched genes_dir='{args.genes_dir}' and standard results/mag_annotation locations. "
                    f"Please ensure mag_annotation has completed and provide its genes directory via --genes-dir."
                )

            with open_file(gff_path, "rt") as in_gff:
                for line in in_gff:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\r\n").split("\t")
                    if len(parts) >= 9:
                        orig_c_id = parts[0]
                        target_c_id = contig_rename_map.get(orig_c_id, orig_c_id)
                        parts[0] = target_c_id

                        # Parse existing attributes
                        attr_str = parts[8]
                        attrs = {}
                        for item in attr_str.split(";"):
                            item = item.strip()
                            if "=" in item:
                                k, v = item.split("=", 1)
                                attrs[k.strip()] = v.strip()
                            elif " " in item:
                                k, v = item.split(" ", 1)
                                attrs[k.strip()] = v.strip().strip('"')

                        gene_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("locus_tag") or f"{genome_id}___{target_c_id}___cds"
                        protein_id = attrs.get("protein_id") or gene_id

                        # Ensure standard comprehensive attributes for featureCounts compatibility
                        parts[8] = f"ID={gene_id};gene_id={gene_id};locus_tag={gene_id};protein_id={protein_id}"
                        out_gff_f.write("\t".join(parts) + "\n")
                        total_genes += 1

    print(f"[combine_mag_reference] Successfully generated combined assets:")
    print(f"  - FASTA: {args.out_fa} ({total_contigs} contigs)")
    print(f"  - GFF3:  {args.out_gff} ({total_genes} CDS features)")


if __name__ == "__main__":
    main()
