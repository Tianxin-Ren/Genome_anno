#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import re
import math
import shutil
import gzip
import json
import glob
import shlex
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Mobile element filtering blacklist
MOBILE_PATTERN = re.compile(
    r"(transposase|integrase|insertion\s*sequence|IS[0-9]+|mobile|phage|recombinase|repeat|transposon|resolvase|invertase)",
    re.IGNORECASE
)

def parse_args():
    parser = argparse.ArgumentParser(description="Stage 0: Three-phase reference selection, consensus anchor discovery, and negative ref validation.")
    parser.add_argument("--target_name", required=True, help="Target genome sample ID")
    parser.add_argument("--target_fa", required=True, help="Target genome assembly FASTA")
    parser.add_argument("--target_faa", default="", help="Target genome predicted proteome FAA (optional)")
    parser.add_argument("--target_taxonomy", required=True, help="Target GTDB taxonomy string or GTDB-Tk summary file")
    parser.add_argument("--target_completeness", type=float, default=100.0, help="CheckM completeness")
    parser.add_argument("--target_contamination", type=float, default=0.0, help="CheckM contamination")
    parser.add_argument("--gtdb_metadata", required=True, help="Path to GTDB metadata TSV (e.g. ar53_metadata_r232.tsv)")
    parser.add_argument("--gtdb_genomes_dir", required=True, help="Directory containing GTDB archaeal genomes (*.fna / *.fna.gz)")
    parser.add_argument("--gtdb_proteins_dir", default="", help="Directory containing pre-annotated GTDB archaeal proteomes (*.faa / *.faa.gz)")
    parser.add_argument("--only_representative", action="store_true", default=False, help="Filter GTDB metadata to representative genomes only")
    parser.add_argument("--query_operon_faa", default="", help="FASTA of query gene cluster / operon proteins (optional if target_kos provided)")
    parser.add_argument("--target_kos", default="", help="Comma-separated list of target KO IDs (e.g. K02117,K02118,...)")
    parser.add_argument("--gtdb_ko_table", default="", help="Path to GTDB archaea KO annotations TSV")
    parser.add_argument("--ko_index_cache", default="", help="Path to pre-filtered target KO index cache TSV")
    parser.add_argument("--file_index_cache", default="", help="Path to pre-built GTDB genome and protein file index JSON cache")
    parser.add_argument("--threads", type=int, default=4, help="CPU threads for parallel tools (Skani)")
    parser.add_argument("--min_af", type=float, default=0.50, help="Minimum aligned fraction threshold (AF_min)")
    parser.add_argument("--adaptive_af", action="store_true", default=True, help="Adaptively scale AF_min by completeness for MAGs")
    parser.add_argument("--min_ani", type=float, default=80.0, help="Minimum ANI for Skani (fallback to taxonomy if lower)")
    parser.add_argument("--max_positive_refs", type=int, default=10, help="Maximum confirmed positive references to select")
    parser.add_argument("--max_negative_refs", type=int, default=2, help="Maximum confirmed negative references to select")
    parser.add_argument("--skani_env", default="", help="Conda environment for skani")
    parser.add_argument("--diamond_env", default="", help="Conda environment for diamond")
    parser.add_argument("--module_id", default="M00159", help="ATP synthase module ID (e.g. M00159, M00157, M00158, M00160)")
    parser.add_argument("--module_name", default="", help="ATP synthase module human-readable name")
    parser.add_argument("--domain", default="", choices=["", "archaea", "bacteria"], help="Organism domain (archaea or bacteria; auto-detected from target_taxonomy if empty)")
    parser.add_argument("--auto_fetch", dest="auto_fetch", action="store_true", default=None, help="Enable on-demand fetching of bacterial references if none found on disk")
    parser.add_argument("--no_auto_fetch", dest="auto_fetch", action="store_false", help="Disable on-demand fetching")
    parser.add_argument("--out_refs_tsv", required=True, help="Output TSV of selected & validated references")
    parser.add_argument("--out_anchors_tsv", required=True, help="Output TSV of consensus orthologous anchors")
    parser.add_argument("--out_seed_operon_faa", default="", help="Output copy/extracted seed query operon FAA")
    parser.add_argument("--out_dir", required=True, help="Output directory to stage reference data")
    return parser.parse_args()

def load_fasta_records(path):
    """Safely loads FASTA records whether uncompressed or gzipped."""
    if not path or not os.path.exists(path):
        return []
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as h:
            return list(SeqIO.parse(h, "fasta"))
    with open(path, "r", errors="replace") as handle:
        return list(SeqIO.parse(handle, "fasta"))

def parse_taxonomy_string(tax_input):
    """Extracts taxonomy dictionary from raw string or GTDB-Tk summary file."""
    tax_str = ""
    if os.path.exists(tax_input):
        try:
            df = pd.read_csv(tax_input, sep="\t")
            if "classification" in df.columns:
                tax_str = str(df["classification"].iloc[0])
            elif "GTDB_taxonomy" in df.columns:
                tax_str = str(df["GTDB_taxonomy"].iloc[0])
        except Exception:
            tax_str = ""
    if not tax_str:
        tax_str = tax_input.strip()
    
    ranks = {"d": "", "p": "", "c": "", "o": "", "f": "", "g": "", "s": ""}
    parts = [p.strip() for p in tax_str.split(";")]
    for p in parts:
        if "__" in p:
            r, name = p.split("__", 1)
            if r in ranks:
                ranks[r] = name
    return ranks, tax_str

def extract_accession_keys(identifier):
    """Generates all possible accession keys and variants for robust dictionary lookup."""
    keys = set()
    s = str(identifier).strip()
    if not s:
        return keys
    keys.add(s)
    clean = s.replace("GB_", "").replace("RS_", "")
    keys.add(clean)
    
    # Check for GCA_ / GCF_ pattern
    m = re.search(r"(GC[AF])_(\d+\.\d+)", s)
    if m:
        prefix, num = m.group(1), m.group(2)
        keys.add(f"{prefix}_{num}")
        keys.add(f"GCA_{num}")
        keys.add(f"GCF_{num}")
        keys.add(f"GB_GCA_{num}")
        keys.add(f"RS_GCF_{num}")
        keys.add(num)
        
        # Also without version number
        num_no_ver = num.split(".")[0]
        keys.add(f"{prefix}_{num_no_ver}")
        keys.add(f"GCA_{num_no_ver}")
        keys.add(f"GCF_{num_no_ver}")
        keys.add(num_no_ver)
    return keys

def index_directory(directory, valid_extensions, cache_file="", cached_all=None):
    """
    Recursively walks a directory (following symlinks) and maps all possible accession keys to the absolute file path.
    If cache_file or cached_all is provided, loads/returns the cached index in O(1) time without filesystem scanning.
    """
    if not directory:
        print(f"Warning: Indexing skipped because directory path is empty.", file=sys.stderr)
        return {}

    # Check in-memory cached_all or on-disk cache_file first
    if cached_all is not None and directory in cached_all and isinstance(cached_all[directory], dict) and cached_all[directory]:
        return cached_all[directory]

    if cache_file and os.path.exists(cache_file) and os.path.getsize(cache_file) > 10:
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                loaded = json.load(f)
            if directory in loaded and isinstance(loaded[directory], dict) and loaded[directory]:
                if cached_all is not None:
                    cached_all[directory] = loaded[directory]
                return loaded[directory]
        except Exception:
            pass

    index = {}
    if not os.path.exists(directory):
        print(f"Warning: Directory does NOT exist on disk: {directory}", file=sys.stderr)
        return index
    if not os.path.isdir(directory):
        print(f"Warning: Path is not a directory: {directory}", file=sys.stderr)
        return index

    total_scanned = 0
    sample_files = []
    for root, _, files in os.walk(directory, followlinks=True):
        for fname in files:
            total_scanned += 1
            if len(sample_files) < 5:
                sample_files.append(fname)
            if any(fname.endswith(ext) for ext in valid_extensions):
                full_p = os.path.join(root, fname)
                keys = extract_accession_keys(fname)
                for k in keys:
                    if k not in index:
                        index[k] = full_p
                # Also index the stem without valid extension
                stem = fname
                for ext in valid_extensions:
                    if stem.endswith(ext):
                        stem = stem[:-len(ext)]
                        break
                index[stem] = full_p
                index[f"GB_{stem}"] = full_p
                index[f"RS_{stem}"] = full_p

    if total_scanned == 0:
        print(f"Warning: Scanned {directory} (followlinks=True), but found 0 files!", file=sys.stderr)
    elif len(index) == 0:
        print(f"Warning: Found {total_scanned} files in {directory}, but NONE matched extensions {valid_extensions}! First few files: {sample_files}", file=sys.stderr)

    # Persist newly built index to cache_file
    if cache_file and index:
        try:
            os.makedirs(os.path.dirname(os.path.abspath(cache_file)), exist_ok=True)
            existing = {}
            if os.path.exists(cache_file) and os.path.getsize(cache_file) > 10:
                try:
                    with open(cache_file, "r", encoding="utf-8") as f:
                        existing = json.load(f)
                except Exception:
                    existing = {}
            existing[directory] = index
            if cached_all is not None:
                cached_all[directory] = index
            tmp_f = cache_file + f".tmp.{os.getpid()}"
            with open(tmp_f, "w", encoding="utf-8") as f:
                json.dump(existing, f)
            os.replace(tmp_f, cache_file)
        except Exception as e:
            print(f"Notice: Could not write file index cache {cache_file}: {e}", file=sys.stderr)

    return index


def find_indexed_file(acc, index_dict, fallback_dir="", valid_extensions=()):
    """Looks up file path from in-memory index or falls back to direct path checks."""
    keys = extract_accession_keys(acc)
    for k in keys:
        if k in index_dict and os.path.exists(index_dict[k]):
            return index_dict[k]
            
    if fallback_dir and os.path.isdir(fallback_dir):
        for k in keys:
            for ext in valid_extensions:
                cand = os.path.join(fallback_dir, f"{k}{ext}")
                if os.path.exists(cand):
                    return cand
    return ""

def index_protein_dir(proteins_dir):
    """Compatibility alias for indexing proteome directory."""
    return index_directory(proteins_dir, (".faa", ".faa.gz"))

def find_protein_file(acc, gtdb_proteins_dir, protein_index=None):
    """Compatibility alias for looking up proteome file."""
    idx = protein_index if protein_index is not None else index_protein_dir(gtdb_proteins_dir)
    return find_indexed_file(acc, idx, gtdb_proteins_dir, (".faa", ".faa.gz"))

def find_genome_file(acc, gtdb_genomes_dir, genome_index=None):
    """Compatibility alias for looking up genome file."""
    idx = genome_index if genome_index is not None else index_directory(gtdb_genomes_dir, (".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz"))
    return find_indexed_file(acc, idx, gtdb_genomes_dir, (".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz"))


def ensure_uncompressed_faa(src_faa_path, out_dir, acc):
    """Ensures a decompressed FAA is available in out_dir if source is gzipped."""
    if not src_faa_path or not os.path.exists(src_faa_path):
        return ""
    if not src_faa_path.endswith(".gz"):
        return src_faa_path

    clean_acc = acc.replace("GB_", "").replace("RS_", "")
    target_faa = os.path.join(out_dir, f"{clean_acc}.faa")
    if not os.path.exists(target_faa) or os.path.getsize(target_faa) == 0:
        with gzip.open(src_faa_path, "rb") as f_in:
            with open(target_faa, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
    return target_faa

def run_diamond_blastp(query_faa, db_faa, out_tsv, diamond_env, evalue=1e-5, identity=30.0, coverage=40.0):
    """Runs DIAMOND blastp between query and db."""
    cmd_prefix = f"conda run -p {diamond_env} " if diamond_env and os.path.isabs(diamond_env) else (f"conda run -n {diamond_env} " if diamond_env else "")
    db_prefix = out_tsv + ".dmnd"
    subprocess.run(
        f"{cmd_prefix}diamond makedb --in {db_faa} -d {db_prefix} --quiet",
        shell=True, check=True
    )
    subprocess.run(
        f"{cmd_prefix}diamond blastp -q {query_faa} -d {db_prefix} -o {out_tsv} -e {evalue} --id {identity} --query-cover {coverage} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --quiet",
        shell=True, check=True
    )
    if os.path.exists(db_prefix + ".dmnd"):
        os.remove(db_prefix + ".dmnd")

def run_skani_search(target_fa, candidate_fastas, out_tsv, skani_env, threads=4):
    """Executes Skani dist to compute ANI and Aligned Fractions with strict fail-fast error checking."""
    if not os.path.exists(target_fa):
        raise FileNotFoundError(f"[select_and_validate_refs] Target FASTA does not exist: {target_fa}")
    valid_candidates = list(dict.fromkeys(p for p in candidate_fastas if p and os.path.exists(p)))
    if not valid_candidates:
        print(f"[{os.path.basename(target_fa)}] Notice: No valid candidate FASTA files found on disk for Skani search.", file=sys.stderr)
        return False
    cmd_prefix = f"conda run -p {skani_env} " if skani_env and os.path.isabs(skani_env) else (f"conda run -n {skani_env} " if skani_env else "")
    t_arg = f"-t {threads}" if threads and threads > 1 else ""
    cand_args = " ".join(shlex.quote(p) for p in valid_candidates)
    cmd = f"{cmd_prefix}skani dist -q {shlex.quote(target_fa)} -r {cand_args} -o {shlex.quote(out_tsv)} {t_arg}"
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        err_msg = res.stderr.strip() if res.stderr else f"exit code {res.returncode}"
        raise RuntimeError(f"Skani execution failed with exit code {res.returncode}: {err_msg}")
    return True

def extract_flanking_anchors(seed_records, operon_pids, coords=None):
    """Identifies conserved flanking genes L2, L1, R1, R2 around the operon on the same contig."""
    if coords:
        op_records = [record for record in seed_records if record.id in operon_pids and record.id in coords]
        contigs = {coords[record.id][0] for record in op_records}
        if len(contigs) == 1 and op_records:
            contig = next(iter(contigs))
            start = min(coords[record.id][1] for record in op_records)
            end = max(coords[record.id][2] for record in op_records)
            flank_records = [record for record in seed_records
                             if record.id not in operon_pids and record.id in coords
                             and coords[record.id][0] == contig
                             and not MOBILE_PATTERN.search(record.description)
                             and not MOBILE_PATTERN.search(record.id)]
            left = sorted((record for record in flank_records if coords[record.id][2] < start),
                          key=lambda record: coords[record.id][2], reverse=True)
            right = sorted((record for record in flank_records if coords[record.id][1] > end),
                           key=lambda record: coords[record.id][1])
            return (left[1].id if len(left) > 1 else "Anchor_L2",
                    left[0].id if left else "Anchor_L1",
                    right[0].id if right else "Anchor_R1",
                    right[1].id if len(right) > 1 else "Anchor_R2")
    def get_contig_prefix(rec_id):
        if "___" in rec_id:
            return rec_id.rsplit("___", 1)[0]
        return rec_id.rsplit("_", 1)[0]

    ctg_to_indices = {}
    for idx, r in enumerate(seed_records):
        if r.id in operon_pids or any(op in r.id or r.id in op for op in operon_pids):
            ctg = get_contig_prefix(r.id)
            ctg_to_indices.setdefault(ctg, []).append(idx)

    if not ctg_to_indices:
        # Fallback if indices couldn't be resolved directly
        non_mobile = [r for r in seed_records if not MOBILE_PATTERN.search(r.description) and not MOBILE_PATTERN.search(r.id)]
        l2 = non_mobile[0].id if len(non_mobile) > 0 else "Anchor_L2"
        l1 = non_mobile[1].id if len(non_mobile) > 1 else "Anchor_L1"
        r1 = non_mobile[2].id if len(non_mobile) > 2 else "Anchor_R1"
        r2 = non_mobile[3].id if len(non_mobile) > 3 else "Anchor_R2"
        return l2, l1, r1, r2

    # Pick the contig that has the maximum number of operon genes (the primary operon contig)
    target_contig = max(ctg_to_indices.keys(), key=lambda c: len(ctg_to_indices[c]))
    ctg_indices = sorted(ctg_to_indices[target_contig])
    min_idx = ctg_indices[0]
    max_idx = ctg_indices[-1]

    # Search upstream for L1, L2 on target_contig
    l1_rec = None
    l2_rec = None
    for idx in range(min_idx - 1, -1, -1):
        r = seed_records[idx]
        if get_contig_prefix(r.id) == target_contig:
            if not MOBILE_PATTERN.search(r.description) and not MOBILE_PATTERN.search(r.id):
                if l1_rec is None:
                    l1_rec = r
                elif l2_rec is None:
                    l2_rec = r
                    break
        else:
            break

    # Search downstream for R1, R2 on target_contig
    r1_rec = None
    r2_rec = None
    for idx in range(max_idx + 1, len(seed_records)):
        r = seed_records[idx]
        if get_contig_prefix(r.id) == target_contig:
            if not MOBILE_PATTERN.search(r.description) and not MOBILE_PATTERN.search(r.id):
                if r1_rec is None:
                    r1_rec = r
                elif r2_rec is None:
                    r2_rec = r
                    break
        else:
            break

    # Fallbacks if operon is at contig edge
    l1 = l1_rec.id if l1_rec else "Anchor_L1"
    l2 = l2_rec.id if l2_rec else (l1_rec.id if l1_rec else "Anchor_L2")
    r1 = r1_rec.id if r1_rec else "Anchor_R1"
    r2 = r2_rec.id if r2_rec else (r1_rec.id if r1_rec else "Anchor_R2")
    return l2, l1, r1, r2

def reference_gene_coordinates(faa_path, source_path=""):
    """Use observed FAA coordinates only; synthetic protein order is not evidence."""
    coords = {}
    for record in load_fasta_records(faa_path):
        parts = record.description.split("#")
        if len(parts) < 4:
            continue
        try:
            start, end, strand = int(parts[1]), int(parts[2]), int(parts[3])
        except ValueError:
            continue
        contig = record.id.rsplit("___", 1)[0] if "___" in record.id else record.id.rsplit("_", 1)[0]
        coords[record.id] = (contig, start, end, "+" if strand >= 0 else "-")
    clean_name = os.path.basename(faa_path).split(".faa")[0]
    base = os.path.dirname(source_path or faa_path)
    candidates = []
    for directory in (base, os.path.join(os.path.dirname(base), "genes"),
                      os.path.join(os.path.dirname(os.path.dirname(base)), "genes")):
        candidates.extend(glob.glob(os.path.join(directory, f"*{clean_name}*.gff*")))
        if source_path:
            original_stem = os.path.basename(source_path).split(".faa")[0]
            candidates.extend(glob.glob(os.path.join(directory, f"*{original_stem}*.gff*")))
    if candidates:
        path = sorted(candidates)[0]
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", errors="replace") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "CDS":
                    continue
                attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
                entry = (fields[0], int(fields[3]), int(fields[4]), fields[6])
                for key in ("protein_id", "ID", "Name", "locus_tag"):
                    if attrs.get(key):
                        coords[attrs[key]] = entry
    return coords

def reference_hits(candidate, query_faa, query_ids, diamond_env, work_dir):
    """Return one distinct protein hit per query gene with observed coordinates."""
    faa = candidate.get("faa_path", "")
    coords = reference_gene_coordinates(faa, candidate.get("faa", ""))
    hits = {}
    if query_faa and os.path.exists(query_faa) and os.path.getsize(query_faa) and faa and os.path.exists(faa):
        out = os.path.join(work_dir, f"{candidate['accession'].replace('/', '_')}.query_hits.tsv")
        run_diamond_blastp(query_faa, faa, out, diamond_env)
        if os.path.exists(out) and os.path.getsize(out):
            frame = pd.read_csv(out, sep="\t", header=None)
            used = set()
            for _, row in frame.sort_values(11, ascending=False).iterrows():
                qid, pid = str(row[0]), str(row[1])
                if qid in hits or pid in used or qid not in query_ids or pid not in coords:
                    continue
                hits[qid] = pid
                used.add(pid)
    elif candidate.get("ko_proteins"):
        for ko, pid in candidate["ko_proteins"].items():
            if ko in query_ids and pid in coords:
                hits[ko] = pid
    return hits, coords

def select_colocalized_hits(hits, coords, minimum, max_span=50000):
    """Select the largest supported module-gene subset in one observed locus.

    Module membership and operon structure are evaluated separately: a module
    may include subunits encoded outside its main locus. This gate therefore
    requires a majority locus cluster, while Stage 2 tests strand/order/gaps.
    """
    by_contig = {}
    for qid, pid in hits.items():
        if pid in coords:
            contig, start, end, strand = coords[pid]
            by_contig.setdefault(contig, []).append((start, end, qid, pid, strand))
    best = {}
    for entries in by_contig.values():
        entries.sort(key=lambda item: (item[0], item[1]))
        for left in range(len(entries)):
            chosen = {}
            seen_loci = set()
            for right in range(left, len(entries)):
                start, end, qid, pid, strand = entries[right]
                if end - entries[left][0] > max_span:
                    break
                locus = (start, end)
                if locus in seen_loci:
                    continue
                seen_loci.add(locus)
                chosen[qid] = pid
                if len(chosen) > len(best):
                    best = dict(chosen)
    return best if len(best) >= minimum else {}

def clustered_hits(hits, coords, minimum, max_span=50000):
    return len(select_colocalized_hits(hits, coords, minimum, max_span)) >= minimum

def reference_locus_threshold(total_query_genes):
    """Require at least 40% of module genes, with a minimum of three when possible."""
    if total_query_genes <= 0:
        return 0
    return min(total_query_genes, max(3, math.ceil(total_query_genes * 0.4)))

def confirm_negative_reference(candidate, anchor_records, seed_faa, diamond_env, work_dir):
    """Require bilateral anchor closure and no translated DNA homolog anywhere."""
    faa = candidate.get("faa_path", "")
    genome = candidate.get("fasta", "")
    if not faa or not genome or not os.path.isfile(faa) or not os.path.isfile(genome):
        return False, "Reference sequence unavailable"
    wanted = {"L1", "R1"}
    if set(anchor_records) != wanted:
        return False, "Validated positive reference lacks bilateral anchors"
    anchor_faa = os.path.join(work_dir, "negative_check_anchors.faa")
    with open(anchor_faa, "w") as handle:
        for key in ("L1", "R1"):
            SeqIO.write(SeqRecord(anchor_records[key].seq, id=key, description=""), handle, "fasta")
    hits, coords = reference_hits(candidate, anchor_faa, wanted, diamond_env, work_dir)
    if set(hits) != wanted:
        return False, "Both flank anchors not detected"
    left, right = coords[hits["L1"]], coords[hits["R1"]]
    if left[0] != right[0] or abs(right[1] - left[2]) > 100000:
        return False, "Flank anchors do not close one locus"
    prefix = f"conda run -p {shlex.quote(diamond_env)} " if diamond_env and os.path.isabs(diamond_env) else (f"conda run -n {shlex.quote(diamond_env)} " if diamond_env else "")
    db = os.path.join(work_dir, "negative_check_seed")
    out = os.path.join(work_dir, f"{candidate['accession'].replace('/', '_')}.dna_hits.tsv")
    for command in (f"{prefix}diamond makedb --in {shlex.quote(seed_faa)} -d {shlex.quote(db)} --quiet",
                    f"{prefix}diamond blastx -q {shlex.quote(genome)} -d {shlex.quote(db)} -o {shlex.quote(out)} -e 1e-4 --subject-cover 30 --quiet"):
        run = subprocess.run(command, shell=True, capture_output=True, text=True)
        if run.returncode:
            return False, f"DNA residue search failed: {run.stderr[:200]}"
    if os.path.getsize(out):
        return False, "Translated genome contains target homologous residue"
    return True, "Bilateral anchor closure with no protein or translated DNA target hit"

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    
    # 1. Parse query operon / target KOs
    target_ko_list = [k.strip() for k in args.target_kos.split(",") if k.strip()] if args.target_kos else []
    query_records = []
    user_provided_query_faa = bool(args.query_operon_faa and os.path.exists(args.query_operon_faa) and os.path.getsize(args.query_operon_faa) > 0)
    if user_provided_query_faa:
        query_records = list(SeqIO.parse(args.query_operon_faa, "fasta"))

    query_gene_ids = [r.id for r in query_records] if user_provided_query_faa else target_ko_list
    total_query_genes = len(query_gene_ids)
    if total_query_genes == 0:
        raise ValueError("Neither valid query_operon_faa nor target_kos was provided!")
    if args.out_seed_operon_faa and not user_provided_query_faa and os.path.abspath(args.out_seed_operon_faa) != os.path.abspath(args.query_operon_faa or ""):
        os.makedirs(os.path.dirname(os.path.abspath(args.out_seed_operon_faa)), exist_ok=True)
        open(args.out_seed_operon_faa, "w").close()
    elif args.out_seed_operon_faa and user_provided_query_faa and os.path.abspath(args.out_seed_operon_faa) != os.path.abspath(args.query_operon_faa):
        os.makedirs(os.path.dirname(os.path.abspath(args.out_seed_operon_faa)), exist_ok=True)
        shutil.copyfile(args.query_operon_faa, args.out_seed_operon_faa)
    
    # Resolve file index cache path if not explicitly provided
    if not args.file_index_cache and args.out_refs_tsv:
        default_cache = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(args.out_refs_tsv)), "..", "..", "..", "cache", "gtdb_file_index.json"))
        if os.path.exists(default_cache):
            args.file_index_cache = default_cache

    genome_exts = (".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")
    genome_index = index_directory(args.gtdb_genomes_dir, genome_exts, cache_file=args.file_index_cache)
    unique_genomes = len(set(genome_index.values()))
    print(f"[{args.target_name}] Indexed {len(genome_index)} genome accession keys across {unique_genomes} unique genome files from {args.gtdb_genomes_dir}", flush=True)

    protein_exts = (".faa", ".faa.gz")
    protein_index = index_directory(args.gtdb_proteins_dir, protein_exts, cache_file=args.file_index_cache)
    if protein_index:
        unique_proteins = len(set(protein_index.values()))
        print(f"[{args.target_name}] Indexed {len(protein_index)} proteome accession keys across {unique_proteins} unique proteome files from {args.gtdb_proteins_dir}", flush=True)

    # 3. Parse target taxonomy
    ranks, tax_str = parse_taxonomy_string(args.target_taxonomy)
    print(f"[{args.target_name}] Target taxonomy parsed: {tax_str}", flush=True)
    
    # Dynamic AF threshold for MAGs
    effective_min_af = args.min_af
    if args.adaptive_af and args.target_completeness < 95.0:
        factor = max(0.5, args.target_completeness / 100.0)
        effective_min_af = round(args.min_af * factor, 3)
        print(f"[{args.target_name}] Adaptive AF active for MAG (completeness={args.target_completeness}%): effective AF_min threshold = {effective_min_af}", flush=True)

    # 4. Load GTDB Metadata
    if not os.path.exists(args.gtdb_metadata):
        raise FileNotFoundError(f"GTDB metadata file not found at: {args.gtdb_metadata}")
        
    print(f"[{args.target_name}] Reading GTDB metadata: {args.gtdb_metadata}...", flush=True)
    with open(args.gtdb_metadata, "r", encoding="utf-8", errors="ignore") as _f:
        header_cols = [c.strip() for c in _f.readline().split("\t")]
    candidate_cols = {
        "accession", "gtdb_taxonomy", "taxonomy", 
        "checkm_completeness", "completeness", 
        "checkm_contamination", "contamination",
        "gtdb_representative", "representative"
    }
    use_cols = [c for c in header_cols if c in candidate_cols]
    if not use_cols or len(use_cols) < 2:
        use_cols = None
    meta = pd.read_csv(args.gtdb_metadata, sep="\t", usecols=use_cols, low_memory=False)
    acc_col = "accession" if "accession" in meta.columns else meta.columns[0]
    tax_col = "gtdb_taxonomy" if "gtdb_taxonomy" in meta.columns else ("taxonomy" if "taxonomy" in meta.columns else meta.columns[1])
    comp_col = "checkm_completeness" if "checkm_completeness" in meta.columns else ("completeness" if "completeness" in meta.columns else "")
    contam_col = "checkm_contamination" if "checkm_contamination" in meta.columns else ("contamination" if "contamination" in meta.columns else "")

    if args.only_representative:
        rep_col = None
        for col in ["gtdb_representative", "representative"]:
            if col in meta.columns:
                rep_col = col
                break
        if rep_col:
            before_cnt = len(meta)
            meta = meta[meta[rep_col].astype(str).str.lower().isin(["t", "true", "1"])].copy()
            print(f"[{args.target_name}] Filtered GTDB metadata by {rep_col} == 't': {before_cnt} -> {len(meta)} representative genomes.", flush=True)

    # 5. Parse GTDB KO table (if in KO-driven mode) with smart micro-index caching
    ko_table_lookup = {}  # acc_key -> set(kos)
    ko_protein_map = {}   # acc_key -> {ko: protein_id}
    
    if target_ko_list and (args.ko_index_cache or (args.gtdb_ko_table and os.path.exists(args.gtdb_ko_table))):
        target_ko_set = {k.strip() for k in target_ko_list if k.strip()}
        
        # Determine cache file path
        cache_file = args.ko_index_cache
        if not cache_file and args.out_refs_tsv:
            potential_cache = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(args.out_refs_tsv)), "..", "..", "cache", "target_kos_gtdb_index.tsv"))
            cache_file = potential_cache

        cache_loaded = False
        if cache_file and os.path.exists(cache_file) and os.path.getsize(cache_file) > 0:
            try:
                print(f"[{args.target_name}] Loading pre-filtered target KO index cache: {cache_file} ...", flush=True)
                df_cache = pd.read_csv(cache_file, sep="\t", low_memory=False)
                g_col = "genome_id" if "genome_id" in df_cache.columns else df_cache.columns[2]
                k_col = "ko" if "ko" in df_cache.columns else "KO"
                p_col = "protein_id" if "protein_id" in df_cache.columns else df_cache.columns[0]
                matched_rows_cnt = 0
                for _, row in df_cache.iterrows():
                    ko_val = str(row[k_col]).strip()
                    if ko_val not in target_ko_set:
                        continue
                    raw_gid = str(row[g_col])
                    pid_val = str(row[p_col]).strip()
                    keys = extract_accession_keys(raw_gid)
                    for k in keys:
                        if k not in ko_table_lookup:
                            ko_table_lookup[k] = set()
                            ko_protein_map[k] = {}
                        ko_table_lookup[k].add(ko_val)
                        ko_protein_map[k][ko_val] = pid_val
                    matched_rows_cnt += 1
                print(f"[{args.target_name}] Fast-loaded {matched_rows_cnt} KO annotations across {len(ko_table_lookup)} GTDB genomes from cache (instant).", flush=True)
                cache_loaded = True
            except Exception as e:
                print(f"Notice: reading cache file failed ({e}), falling back to master KO table.", file=sys.stderr)
                cache_loaded = False

        if not cache_loaded and args.gtdb_ko_table and os.path.exists(args.gtdb_ko_table):
            print(f"[{args.target_name}] Fast-streaming GTDB KO summary: {args.gtdb_ko_table} for {len(target_ko_list)} target KOs...", flush=True)
            matched_cache_lines = []
            header_saved = ""
            opener = gzip.open if args.gtdb_ko_table.endswith(".gz") else open
            try:
                with opener(args.gtdb_ko_table, "rt", encoding="utf-8", errors="ignore") as f:
                    header_line = f.readline()
                    header_saved = header_line
                    header = [c.strip().lower() for c in header_line.strip().split("\t")]
                    prot_idx = header.index("protein_id") if "protein_id" in header else 0
                    genome_idx = header.index("genome_id") if "genome_id" in header else (2 if len(header) > 2 else 0)
                    ko_idx = header.index("ko") if "ko" in header else 1
                    
                    target_ko_tuple = tuple(target_ko_set)
                    for line in f:
                        if any(k in line for k in target_ko_tuple):
                            parts = line.strip().split("\t")
                            if len(parts) > max(prot_idx, genome_idx, ko_idx):
                                ko_val = parts[ko_idx].strip()
                                if ko_val in target_ko_set:
                                    matched_cache_lines.append(line)
                                    raw_gid = parts[genome_idx].strip()
                                    pid_val = parts[prot_idx].strip()
                                    keys = extract_accession_keys(raw_gid)
                                    for k in keys:
                                        if k not in ko_table_lookup:
                                            ko_table_lookup[k] = set()
                                            ko_protein_map[k] = {}
                                        ko_table_lookup[k].add(ko_val)
                                        ko_protein_map[k][ko_val] = pid_val

                print(f"[{args.target_name}] Successfully mapped {len(matched_cache_lines)} KO hits across {len(ko_table_lookup)} GTDB genomes.", flush=True)
                
                # Automatically save cache for subsequent samples
                if cache_file and matched_cache_lines:
                    try:
                        os.makedirs(os.path.dirname(os.path.abspath(cache_file)), exist_ok=True)
                        tmp_cache = cache_file + f".tmp.{os.getpid()}"
                        with open(tmp_cache, "w", encoding="utf-8") as out_c:
                            out_c.write(header_saved)
                            for ml in matched_cache_lines:
                                out_c.write(ml if ml.endswith("\n") else ml + "\n")
                        os.replace(tmp_cache, cache_file)
                        print(f"[{args.target_name}] Saved target KO micro-index cache to: {cache_file} ({len(matched_cache_lines)} rows)", flush=True)
                    except Exception as ce:
                        print(f"Notice: could not save cache to {cache_file}: {ce}", file=sys.stderr)
            except Exception as e:
                print(f"Warning parsing gtdb_ko_table: {e}", file=sys.stderr)

    # 6. Evaluate all GTDB genomes that exist on disk (matching either genome or proteome)
    matched_meta_rows = []
    for _, row in meta.iterrows():
        acc = str(row[acc_col])
        f_path = find_indexed_file(acc, genome_index, args.gtdb_genomes_dir, genome_exts)
        p_path = find_indexed_file(acc, protein_index, args.gtdb_proteins_dir, protein_exts)
        if not f_path and not p_path:
            continue

        
        comp = float(row[comp_col]) if comp_col and pd.notnull(row[comp_col]) else 100.0
        contam = float(row[contam_col]) if contam_col and pd.notnull(row[contam_col]) else 0.0
        quality_score = comp - 5.0 * contam

        # Gather KO annotations
        matched_kos = set()
        matched_proteins = {}
        for k in extract_accession_keys(acc):
            if k in ko_table_lookup:
                matched_kos.update(ko_table_lookup[k])
                matched_proteins.update(ko_protein_map[k])

        matched_meta_rows.append({
            "accession": acc,
            "taxonomy": str(row[tax_col]),
            "fasta": f_path,
            "faa": p_path,
            "completeness": comp,
            "contamination": contam,
            "quality_score": quality_score,
            "detected_count": len(matched_kos),
            "matched_kos": list(matched_kos),
            "ko_proteins": matched_proteins
        })

    print(f"[{args.target_name}] Matched {len(matched_meta_rows)} GTDB reference genomes on disk.", flush=True)

    # Determine organism domain
    detected_domain = args.domain.lower() if args.domain else ""
    if not detected_domain:
        if "d__Archaea" in tax_str:
            detected_domain = "archaea"
        else:
            detected_domain = "bacteria"
    
    # Auto-fetch flag resolution:
    # If not explicitly specified on CLI, enable for bacteria, disable for archaea
    should_auto_fetch = args.auto_fetch if args.auto_fetch is not None else (detected_domain == "bacteria")

    # Safety fallback if disk lookup failed completely (Bacteria on-demand mode only)
    if not matched_meta_rows and detected_domain == "bacteria" and should_auto_fetch:
        fetcher_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fetch_bacterial_references.py")
        if os.path.exists(fetcher_script) and (args.query_operon_faa or target_ko_list):
            print(f"[{args.target_name}] [Auto-Fetch] 0 local references matched for bacterial target. Triggering on-demand fetch_bacterial_references.py...", flush=True)
            fetch_cmd = [
                sys.executable,
                fetcher_script,
                "--target_taxonomy", tax_str,
                "--gtdb_metadata", args.gtdb_metadata,
                "--out_genomes_dir", args.gtdb_genomes_dir,
                "--out_proteins_dir", args.gtdb_proteins_dir,
                "--max_positive_refs", str(args.max_positive_refs),
                "--max_negative_refs", str(args.max_negative_refs)
            ]
            if args.query_operon_faa:
                fetch_cmd.extend(["--query_operon_faa", args.query_operon_faa])
            if args.diamond_env:
                fetch_cmd.extend(["--conda_env", args.diamond_env])
            try:
                subprocess.run(fetch_cmd, check=False)
                genome_index = index_directory(args.gtdb_genomes_dir, genome_exts, cache_file=args.file_index_cache)
                protein_index = index_directory(args.gtdb_proteins_dir, protein_exts, cache_file=args.file_index_cache)
                for _, row in meta.iterrows():
                    acc = str(row[acc_col])
                    f_path = find_indexed_file(acc, genome_index, args.gtdb_genomes_dir, genome_exts)
                    p_path = find_indexed_file(acc, protein_index, args.gtdb_proteins_dir, protein_exts)
                    if not f_path and not p_path:
                        continue
                    comp = float(row[comp_col]) if comp_col and pd.notnull(row[comp_col]) else 100.0
                    contam = float(row[contam_col]) if contam_col and pd.notnull(row[contam_col]) else 0.0
                    quality_score = comp - 5.0 * contam
                    matched_kos = set()
                    matched_proteins = {}
                    for k in extract_accession_keys(acc):
                        if k in ko_table_lookup:
                            matched_kos.update(ko_table_lookup[k])
                            matched_proteins.update(ko_protein_map[k])
                    matched_meta_rows.append({
                        "accession": acc,
                        "taxonomy": str(row[tax_col]),
                        "fasta": f_path,
                        "faa": p_path,
                        "completeness": comp,
                        "contamination": contam,
                        "quality_score": quality_score,
                        "detected_count": len(matched_kos),
                        "matched_kos": list(matched_kos),
                        "ko_proteins": matched_proteins
                    })
                print(f"[{args.target_name}] Auto-fetch completed: now matched {len(matched_meta_rows)} GTDB reference genomes on disk.", flush=True)
            except Exception as e:
                print(f"[{args.target_name}] Auto-fetch failed or offline: {e}", file=sys.stderr)

    if not matched_meta_rows:
        sample_meta_accs = meta[acc_col].head(5).tolist() if not meta.empty else []
        print(f"ERROR: 0 GTDB reference genomes matched on disk!", file=sys.stderr)
        print(f"  -> domain: {detected_domain}", file=sys.stderr)
        print(f"  -> gtdb_genomes_dir: {args.gtdb_genomes_dir} (indexed keys: {len(genome_index)})", file=sys.stderr)
        print(f"  -> GTDB metadata rows: {len(meta)}, sample accessions from meta: {sample_meta_accs}", file=sys.stderr)
        if detected_domain == "archaea":
            print(f"  -> Hint (Archaea mode): Please ensure archaea reference genomes (*.fna) exist in {args.gtdb_genomes_dir} or check ar53 metadata accession matching.", file=sys.stderr)
        else:
            print(f"  -> Hint (Bacteria mode): Run 'python3 scripts/fetch_bacterial_references.py' to acquire real reference genomes.", file=sys.stderr)
        print("No reference will be fabricated from the target assembly.", file=sys.stderr)


    # 7. Stratify into Positive Pool and Negative Pool
    # Search broadly enough to retain modules with peripheral subunits; Stage 0
    # later requires this many genes to have real coordinates in one locus.
    positive_threshold = reference_locus_threshold(total_query_genes)
    positive_pool = []
    negative_pool = []

    if ko_table_lookup:
        for r in matched_meta_rows:
            if r["detected_count"] >= positive_threshold:
                positive_pool.append(r)
            elif r["detected_count"] == 0:
                negative_pool.append(r)
    else:
        # If KO table not provided, all matched rows are candidates for Diamond screening
        positive_pool = matched_meta_rows.copy()
        negative_pool = matched_meta_rows.copy()

    # 8. Hierarchical Selection of Confirmed Positives (Phylogenetic Nearest Neighbors with Operon)
    # Tiers: Genus -> Family -> Order -> Class -> Phylum -> Domain
    tax_hierarchy = [
        ("Tier_A_Genus", f"g__{ranks['g']};" if ranks["g"] else ""),
        ("Tier_B_Family", f"f__{ranks['f']};" if ranks["f"] else ""),
        ("Tier_C_Order", f"o__{ranks['o']};" if ranks["o"] else ""),
        ("Tier_D_Class", f"c__{ranks['c']};" if ranks["c"] else ""),
        ("Tier_E_Phylum", f"p__{ranks['p']};" if ranks["p"] else ""),
        ("Tier_F_Domain", f"d__{ranks['d']};" if ranks["d"] else "d__Archaea;"),
    ]

    selected_positives = []
    positive_screen_cap = max(args.max_positive_refs, 5 * args.max_positive_refs)
    found_positive_tier = ""

    for tier_name, tier_pattern in tax_hierarchy:
        if not tier_pattern:
            continue
        tier_candidates = [r for r in positive_pool if tier_pattern in r["taxonomy"]]
        if tier_candidates:
            found_positive_tier = tier_name
            tier_candidates.sort(key=lambda x: (x.get("detected_count", 0), x.get("quality_score", 0)), reverse=True)
            for cand in tier_candidates:
                if cand not in selected_positives:
                    cand["tier"] = tier_name
                    selected_positives.append(cand)
                    if len(selected_positives) >= positive_screen_cap:
                        break
        if len(selected_positives) >= positive_screen_cap:
            break

    # If still fewer than max_positive_refs, gather remaining top positives across domain
    if len(selected_positives) < positive_screen_cap and positive_pool:
        positive_pool.sort(key=lambda x: (x.get("detected_count", 0), x.get("quality_score", 0)), reverse=True)
        for cand in positive_pool:
            if cand not in selected_positives:
                cand["tier"] = cand.get("tier", "Tier_F_Domain")
                selected_positives.append(cand)
                if len(selected_positives) >= positive_screen_cap:
                    break

    print(f"[{args.target_name}] Selected {len(selected_positives)} positive references (closest tier: {found_positive_tier or 'None'}).", flush=True)

    # 9. Hierarchical Selection of Confirmed Negatives (Closest relatives lacking the operon)
    selected_negatives = []
    negative_screen_cap = max(args.max_negative_refs, 5 * args.max_negative_refs)
    for tier_name, tier_pattern in tax_hierarchy[:3]:
        if not tier_pattern:
            continue
        tier_candidates = [r for r in negative_pool if tier_pattern in r["taxonomy"] and r not in selected_positives]
        if tier_candidates:
            tier_candidates.sort(key=lambda x: x.get("quality_score", 0), reverse=True)
            for cand in tier_candidates:
                if cand not in selected_negatives:
                    cand["tier"] = tier_name
                    selected_negatives.append(cand)
                    if len(selected_negatives) >= negative_screen_cap:
                        break
        if len(selected_negatives) >= negative_screen_cap:
            break

    if len(selected_negatives) < negative_screen_cap and negative_pool:
        for cand in negative_pool:
            if cand not in selected_positives and cand not in selected_negatives:
                cand["tier"] = cand.get("tier", "Tier_C_Order")
                selected_negatives.append(cand)
                if len(selected_negatives) >= negative_screen_cap:
                    break

    print(f"[{args.target_name}] Selected {len(selected_negatives)} negative references.", flush=True)

    # 10. Ensure proteome FAA files are uncompressed and ready
    all_selected = selected_positives + selected_negatives
    for c in all_selected:
        acc = c["accession"]
        fa_file = c["fasta"]
        src_faa = c.get("faa", "")
        if not src_faa:
            src_faa = find_indexed_file(acc, protein_index, args.gtdb_proteins_dir, protein_exts)
            
        if src_faa:
            faa_file = ensure_uncompressed_faa(src_faa, args.out_dir, acc)
        else:
            faa_file = fa_file.replace(".fna", ".faa").replace(".fasta", ".faa").replace(".fa", ".faa")
            if not os.path.exists(faa_file) or os.path.getsize(faa_file) == 0:
                faa_file = os.path.join(args.out_dir, f"{acc}.faa")
                if not os.path.exists(faa_file) and os.path.exists(fa_file) and fa_file != args.target_fa:
                    cmd_dmnd = f"conda run -p {args.diamond_env} " if args.diamond_env and os.path.isabs(args.diamond_env) else (f"conda run -n {args.diamond_env} " if args.diamond_env else "")
                    cmd_prod = f"{cmd_dmnd}prodigal -i {fa_file} -a {faa_file} -o /dev/null -p single -q"
                    subprocess.run(cmd_prod, shell=True)
        c["faa_path"] = faa_file

    # 11. Run Skani search on selected candidates
    fastas = [c["fasta"] for c in all_selected if c.get("fasta") and os.path.exists(c["fasta"]) and c["fasta"] != args.target_fa]
    skani_out = os.path.join(args.out_dir, "skani_search_results.tsv")
    skani_map = {}
    if fastas:
        run_skani_search(args.target_fa, fastas, skani_out, args.skani_env, threads=args.threads)
        if os.path.exists(skani_out) and os.path.getsize(skani_out) > 0:
            try:
                sk_df = pd.read_csv(skani_out, sep="\t")
                for _, r in sk_df.iterrows():
                    f_ref = str(r.get("Ref_file", r.iloc[0] if len(r) > 0 else ""))
                    f_q = str(r.get("Query_file", ""))
                    try:
                        ani = float(r.get("ANI", 0.0))
                    except (ValueError, TypeError):
                        ani = 0.0
                    try:
                        raw_af_ref = float(r.get("Align_fraction_ref", 0.0))
                    except (ValueError, TypeError):
                        raw_af_ref = 0.0
                    try:
                        raw_af_q = float(r.get("Align_fraction_query", 0.0))
                    except (ValueError, TypeError):
                        raw_af_q = 0.0

                    af_ref_val = raw_af_ref / 100.0 if raw_af_ref > 1.0 else raw_af_ref
                    af_q_val = raw_af_q / 100.0 if raw_af_q > 1.0 else raw_af_q

                    for c in all_selected:
                        acc = c["accession"]
                        c_fa = c.get("fasta", "")
                        c_base = os.path.basename(c_fa) if c_fa else ""
                        clean_acc = acc.replace("GB_", "").replace("RS_", "")
                        acc_no_ver = clean_acc.split('.')[0] if '.' in clean_acc else clean_acc

                        matched_as_ref = (c_fa and (c_fa in f_ref or f_ref in c_fa)) or (c_base and c_base in f_ref) or (clean_acc and clean_acc in f_ref) or (len(acc_no_ver) > 6 and acc_no_ver in f_ref)
                        matched_as_query = (c_fa and (c_fa in f_q or f_q in c_fa)) or (c_base and c_base in f_q) or (clean_acc and clean_acc in f_q) or (len(acc_no_ver) > 6 and acc_no_ver in f_q)

                        if matched_as_ref:
                            skani_map[acc] = {
                                "ANI": round(ani, 2),
                                "AF_ref": round(af_ref_val, 4),
                                "AF_query": round(af_q_val, 4),
                                "AF_min": round(min(af_ref_val, af_q_val), 4)
                            }
                        elif matched_as_query:
                            skani_map[acc] = {
                                "ANI": round(ani, 2),
                                "AF_ref": round(af_q_val, 4),
                                "AF_query": round(af_ref_val, 4),
                                "AF_min": round(min(af_ref_val, af_q_val), 4)
                            }
            except Exception as e:
                raise RuntimeError(f"Error parsing Skani results {skani_out}: {e}")

    # Assign real Skani metrics or explicit 'NA' (beyond Skani limit / no alignment).
    # NEVER fabricate synthetic ANI/AF estimates or defensive constants!
    for c in all_selected:
        acc = c["accession"]
        if acc in skani_map:
            c.update(skani_map[acc])
        else:
            c.update({
                "ANI": "NA",
                "AF_ref": "NA",
                "AF_query": "NA",
                "AF_min": "NA"
            })

    # Re-sort selected positives by combined genomic distance (Skani ANI * AF_min) and gene completeness.
    # Genomes without Skani hits (ANI < 80% or unaligned) have distance product 0.0 and rank by detected_count.
    def _pos_rank_key(x):
        ani = x.get("ANI")
        af = x.get("AF_min")
        ani_val = float(ani) if (ani is not None and ani != "NA") else 0.0
        af_val = float(af) if (af is not None and af != "NA") else 0.0
        return (
            ani_val * af_val,
            int(x.get("detected_count", 0)),
            float(x.get("quality_score", 0.0))
        )

    selected_positives.sort(key=_pos_rank_key, reverse=True)
    print(f"[{args.target_name}] Re-ranked positive references by genomic proximity (ANI * AF): {[p['accession'] for p in selected_positives]}", flush=True)

    # A KO count alone does not establish a conserved locus. Validate every
    # selected positive on observed genomic coordinates before choosing a seed.
    validated_positives = []
    provisional_positives = []
    # Stage 0 establishes a reference locus, not a complete operon. Some KEGG
    # modules include subunits encoded outside their main locus.
    min_locus_genes = positive_threshold
    for candidate in selected_positives:
        try:
            hits, coords = reference_hits(candidate,
                                           args.query_operon_faa if user_provided_query_faa else "",
                                           set(query_gene_ids), args.diamond_env, args.out_dir)
            cluster_hits = select_colocalized_hits(hits, coords, min_locus_genes)
            candidate["cluster_query_hits"] = cluster_hits
            candidate["cluster_pids"] = set(cluster_hits.values())
            candidate["detected_count"] = len(hits)
            candidate["locus_cluster_gene_count"] = len(cluster_hits)
            if len(cluster_hits) >= min_locus_genes:
                candidate["final_validation_status"] = "CONFIRMED_POSITIVE"
                candidate["validation_note"] = (f"{len(cluster_hits)}/{total_query_genes} module genes co-localized within 50 kb; "
                                                 "operon structure assessed separately in Stage 2")
                validated_positives.append(candidate)
            else:
                candidate["final_validation_status"] = "PROVISIONAL_POSITIVE"
                candidate["validation_note"] = f"Only {len(cluster_hits)}/{total_query_genes} module genes co-localized; need {min_locus_genes}"
                provisional_positives.append(candidate)
        except Exception as exc:
            candidate["final_validation_status"] = "PROVISIONAL_POSITIVE"
            candidate["validation_note"] = f"Reference validation failed: {exc}"
            provisional_positives.append(candidate)
    selected_positives = validated_positives[:args.max_positive_refs]
    print(f"[{args.target_name}] Positive-reference audit: scanned={len(matched_meta_rows)}, "
          f"gene-screen_candidates={len(positive_pool)}, selected_for_coordinate_check={len(validated_positives) + len(provisional_positives)}, "
          f"confirmed_loci={len(selected_positives)}, minimum_locus_genes={min_locus_genes}.", flush=True)

    # 12. Extract Authentic Seed Query Operon FAA from Top Confirmed Positive Reference
    seed_pids = set()
    if selected_positives:
        seed_p = selected_positives[0]
        seed_proteins = seed_p.get("cluster_query_hits", {})
        seed_faa = seed_p.get("faa_path", "")
        extracted_recs = []

        if not user_provided_query_faa and seed_proteins and os.path.exists(seed_faa):
            seed_pids = set(seed_p.get("cluster_pids", set()))
            seed_coords = reference_gene_coordinates(seed_faa, seed_p.get("faa", ""))
            for r in sorted(load_fasta_records(seed_faa),
                            key=lambda rec: (seed_coords.get(rec.id, ("", float("inf")))[0],
                                             seed_coords.get(rec.id, ("", float("inf")))[1])):
                for ko, pid in seed_proteins.items():
                    if pid == r.id or pid in r.id or r.id in pid:
                        r_copy = SeqRecord(r.seq, id=f"{ko}_{r.id}", description=f"KO:{ko} | from {seed_p['accession']}")
                        extracted_recs.append(r_copy)
                        break

        if extracted_recs:
            if args.query_operon_faa:
                os.makedirs(os.path.dirname(os.path.abspath(args.query_operon_faa)), exist_ok=True)
                with open(args.query_operon_faa, "w") as f_qout:
                    SeqIO.write(extracted_recs, f_qout, "fasta")
            if args.out_seed_operon_faa:
                os.makedirs(os.path.dirname(os.path.abspath(args.out_seed_operon_faa)), exist_ok=True)
                with open(args.out_seed_operon_faa, "w") as f_qout:
                    SeqIO.write(extracted_recs, f_qout, "fasta")
            print(f"[{args.target_name}] Extracted {len(extracted_recs)} authentic seed operon proteins from {seed_p['accession']}.", flush=True)

    # A KO identifier is not a protein sequence. Without a validated positive
    # reference, keep the seed empty unless the user supplied real proteins.
    if args.out_seed_operon_faa and not selected_positives and not user_provided_query_faa:
        os.makedirs(os.path.dirname(os.path.abspath(args.out_seed_operon_faa)), exist_ok=True)
        open(args.out_seed_operon_faa, "w").close()
    elif args.out_seed_operon_faa and (not os.path.exists(args.out_seed_operon_faa) or os.path.getsize(args.out_seed_operon_faa) == 0):
        os.makedirs(os.path.dirname(os.path.abspath(args.out_seed_operon_faa)), exist_ok=True)
        if user_provided_query_faa:
            shutil.copyfile(args.query_operon_faa, args.out_seed_operon_faa)
        else:
            open(args.out_seed_operon_faa, "w").close()

    # 13. Extract Consensus Anchors from Top Confirmed Positive Reference
    consensus_anchors = []
    anchor_records = {}
    if selected_positives:
        seed_ref = selected_positives[0]
        seed_faa = seed_ref.get("faa_path", "")
        seed_records = load_fasta_records(seed_faa)
        if not seed_pids:
            seed_pids = set(seed_ref.get("cluster_pids", set()))

        l2, l1, r1, r2 = extract_flanking_anchors(
            seed_records, seed_pids, reference_gene_coordinates(seed_faa, seed_ref.get("faa", "")))
        seed_by_id = {record.id: record for record in seed_records}
        anchor_records = {key: seed_by_id[pid] for key, pid in (("L1", l1), ("R1", r1)) if pid in seed_by_id}
        consensus_anchors = [
            {"anchor_id": "L2", "side": "left", "rank": 2, "gene_symbol": l2, "product": f"conserved_anchor_L2_from_{seed_ref['accession']}"},
            {"anchor_id": "L1", "side": "left", "rank": 1, "gene_symbol": l1, "product": f"conserved_anchor_L1_from_{seed_ref['accession']}"},
            {"anchor_id": "R1", "side": "right", "rank": 1, "gene_symbol": r1, "product": f"conserved_anchor_R1_from_{seed_ref['accession']}"},
            {"anchor_id": "R2", "side": "right", "rank": 2, "gene_symbol": r2, "product": f"conserved_anchor_R2_from_{seed_ref['accession']}"},
        ]
    else:
        consensus_anchors = [
            {"anchor_id": "L2", "side": "left", "rank": 2, "gene_symbol": "Anchor_L2", "product": "conserved_flanking_protein_L2"},
            {"anchor_id": "L1", "side": "left", "rank": 1, "gene_symbol": "Anchor_L1", "product": "conserved_flanking_protein_L1"},
            {"anchor_id": "R1", "side": "right", "rank": 1, "gene_symbol": "Anchor_R1", "product": "conserved_flanking_protein_R1"},
            {"anchor_id": "R2", "side": "right", "rank": 2, "gene_symbol": "Anchor_R2", "product": "conserved_flanking_protein_R2"},
        ]

    pd.DataFrame(consensus_anchors).to_csv(args.out_anchors_tsv, sep="\t", index=False)
    print(f"[{args.target_name}] Consensus Anchors defined and written to: {args.out_anchors_tsv}", flush=True)

    # Negative references remain provisional until anchor closure and DNA-level
    # residue exclusion have both been demonstrated.
    for neg in selected_negatives:
        neg["final_validation_status"] = "PROVISIONAL_NEGATIVE"
        neg["validation_note"] = "No validated bilateral anchor closure and DNA-level absence evidence"
        if selected_positives and neg.get("detected_count", 0) == 0:
            try:
                confirmed, note = confirm_negative_reference(
                    neg, anchor_records, args.out_seed_operon_faa,
                    args.diamond_env, args.out_dir)
                neg["validation_note"] = note
                if confirmed:
                    neg["final_validation_status"] = "CONFIRMED_NEGATIVE"
            except Exception as exc:
                neg["validation_note"] = f"Negative validation failed: {exc}"
    selected_negatives = [neg for neg in selected_negatives if neg["final_validation_status"] == "CONFIRMED_NEGATIVE"][:args.max_negative_refs]

    # Calculate positive_taxonomic_group_count
    tax_groups = set()
    for p in selected_positives:
        parts = p["taxonomy"].split(";")
        g = [x for x in parts if "g__" in x or "f__" in x or "o__" in x]
        tax_groups.add(g[0] if g else p["accession"])
        
    positive_taxonomic_group_count = len(tax_groups)

    summary_rows = []
    seen_summary_accs = set()
    for r in (selected_positives + provisional_positives + selected_negatives):
        acc = str(r["accession"])
        clean_acc = acc.replace("GB_", "").replace("RS_", "")
        if clean_acc in seen_summary_accs:
            continue
        seen_summary_accs.add(clean_acc)
        summary_rows.append({
            "target_id": args.target_name,
            "module_id": args.module_id,
            "module_name": args.module_name or args.module_id,
            "accession": clean_acc,
            "tier": r.get("tier", "Tier_C"),
            "validation_status": r["final_validation_status"],
            "validation_note": r.get("validation_note", ""),
            "detected_query_genes": r.get("detected_count", 0),
            "locus_cluster_gene_count": r.get("locus_cluster_gene_count", 0),
            "expected_query_genes": total_query_genes,
            "ANI": r.get("ANI", "NA"),
            "AF_query": r.get("AF_query", "NA"),
            "AF_ref": r.get("AF_ref", "NA"),
            "AF_min": r.get("AF_min", "NA"),
            "taxonomy": r["taxonomy"],
            "genome_fasta": r["fasta"],
            "proteome_faa": r.get("faa_path", ""),
            "original_proteome_faa": r.get("faa", ""),
            "ko_proteins": json.dumps(r.get("cluster_query_hits", {})),
            "positive_taxonomic_group_count": positive_taxonomic_group_count
        })

    if not any(row.get("validation_status") == "CONFIRMED_POSITIVE" for row in summary_rows):
        summary_rows.append({
            "target_id": args.target_name,
            "module_id": args.module_id,
            "module_name": args.module_name or args.module_id,
            "accession": "NONE",
            "tier": "NONE",
            "validation_status": "UNRESOLVED_NO_POSITIVE_REFERENCE",
            "validation_note": (f"No validated positive reference locus found for module {args.module_id}; "
                                f"scanned {len(matched_meta_rows)} local references, "
                                f"{len(positive_pool)} passed the gene-count screen, and "
                                f"{len(provisional_positives)} failed coordinate validation"),
            "detected_query_genes": 0,
            "expected_query_genes": total_query_genes,
            "ANI": "NA",
            "AF_query": "NA",
            "AF_ref": "NA",
            "AF_min": "NA",
            "taxonomy": args.target_taxonomy,
            "genome_fasta": args.target_fa,
            "proteome_faa": args.target_faa if args.target_faa else "",
            "original_proteome_faa": "",
            "ko_proteins": "{}",
            "positive_taxonomic_group_count": 0
        })

    out_df = pd.DataFrame(summary_rows)
    out_df.to_csv(args.out_refs_tsv, sep="\t", index=False)
    print(f"[{args.target_name}] Stage 0 Reference Selection Completed for {args.module_id}! Wrote {len(out_df)} validated references to {args.out_refs_tsv}", flush=True)

if __name__ == "__main__":
    main()
