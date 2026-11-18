#!/usr/bin/env python3
"""
run_pyrodigal.py
Deterministic gene calling using Pyrodigal (or biopython/prodigal fallback) with:
- Strict coordinate-based gene ordering: (contig_id, start, strand)
- Canonical ID scheme: {genome_id}___{contig_id}___gene{order:05d}
- Simultaneous output of:
  1. proteins.faa.gz
  2. genes.gff.gz
  3. protein_metadata.tsv.gz (conforming to schemas/v1/protein_metadata.schema.json)
"""

import argparse
import gzip
import os
import sys
import re
import pickle
import tempfile

def open_file(path, mode="rt"):
    # Robustly detect gzip by magic bytes \x1f\x8b or .gz extension
    is_gzip = path.endswith(".gz")
    if not is_gzip and os.path.exists(path) and os.path.getsize(path) >= 2:
        try:
            with open(path, "rb") as test_f:
                if test_f.read(2) == b"\x1f\x8b":
                    is_gzip = True
        except Exception:
            pass
    if is_gzip:
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None, errors="replace" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None, errors="replace" if "t" in mode else None)

def read_fasta_records(fpath):
    records = []
    curr_header = None
    curr_seq = []
    with open_file(fpath, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if curr_header is not None:
                    records.append((curr_header, "".join(curr_seq)))
                curr_header = line[1:].split()[0]
                curr_seq = []
            else:
                curr_seq.append(line)
        if curr_header is not None:
            records.append((curr_header, "".join(curr_seq)))
    return records

def sanitize_sequence(seq_str):
    # Remove whitespace, newlines, carriage returns, tabs
    s = re.sub(r'[\s\r\n\t]+', '', seq_str)
    # Convert non-ACGT letters to 'N' and uppercase
    s = re.sub(r'[^ACGTacgt]', 'N', s).upper()
    return s


# Header markers used by the optional ``mode=auto`` policy.  We only switch
# from metagenomic mode when every retained contig explicitly declares a
# closed/complete topology.  This keeps the production default conservative
# for ordinary draft MAGs.
_CLOSED_TOPOLOGY_PATTERNS = (
    re.compile(r"(?:^|[^a-z0-9])circular[-_=](?:yes|true|1)(?=$|[^a-z0-9])", re.I),
    re.compile(r"(?:^|[^a-z0-9])closed[-_=](?:yes|true|1)(?=$|[^a-z0-9])", re.I),
    re.compile(r"(?:^|[^a-z0-9])complete[-_=](?:yes|true|1)(?=$|[^a-z0-9])", re.I),
    re.compile(r"(?:^|[^a-z0-9])topology[=:]circular(?=$|[^a-z0-9])", re.I),
)


def infer_auto_gene_calling_mode(clean_records, complete_mode="single", fallback_mode="meta"):
    """Select a safe mode for the optional automatic gene-calling policy.

    ``clean_records`` contains FASTA identifiers after sequence filtering.  A
    complete/closed decision is made only when all retained contigs carry an
    explicit topology marker in their identifier.  The returned reason is
    written to the run log so the per-genome decision is auditable.
    """
    if complete_mode not in ("meta", "single"):
        raise ValueError(f"invalid complete_mode: {complete_mode}")
    if fallback_mode not in ("meta", "single"):
        raise ValueError(f"invalid fallback_mode: {fallback_mode}")

    if not clean_records:
        return fallback_mode, "no retained contigs"

    unmarked = []
    for contig_id, _ in clean_records:
        header = str(contig_id)
        if not any(pattern.search(header) for pattern in _CLOSED_TOPOLOGY_PATTERNS):
            unmarked.append(header)

    if not unmarked:
        return complete_mode, f"all {len(clean_records)} retained contig(s) have explicit closed/complete markers"

    preview = ", ".join(unmarked[:3])
    if len(unmarked) > 3:
        preview += f", ... (+{len(unmarked) - 3} more)"
    return fallback_mode, f"{len(unmarked)}/{len(clean_records)} contig(s) lack closed/complete markers: {preview}"


_CODON_TABLE_11 = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _scan_simple_orfs(seq, strand, trans_table):
    """Small dependency-free emergency fallback, not the production caller."""
    starts = {"ATG", "GTG", "TTG"}
    stops = {"TAA", "TAG", "TGA"}
    records = []
    for frame in range(3):
        pos = frame
        while pos + 3 <= len(seq):
            codon = seq[pos:pos + 3]
            if codon not in starts:
                pos += 3
                continue
            stop_pos = pos + 3
            while stop_pos + 3 <= len(seq) and seq[stop_pos:stop_pos + 3] not in stops:
                stop_pos += 3
            if stop_pos + 3 > len(seq):
                pos += 3
                continue
            nt_len = stop_pos + 3 - pos
            if nt_len >= 60:
                translated = "".join(
                    _CODON_TABLE_11.get(seq[i:i + 3], "X")
                    for i in range(pos, stop_pos + 3, 3)
                ).rstrip("*")
                if strand == "+":
                    begin, end = pos + 1, stop_pos + 3
                else:
                    begin, end = len(seq) - (stop_pos + 3) + 1, len(seq) - pos
                records.append((begin, end, strand, translated))
            pos = stop_pos + 3
    return records


def run_fallback_prediction(clean_records, trans_table):
    """Return deterministic emergency ORFs when Pyrodigal is unavailable."""
    if trans_table != 11:
        sys.stderr.write(
            f"[run_pyrodigal WARNING] Built-in fallback supports translation table 11; received {trans_table}.\n"
        )
    genes = []
    for contig_id, seq in clean_records:
        for begin, end, strand, translation in _scan_simple_orfs(seq, "+", trans_table):
            genes.append({
                "contig_id": contig_id,
                "start": begin,
                "end": end,
                "strand": strand,
                "trans_table": trans_table,
                "translation": translation,
                "partial": False,
            })
        rc = _reverse_complement(seq)
        for begin, end, strand, translation in _scan_simple_orfs(rc, "-", trans_table):
            genes.append({
                "contig_id": contig_id,
                "start": begin,
                "end": end,
                "strand": strand,
                "trans_table": trans_table,
                "translation": translation,
                "partial": False,
            })
    return genes, "internal_fallback", "builtin"



def _pyrodigal_worker(clean_records, mode, trans_table, genome_id, parent_conn, child_conn, closed=False, mask=False):
    try:
        parent_conn.close()
    except Exception:
        pass
    import faulthandler
    faulthandler.enable()
    tmp_path = None
    try:
        import pyrodigal
        tool_version = getattr(pyrodigal, "__version__", "3.x")
        if mode == "meta":
            orf_finder = pyrodigal.GeneFinder(meta=True, closed=closed, mask=mask)
        else:
            orf_finder = pyrodigal.GeneFinder(meta=False, closed=closed, mask=mask)
            training_sequences = [seq.encode("ascii", errors="replace") for _, seq in clean_records]
            total_length = sum(len(seq) for seq in training_sequences)
            if total_length >= 20000 and training_sequences:
                # Passing contigs separately is the Pyrodigal/Bakta training
                # API.  Pyrodigal joins them with linkers internally and does
                # not mistake an artificial raw concatenation boundary for
                # biological sequence.
                orf_finder.train(
                    training_sequences[0],
                    *training_sequences[1:],
                    translation_table=trans_table,
                )
            else:
                # Preserve the existing short-input fallback, while keeping
                # the selected edge/masking policy.
                orf_finder = pyrodigal.GeneFinder(meta=True, closed=closed, mask=mask)

        worker_genes = []
        for contig_id, seq_clean in clean_records:
            seq_bytes = seq_clean.encode("ascii", errors="replace")
            try:
                genes = orf_finder.find_genes(seq_bytes)
            except Exception as e:
                sys.stderr.write(f"[run_pyrodigal WARNING] Failed to find genes on {contig_id} in {genome_id}: {e}\n")
                continue
            for g in genes:
                strand_str = "+" if g.strand == 1 else "-"
                worker_genes.append({
                    "contig_id": contig_id,
                    "start": g.begin,
                    "end": g.end,
                    "strand": strand_str,
                    "trans_table": trans_table,
                    "translation": g.translate(),
                    "partial": g.partial_begin or g.partial_end
                })

        # Write to temporary file to avoid OS pipe buffer deadlock on large genomes (>64 KB)
        tmp_fd, tmp_path = tempfile.mkstemp(prefix=f"pyrodigal_{genome_id}_", suffix=".pkl")
        with os.fdopen(tmp_fd, "wb") as f:
            pickle.dump(worker_genes, f, protocol=pickle.HIGHEST_PROTOCOL)
        child_conn.send(("OK", tmp_path, tool_version))
    except Exception as e:
        if tmp_path and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        try:
            child_conn.send(("ERROR", str(e), "unknown"))
        except Exception:
            pass
    finally:
        try:
            child_conn.close()
        except Exception:
            pass

def run_prediction(
    genome_fasta,
    genome_id,
    out_prefix="",
    mode="meta",
    trans_table=11,
    out_faa=None,
    out_gff=None,
    out_meta=None,
    auto_complete_mode="single",
    auto_fallback_mode="meta",
):
    import multiprocessing as mp
    import faulthandler
    faulthandler.enable()

    tool_name = "pyrodigal"
    tool_version = "unknown"
    use_pyrodigal = False
    
    try:
        import pyrodigal
        tool_version = getattr(pyrodigal, "__version__", "3.x")
        use_pyrodigal = True
    except ImportError:
        pass

    raw_records = read_fasta_records(genome_fasta)
    
    # Defensive sequence sanitization:
    # 1. Strip whitespace/newlines/CR
    # 2. Convert non-ACGT to N and uppercase
    # 3. Filter out contigs < 60 bp (Prodigal dynamic programming window limit)
    # 4. Filter out contigs with no valid ACGT bases
    clean_records = []
    for contig_id, raw_seq in raw_records:
        s = sanitize_sequence(raw_seq)
        if len(s) < 60:
            continue
        if not re.search(r'[ACGT]', s):
            continue
        clean_records.append((contig_id, s))

    requested_mode = mode
    auto_reason = ""
    if mode == "auto":
        mode, auto_reason = infer_auto_gene_calling_mode(
            clean_records,
            complete_mode=auto_complete_mode,
            fallback_mode=auto_fallback_mode,
        )
        sys.stderr.write(
            f"[run_pyrodigal] auto mode selected '{mode}' for {genome_id}: {auto_reason}\n"
        )

    # A marker-driven auto decision reproduces Bakta's complete/circular CDS
    # edge policy: open ends plus masking of ambiguous sequence. Explicit
    # meta/single runs retain their historical defaults for compatibility.
    closed = False
    mask = False
    if requested_mode == "auto":
        closed = False
        mask = True

    all_genes = []

    if use_pyrodigal:
        tmp_payload_path = None
        parent_conn, child_conn = None, None
        try:
            ctx = mp.get_context()
            parent_conn, child_conn = ctx.Pipe()
            p = ctx.Process(
                target=_pyrodigal_worker,
                args=(clean_records, mode, trans_table, genome_id, parent_conn, child_conn, closed, mask),
            )
            p.start()
            child_conn.close()
            child_conn = None

            # Poll with timeout (300s) to avoid OS pipe buffer deadlocks and handle hung workers
            if parent_conn.poll(timeout=300):
                try:
                    status, res, ver = parent_conn.recv()
                except EOFError:
                    status, res, ver = "EOF", "Worker process terminated unexpectedly (EOF)", "unknown"
                p.join()

                if status == "OK":
                    tmp_payload_path = res
                    with open(tmp_payload_path, "rb") as pf:
                        all_genes = pickle.load(pf)
                    tool_name = "pyrodigal"
                    tool_version = ver
                elif status == "EOF":
                    sig_desc = f"signal {abs(p.exitcode)}" if p.exitcode < 0 else f"code {p.exitcode}"
                    sys.stderr.write(f"[run_pyrodigal WARNING] pyrodigal C-core terminated unexpectedly with {sig_desc} on {genome_id}. Switching to safe fallback engine...\n")
                    all_genes, tool_name, tool_version = run_fallback_prediction(clean_records, trans_table)
                    tool_name = "internal_fallback_after_crash"
                else:
                    sys.stderr.write(f"[run_pyrodigal WARNING] pyrodigal worker error on {genome_id}: {res}. Switching to fallback engine...\n")
                    all_genes, tool_name, tool_version = run_fallback_prediction(clean_records, trans_table)
                    tool_name = "internal_fallback_after_error"
            else:
                p.terminate()
                p.join()
                sys.stderr.write(f"[run_pyrodigal WARNING] pyrodigal timed out after 300s on {genome_id}. Switching to internal fallback engine...\n")
                all_genes, tool_name, tool_version = run_fallback_prediction(clean_records, trans_table)
                tool_name = "internal_fallback_after_timeout"
        except Exception as e:
            sys.stderr.write(f"[run_pyrodigal WARNING] pyrodigal process runner failed ({e}). Switching to fallback engine...\n")
            all_genes, tool_name, tool_version = run_fallback_prediction(clean_records, trans_table)
            tool_name = "internal_fallback"
        finally:
            if child_conn:
                try:
                    child_conn.close()
                except Exception:
                    pass
            if parent_conn:
                try:
                    parent_conn.close()
                except Exception:
                    pass
            if tmp_payload_path and os.path.exists(tmp_payload_path):
                try:
                    os.remove(tmp_payload_path)
                except OSError:
                    pass
    else:
        all_genes, tool_name, tool_version = run_fallback_prediction(clean_records, trans_table)

    # Strict deterministic coordinate sorting: (contig_id, start, strand)
    all_genes.sort(key=lambda x: (x["contig_id"], x["start"], x["strand"]))

    faa_path = out_faa or f"{out_prefix}.faa.gz"
    gff_path = out_gff or f"{out_prefix}.gff.gz"
    meta_path = out_meta or f"{out_prefix}_protein_metadata.tsv.gz"

    for p in (faa_path, gff_path, meta_path):
        d = os.path.dirname(os.path.abspath(p))
        if d:
            os.makedirs(d, exist_ok=True)

    with gzip.open(faa_path, "wt", encoding="utf-8") as faa_out, \
         gzip.open(gff_path, "wt", encoding="utf-8") as gff_out, \
         gzip.open(meta_path, "wt", encoding="utf-8") as meta_out:

        # Header for GFF3
        gff_out.write("##gff-version 3\n")
        gff_out.write(f"##genome-id {genome_id}\n")
        gff_out.write(f"##gene-calling-requested-mode {requested_mode}\n")
        gff_out.write(f"##gene-calling-resolved-mode {mode}\n")
        if auto_reason:
            gff_out.write(f"##gene-calling-auto-reason {auto_reason}\n")

        # Header for protein metadata
        meta_out.write("genome_id\tcontig_id\tgene_id\tprotein_id\taa_length\tstart\tend\tstrand\ttranslation_table\tgene_calling_tool\tgene_calling_version\tgene_calling_mode\n")

        order_by_contig = {}
        for gene in all_genes:
            c_id = gene["contig_id"]
            order_by_contig[c_id] = order_by_contig.get(c_id, 0) + 1
            order = order_by_contig[c_id]

            gene_id = f"{genome_id}___{c_id}___gene{order:05d}"
            protein_id = f"{genome_id}___{c_id}___prot{order:05d}"
            aa_seq = gene["translation"]
            aa_len = len(aa_seq)

            # Write FAA
            faa_out.write(f">{protein_id} {gene_id} {genome_id} [{gene['start']} - {gene['end']}] [strand={gene['strand']}]\n")
            # Wrap fasta at 60 characters
            for i in range(0, len(aa_seq), 60):
                faa_out.write(aa_seq[i:i+60] + "\n")

            # Write GFF3
            gff_out.write(f"{c_id}\t{tool_name}\tCDS\t{gene['start']}\t{gene['end']}\t.\t{gene['strand']}\t0\tID={gene_id};protein_id={protein_id};transl_table={gene['trans_table']}\n")

            # Write Metadata
            meta_out.write(f"{genome_id}\t{c_id}\t{gene_id}\t{protein_id}\t{aa_len}\t{gene['start']}\t{gene['end']}\t{gene['strand']}\t{gene['trans_table']}\t{tool_name}\t{tool_version}\t{mode}\n")

    return faa_path, gff_path, meta_path

def main():
    parser = argparse.ArgumentParser(description="Deterministic gene calling with GFF and metadata extraction.")
    parser.add_argument("-i", "--genome_fasta", required=True, help="Path to input genome FASTA (.fna/.fa/.fasta, can be .gz)")
    parser.add_argument("-g", "--genome_id", required=True, help="Identifier for the genome")
    parser.add_argument("-o", "--out_prefix", default="", help="Prefix for output files (.faa.gz, .gff.gz, etc.)")
    parser.add_argument("--out_faa", default="", help="Direct destination path for output .faa.gz")
    parser.add_argument("--out_gff", default="", help="Direct destination path for output .gff.gz")
    parser.add_argument("--out_meta", default="", help="Direct destination path for output _protein_metadata.tsv.gz")
    parser.add_argument("-m", "--mode", default="meta", choices=["meta", "single", "auto"], help="Pyrodigal mode (default: meta; auto selects single only for explicitly closed/complete FASTA headers)")
    parser.add_argument("-t", "--trans_table", type=int, default=11, help="Translation table (default: 11 for Archaea/Bacteria)")
    parser.add_argument("--auto-complete-mode", choices=["meta", "single"], default="single", help="Mode used by --mode auto when all retained contigs have closed/complete markers (default: single)")
    parser.add_argument("--auto-fallback-mode", choices=["meta", "single"], default="meta", help="Mode used by --mode auto when topology markers are missing (default: meta)")

    args = parser.parse_args()
    if not args.out_prefix and not (args.out_faa and args.out_gff and args.out_meta):
        parser.error("Either -o/--out_prefix or all of (--out_faa, --out_gff, --out_meta) must be specified.")

    faa, gff, meta = run_prediction(
        args.genome_fasta,
        args.genome_id,
        args.out_prefix,
        args.mode,
        args.trans_table,
        out_faa=args.out_faa,
        out_gff=args.out_gff,
        out_meta=args.out_meta,
        auto_complete_mode=args.auto_complete_mode,
        auto_fallback_mode=args.auto_fallback_mode,
    )
    print(f"[run_pyrodigal] Completed gene calling for {args.genome_id}:")
    print(f"  - Proteins: {faa}")
    print(f"  - GFF:      {gff}")
    print(f"  - Metadata: {meta}")

if __name__ == "__main__":
    main()
