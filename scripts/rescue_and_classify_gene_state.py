#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import re
import gzip
import pandas as pd
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    class DummyRecord:
        def __init__(self, record_id, seq):
            self.id = record_id
            self.seq = seq
    class SeqIO:
        @staticmethod
        def parse(handle_or_file, fmt="fasta"):
            records = []
            cur_id = None
            cur_seq = []
            if isinstance(handle_or_file, str):
                f = open(handle_or_file, "r", errors="replace")
                should_close = True
            else:
                f = handle_or_file
                should_close = False
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if cur_id is not None:
                        records.append(DummyRecord(cur_id, "".join(cur_seq)))
                    cur_id = line[1:].split()[0]
                    cur_seq = []
                else:
                    cur_seq.append(line)
            if cur_id is not None:
                records.append(DummyRecord(cur_id, "".join(cur_seq)))
            if should_close:
                f.close()
            return records

        @staticmethod
        def write(record, handle, fmt="fasta"):
            handle.write(f">{record.id}\n{str(record.seq)}\n")

MODULE_TARGET_KOS = {
    "M00595": {"K17224", "K17222", "K17223", "K17227", "K17226", "K17225", "K22622"},
    "M00159": {"K02117", "K02118", "K02119", "K02120", "K02121", "K02122", "K02107", "K02123", "K02124"},
    "M00157": {"K02111", "K02112", "K02113", "K02114", "K02115", "K02108", "K02109", "K02110"},
    "M00158": {"K02126", "K02127", "K02128", "K02129", "K02130", "K02131", "K02132", "K02133", "K02134", "K02135", "K02136", "K02137", "K02138"},
    "M00160": {"K02145", "K02146", "K02147", "K02148", "K02149", "K02150", "K02151", "K02152", "K02153", "K02154", "K02155"},
}

ALT_MODULE_MAP = {
    "M00159": "M00157",
    "M00157": "M00159",
    "M00158": "M00160",
    "M00160": "M00158",
}

def parse_args():
    parser = argparse.ArgumentParser(description="Stage 2: Multi-layer sequence rescue with frameshift and archaeal recoding guards, gene-level state classification, OFI metrics and cluster aggregation.")
    parser.add_argument("--target_name", required=True, help="Target genome sample ID")
    parser.add_argument("--target_locus_fa", required=True, help="Target locus DNA FASTA")
    parser.add_argument("--target_fa", default="", help="Whole target genome FASTA for global DNA-level absence search")
    parser.add_argument("--target_faa", required=True, help="Target proteome FAA")
    parser.add_argument("--query_operon_faa", required=True, help="Query operon / gene cluster proteins FAA")
    parser.add_argument("--topology_tsv", default="", help="Optional Stage 1 anchor topology TSV")
    parser.add_argument("--target_annotation_tsv", default="", help="Optional gene annotation TSV (with KO or gene product)")
    parser.add_argument("--target_gff", default="", help="Optional target GFF file")
    parser.add_argument("--diamond_env", default="", help="Conda environment for diamond")
    parser.add_argument("--threads", type=int, default=4, help="Number of CPU threads")
    parser.add_argument("--module_id", default="M00159", help="ATP synthase module ID")
    parser.add_argument("--module_name", default="", help="ATP synthase module human-readable name")
    parser.add_argument("--blastp_identity", type=float, default=25.0, help="Minimum identity for Level 2 protein match")
    parser.add_argument("--blastp_coverage", type=float, default=30.0, help="Minimum coverage for Level 2 protein match")
    parser.add_argument("--blastx_evalue", type=float, default=1e-3, help="E-value cutoff for Level 3 DNA search")
    parser.add_argument("--max_operon_intergenic_bp", type=int, default=300, help="Largest adjacent CDS gap for a candidate operon")
    parser.add_argument("--alternative_module", default="", help="Alternative/isoform module ID (e.g. M00157)")
    parser.add_argument("--alternative_kos", default="", help="Comma-separated list of known alternative/isoform KOs")
    parser.add_argument("--out_gene_tsv", required=True, help="Output TSV with per-gene detailed rescue status")
    parser.add_argument("--out_cluster_tsv", required=True, help="Output TSV with aggregated cluster metrics")
    return parser.parse_args()

def run_cmd(cmd):
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed ({res.returncode}): {cmd}\n{res.stderr}")
    return res

def run_locus_blastx_compatible(modern_command, legacy_command):
    """Use translated-query output when DIAMOND supports it (v2.0.7+)."""
    try:
        run_cmd(modern_command)
        return True
    except RuntimeError as exc:
        message = str(exc)
        if "qseq_translated" not in message or "Invalid output field" not in message:
            raise
        print("DIAMOND lacks qseq_translated (available since v2.0.7); retrying without translated-query diagnostics.", flush=True)
        run_cmd(legacy_command)
        return False

def load_gene_coordinates(gff_file, faa_file):
    """Return real CDS coordinates; never infer a locus from a protein hit alone."""
    coords = {}
    if gff_file and os.path.exists(gff_file):
        opener = gzip.open if gff_file.endswith(".gz") else open
        with opener(gff_file, "rt", errors="replace") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "CDS":
                    continue
                attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
                entry = (fields[0], int(fields[3]) - 1, int(fields[4]), fields[6])
                for key in ("protein_id", "ID", "Name", "locus_tag"):
                    if attrs.get(key):
                        coords[attrs[key]] = entry
    if faa_file and os.path.exists(faa_file):
        opener = gzip.open if faa_file.endswith(".gz") else open
        with opener(faa_file, "rt", errors="replace") as faa_handle:
            for record in SeqIO.parse(faa_handle, "fasta"):
                if record.id in coords:
                    continue
                parts = record.description.split("#")
                if len(parts) >= 4:
                    try:
                        start, end, strand = int(parts[1]), int(parts[2]), int(parts[3])
                        contig = record.id.rsplit("___", 1)[0] if "___" in record.id else record.id.rsplit("_", 1)[0]
                        coords[record.id] = (contig, start - 1, end, "+" if strand >= 0 else "-")
                    except ValueError:
                        pass
    return coords

def hit_in_interval(hit_id, coords, target_contig, interval_start, interval_end):
    hit = coords.get(hit_id)
    if hit is None:
        return None
    contig, start, end, _ = hit
    if contig != target_contig and contig.split("___")[-1] != target_contig.split("___")[-1]:
        return False
    return start < interval_end and end > interval_start

def classify_operon_structure(gene_results, coords, query_order, intervening_genes, max_gap=300):
    """Predict an operon from observed CDS order, orientation and spacing."""
    result = {"operon_structure_state": "UNRESOLVED", "operon_prediction": "UNRESOLVED",
              "same_strand": "NA", "order_consistent": "NA", "max_intergenic_bp": "NA",
              "intervening_non_target_count": "NA"}
    by_query = {gene["query_gene"]: gene for gene in gene_results}
    if any(by_query[q]["gene_state"] not in ("PRESENT", "ANNOTATION_MISS") for q in query_order):
        result["operon_structure_state"] = "INCOMPLETE_CLUSTER"
        return result
    pids = [by_query[q]["hit_id"] for q in query_order]
    if len(set(pids)) != len(pids) or any(pid not in coords for pid in pids):
        return result
    locations = [coords[pid] for pid in pids]
    if len(set(locations)) != len(locations):
        result.update(operon_structure_state="DUPLICATED_HIT", operon_prediction="NOT_SUPPORTED")
        return result
    if len({loc[0] for loc in locations}) != 1:
        result.update(operon_structure_state="DISPERSED", operon_prediction="NOT_SUPPORTED")
        return result
    same_strand = len({loc[3] for loc in locations}) == 1
    result["same_strand"] = "YES" if same_strand else "NO"
    if not same_strand:
        result.update(operon_structure_state="DISRUPTED_ORIENTATION", operon_prediction="NOT_SUPPORTED")
        return result
    starts = [loc[1] for loc in locations]
    ordered = starts == sorted(starts) or starts == sorted(starts, reverse=True)
    result["order_consistent"] = "YES" if ordered else "NO"
    if not ordered:
        result.update(operon_structure_state="DISRUPTED_ORDER", operon_prediction="NOT_SUPPORTED")
        return result
    physical = sorted(locations, key=lambda loc: loc[1])
    result["max_intergenic_bp"] = max((right[1] - left[2] for left, right in zip(physical, physical[1:])), default=0)
    target_locations = set(locations)
    other_genes = set(intervening_genes) - set(pids)
    contig = physical[0][0]
    for pid, location in coords.items():
        if pid not in pids and location not in target_locations and location[0] == contig \
                and location[1] < physical[-1][2] and location[2] > physical[0][1]:
            other_genes.add(pid)
    unrelated = len(other_genes)
    result["intervening_non_target_count"] = unrelated
    if unrelated:
        result["operon_structure_state"] = "INTERRUPTED_CLUSTER"
    elif result["max_intergenic_bp"] > max_gap:
        result["operon_structure_state"] = "SPACED_CLUSTER"
    else:
        result.update(operon_structure_state="CANDIDATE_OPERON", operon_prediction="CANDIDATE_OPERON")
    return result

def safe_int(val, default=0):
    if pd.isna(val) or val is None or val == "":
        return default
    try:
        return int(val)
    except Exception:
        return default

def safe_str(val, default=""):
    if pd.isna(val) or val is None:
        return default
    return str(val)

def load_gene_annotations(annotation_tsv, gff_file):
    """Parses optional annotation TSV or GFF to extract KO or product annotations."""
    gene_info = {}
    if annotation_tsv and os.path.exists(annotation_tsv) and os.path.getsize(annotation_tsv) > 0:
        try:
            df_ann = pd.read_csv(annotation_tsv, sep="\t")
            # Try to identify protein ID column and KO / function columns
            pid_col = None
            for c in ["protein_id", "gene_id", "query", "CDS_ID", "ID"]:
                if c in df_ann.columns:
                    pid_col = c
                    break
            ko_col = None
            for c in ["KO", "ko", "ko_id", "KEGG_ko"]:
                if c in df_ann.columns:
                    ko_col = c
                    break
            prod_col = None
            for c in ["product", "description", "function", "Preferred_name", "ko_definition"]:
                if c in df_ann.columns:
                    prod_col = c
                    break

            if pid_col:
                for _, r in df_ann.iterrows():
                    pid = str(r[pid_col])
                    ko = str(r[ko_col]) if ko_col and pd.notna(r[ko_col]) else ""
                    prod = str(r[prod_col]) if prod_col and pd.notna(r[prod_col]) else ""
                    gene_info[pid] = {"KO": ko, "product": prod}
        except Exception as e:
            print(f"Notice: Failed parsing annotation TSV {annotation_tsv}: {e}", file=sys.stderr)

    if not gene_info and gff_file and os.path.exists(gff_file) and os.path.getsize(gff_file) > 0:
        try:
            opener = gzip.open if gff_file.endswith(".gz") else open
            with opener(gff_file, "rt", errors="replace") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 9 and parts[2] == "CDS":
                        attr = parts[8]
                        m_id = re.search(r'ID=([^;]+)', attr)
                        pid = m_id.group(1) if m_id else ""
                        m_prod = re.search(r'product=([^;]+)', attr)
                        prod = m_prod.group(1) if m_prod else ""
                        m_ko = re.search(r'ko=([^;]+)', attr)
                        ko = m_ko.group(1) if m_ko else ""
                        if pid:
                            gene_info[pid] = {"KO": ko, "product": prod}
        except Exception as e:
            print(f"Notice: Failed parsing GFF {gff_file}: {e}", file=sys.stderr)

    return gene_info

def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.out_gene_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_cluster_tsv), exist_ok=True)

    cmd_prefix = f"conda run -p {args.diamond_env} " if args.diamond_env and os.path.isabs(args.diamond_env) else (f"conda run -n {args.diamond_env} " if args.diamond_env else "")

    # 1. Parse query operon
    if os.path.exists(args.query_operon_faa):
        with open(args.query_operon_faa, "r", errors="replace") as query_handle:
            query_records = list(SeqIO.parse(query_handle, "fasta"))
    else:
        query_records = []
    query_dict = {r.id: r for r in query_records}
    total_query_genes = len(query_records)
    if total_query_genes == 0:
        print(f"[{args.target_name}] Seed query operon FASTA is empty; state is unresolved.", flush=True)
        pd.DataFrame(columns=[
            "target_id", "query_gene_id", "gene_symbol", "product", "gene_state",
            "matched_protein_id", "identity", "coverage", "frameshifts", "premature_stops",
            "mutation_evidence_status", "note"
        ]).to_csv(args.out_gene_tsv, sep="\t", index=False)

        cluster_row = {
            "target_id": args.target_name,
            "expected_gene_count": 0,
            "protein_detected_count": 0,
            "location_unknown_count": 0,
            "blastx_rescued_count": 0,
            "intact_rescued_count": 0,
            "pseudogene_count": 0,
            "frameshift_count": 0,
            "premature_stop_count": 0,
            "mutation_evidence_status": "NOT_ASSESSED",
            "annotation_issue": "NO",
            "pseudogene_detected": "NO",
            "member_gene_states": "",
            "cluster_state": "UNRESOLVED_NO_QUERY",
            "T_local": 0.0,
            "T_global": 0.0,
            "R_known": 0.0,
            "I_novel": 0.0,
            "A_local": "NA",
            "A_global": "NA"
        }
        pd.DataFrame([cluster_row]).to_csv(args.out_cluster_tsv, sep="\t", index=False)
        return

    # 2. Parse Stage 1 topology TSV if available
    topo_path = args.topology_tsv
    if not topo_path or not os.path.exists(topo_path):
        cand_topo = os.path.join(os.path.dirname(args.target_locus_fa), "anchor_topology.tsv")
        if os.path.exists(cand_topo):
            topo_path = cand_topo

    target_slice_start = 0
    target_slice_end = 0
    ancestral_interval_start = 0
    ancestral_interval_end = 0
    target_contig = ""
    intervening_genes_list = []
    same_contig = "YES"
    syntenic_geometry = "BILATERAL_SAME_CONTIG"

    if topo_path and os.path.exists(topo_path) and os.path.getsize(topo_path) > 0:
        try:
            topo_df = pd.read_csv(topo_path, sep="\t")
            if not topo_df.empty:
                topo_row = topo_df.iloc[0].to_dict()
                target_slice_start = safe_int(topo_row.get("target_slice_start", 0))
                target_slice_end = safe_int(topo_row.get("target_slice_end", 0))
                ancestral_interval_start = safe_int(topo_row.get("ancestral_interval_start", 0))
                ancestral_interval_end = safe_int(topo_row.get("ancestral_interval_end", 0))
                target_contig = safe_str(topo_row.get("target_contig", ""))
                same_contig = safe_str(topo_row.get("same_contig", "YES"))
                syntenic_geometry = safe_str(topo_row.get("syntenic_geometry", "BILATERAL_SAME_CONTIG"))
                interv_raw = safe_str(topo_row.get("intervening_genes", ""))
                if interv_raw and interv_raw != "CROSS_CONTIG_SPLIT":
                    intervening_genes_list = [g.strip() for g in interv_raw.split(",") if g.strip()]
        except Exception as e:
            print(f"Warning parsing topology TSV {topo_path}: {e}", file=sys.stderr)

    intervening_genes_set = set(intervening_genes_list)
    has_bounded_interval = (ancestral_interval_start >= 0 and ancestral_interval_end > ancestral_interval_start)
    gene_coordinates = load_gene_coordinates(args.target_gff, args.target_faa)

    # 3. Load gene annotations for replacement detection
    target_annot = args.target_annotation_tsv
    if not target_annot or not os.path.exists(target_annot):
        cand_annot = os.path.join(os.path.dirname(args.target_locus_fa), f"{args.target_name}_100kb_locus_proteins.tsv")
        if os.path.exists(cand_annot):
            target_annot = cand_annot
    gene_annotations = load_gene_annotations(target_annot, args.target_gff)

    # --------------------------------------------------------------------------
    # Level 2: DIAMOND blastp (Query vs Target Proteome)
    # --------------------------------------------------------------------------
    l2_hits = {}  # q_id -> best non-overlapping protein match
    target_faa_searched = False
    if os.path.exists(args.target_faa) and os.path.getsize(args.target_faa) > 0:
        target_faa_searched = True
        tmp_dmnd = args.out_gene_tsv + ".l2.dmnd"
        tmp_blastp = args.out_gene_tsv + ".l2.tsv"
        run_cmd(f"{cmd_prefix}diamond makedb --in {args.target_faa} -d {tmp_dmnd} --quiet")
        run_cmd(f"{cmd_prefix}diamond blastp -q {args.query_operon_faa} -d {tmp_dmnd} -o {tmp_blastp} -e 1e-4 --id {args.blastp_identity} --query-cover {args.blastp_coverage} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --quiet")

        if os.path.exists(tmp_blastp) and os.path.getsize(tmp_blastp) > 0:
            df_l2 = pd.read_csv(tmp_blastp, sep="\t", header=None)
            for _, row in df_l2.iterrows():
                qid = str(row[0])
                sid = str(row[1])
                if qid not in l2_hits and sid not in {h["sseqid"] for h in l2_hits.values()}:
                    cov = float(row[3]) / float(row[12]) * 100.0 if float(row[12]) > 0 else 0.0
                    is_loc = hit_in_interval(sid, gene_coordinates, target_contig,
                                             ancestral_interval_start, ancestral_interval_end) if has_bounded_interval and target_contig else None
                    l2_hits[qid] = {
                        "sseqid": sid,
                        "pident": float(row[2]),
                        "coverage": cov,
                        "evalue": float(row[10]),
                        "is_local": is_loc
                    }
        # cleanup
        if os.path.exists(tmp_dmnd + ".dmnd"): os.remove(tmp_dmnd + ".dmnd")
        if os.path.exists(tmp_blastp): os.remove(tmp_blastp)

    # --------------------------------------------------------------------------
    # Level 3A & 3B: Frameshift-aware blastx against Target Locus DNA
    # [P0-1] Bounded strictly to ancestral_interval
    # --------------------------------------------------------------------------
    l3_hits = {}  # q_id -> {found, pident, coverage, evalue, frameshifts, premature_stops, note, state}
    missing_queries = [qid for qid in query_dict.keys() if (qid not in l2_hits or not l2_hits[qid]["is_local"])]
    target_locus_searched = False
    translated_diagnostics = None

    if missing_queries and os.path.exists(args.target_locus_fa) and os.path.getsize(args.target_locus_fa) > 0:
        target_locus_searched = False
        tmp_q_faa = args.out_gene_tsv + ".missing_queries.faa"
        with open(tmp_q_faa, "w") as f_mq:
            for qid in missing_queries:
                SeqIO.write(query_dict[qid], f_mq, "fasta")

        tmp_dmnd_dna = args.out_gene_tsv + ".locus.dmnd"
        tmp_blastx = args.out_gene_tsv + ".l3_blastx.tsv"

        # blastx translates the DNA query against a protein database.
        run_cmd(f"{cmd_prefix}diamond makedb --in {tmp_q_faa} -d {tmp_dmnd_dna} --quiet")
        blastx_base = (f"{cmd_prefix}diamond blastx -q {args.target_locus_fa} -d {tmp_dmnd_dna} -o {tmp_blastx} "
                       f"-e {args.blastx_evalue} -F 15 --range-culling --outfmt 6 "
                       "qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
                       "evalue bitscore qlen slen ")
        translated_diagnostics = run_locus_blastx_compatible(
            blastx_base + "qseq_translated sseq --quiet",
            blastx_base + "sseq --quiet")
        target_locus_searched = True

        if os.path.exists(tmp_blastx) and os.path.getsize(tmp_blastx) > 0:
            try:
                df_l3 = pd.read_csv(tmp_blastx, sep="\t", header=None)
                for _, row in df_l3.iterrows():
                    qid = str(row[1])

                    # Local nucleotide coordinates on target_locus_fa (1-based)
                    sstart = int(row[6])
                    send = int(row[7])
                    s_min = min(sstart, send)
                    s_max = max(sstart, send)

                    # Genomic coordinates on target contig
                    genomic_s = target_slice_start + s_min - 1
                    genomic_e = target_slice_start + s_max

                    # [P0-1] Ancestral Interval Coordinate Check:
                    # Filter out any hits outside [ancestral_interval_start, ancestral_interval_end]
                    if has_bounded_interval:
                        overlap = (genomic_s < ancestral_interval_end and genomic_e > ancestral_interval_start)
                        if not overlap:
                            continue

                    qlen = float(row[13])
                    aln_len = float(row[3])
                    cov = (aln_len / qlen * 100.0) if qlen > 0 else 0.0
                    pident = float(row[2])
                    evalue = float(row[10])
                    sseq = str(row[14]) if translated_diagnostics and len(row) > 15 else ""

                    # Detect frameshifts: Diamond outputs '\' or '/' for frameshifts
                    frameshifts = sseq.count("\\") + sseq.count("/")
                    # Detect premature stops: '*' in translated alignment before the end
                    stops = sseq[:-1].count("*")

                    # A stop cannot be dismissed as Sec/Pyl without codon and machinery evidence.
                    recoding_flag = False

                    # [Core-3] Restore ANNOTATION_MISS for intact ORF without disruptive mutations
                    if translated_diagnostics is False:
                        gene_state = "AMBIGUOUS"
                        note = (f"Homologous DNA match detected (id={pident:.1f}%, cov={cov:.1f}%); "
                                "frameshift/stop evidence unavailable with this DIAMOND version")
                    elif cov >= 60.0 and frameshifts == 0 and stops == 0:
                        gene_state = "ANNOTATION_MISS"
                        note = f"Unannotated intact ORF rescued (id={pident:.1f}%, cov={cov:.1f}%)"
                    elif frameshifts > 0 or stops > 0:
                        if frameshifts >= 2 or stops >= 2 or cov >= 70.0:
                            gene_state = "PSEUDOGENIZED_HIGH_CONF"
                        else:
                            gene_state = "PSEUDOGENIZED_PUTATIVE"
                        note = f"Disruptive mutation (frameshifts={frameshifts}, stops={stops}, cov={cov:.1f}%)"
                    elif cov < 30.0:
                        gene_state = "NO_DETECTABLE_HOMOLOGOUS_RESIDUE"
                        note = f"Negligible fragment below cutoff (cov={cov:.1f}%)"
                    else:
                        gene_state = "AMBIGUOUS"
                        note = f"Weak match (id={pident:.1f}%, cov={cov:.1f}%)"

                    l3_hits[qid] = {
                        "found": True,
                        "pident": pident,
                        "coverage": cov,
                        "evalue": evalue,
                        "frameshifts": frameshifts,
                        "premature_stops": stops,
                        "mutation_evidence_status": "ASSESSED" if translated_diagnostics else "VERSION_LIMITED",
                        "bitscore": float(row[11]),
                        "recoding_detected": recoding_flag,
                        "genomic_start": genomic_s,
                        "genomic_end": genomic_e,
                        "strand": "+" if sstart <= send else "-",
                        "state": gene_state,
                        "note": note
                    } if qid not in l3_hits or float(row[11]) > l3_hits[qid]["bitscore"] else l3_hits[qid]
            except Exception as e:
                raise RuntimeError(f"Cannot interpret locus blastx evidence: {e}") from e

        # Cleanup
        if os.path.exists(tmp_q_faa): os.remove(tmp_q_faa)
        if os.path.exists(tmp_dmnd_dna + ".dmnd"): os.remove(tmp_dmnd_dna + ".dmnd")
        if os.path.exists(tmp_blastx): os.remove(tmp_blastx)
    elif os.path.exists(args.target_locus_fa) and os.path.getsize(args.target_locus_fa) > 0:
        target_locus_searched = True

    # --------------------------------------------------------------------------
    # Gene-Level Classification
    # --------------------------------------------------------------------------
    gene_results = []
    member_states = []

    total_frameshifts = 0
    total_premature_stops = 0

    for qid in query_dict.keys():
        if qid in l2_hits and l2_hits[qid]["is_local"]:
            state = "PRESENT"
            note = f"Protein detected locally ({l2_hits[qid]['sseqid']}, id={l2_hits[qid]['pident']:.1f}%, cov={l2_hits[qid]['coverage']:.1f}%)"
            pid_hit = l2_hits[qid]["sseqid"]
            pident = l2_hits[qid]["pident"]
            cov = l2_hits[qid]["coverage"]
            fs = 0
            stops = 0
        elif qid in l3_hits and l3_hits[qid]["state"] in ["ANNOTATION_MISS", "PSEUDOGENIZED_HIGH_CONF", "PSEUDOGENIZED_PUTATIVE", "AMBIGUOUS"]:
            hit = l3_hits[qid]
            state = hit["state"]
            fs = hit["frameshifts"]
            stops = hit["premature_stops"]
            total_frameshifts += fs
            total_premature_stops += stops
            pid_hit = f"DNA_ORF_{qid}"
            if target_contig:
                gene_coordinates[pid_hit] = (target_contig, hit["genomic_start"],
                                             hit["genomic_end"], hit["strand"])
            pident = hit["pident"]
            cov = hit["coverage"]
            note = hit["note"]
        elif qid in l2_hits:
            # Hit exists in genome, but not in ancestral interval
            state = "PRESENT_GLOBAL_ONLY" if l2_hits[qid]["is_local"] is False else "PRESENT_LOCATION_UNKNOWN"
            note = f"Protein detected outside locus or coordinates unknown ({l2_hits[qid]['sseqid']}, id={l2_hits[qid]['pident']:.1f}%, cov={l2_hits[qid]['coverage']:.1f}%)"
            pid_hit = l2_hits[qid]["sseqid"]
            pident = l2_hits[qid]["pident"]
            cov = l2_hits[qid]["coverage"]
            fs = 0
            stops = 0
        else:
            state = "NO_DETECTABLE_HOMOLOGOUS_RESIDUE"
            note = "No significant protein or DNA homology"
            pid_hit = "-"
            pident = 0.0
            cov = 0.0
            fs = 0
            stops = 0

        member_states.append(f"{qid}:{state}")
        gene_results.append({
            "target_id": args.target_name,
            "query_gene": qid,
            "gene_state": state,
            "hit_id": pid_hit,
            "identity": pident,
            "coverage": cov,
            "frameshifts": fs,
            "premature_stops": stops,
            "mutation_evidence_status": l3_hits.get(qid, {}).get("mutation_evidence_status", "NOT_ASSESSED"),
            "note": note
        })

    # Save detailed gene-level table
    pd.DataFrame(gene_results).to_csv(args.out_gene_tsv, sep="\t", index=False)

    # Search the complete assembly before claiming genome-wide absence. A
    # missing protein annotation or local locus hit is insufficient evidence.
    genome_dna_searched = False
    genome_dna_hit = False
    if args.target_fa and os.path.isfile(args.target_fa) and os.path.getsize(args.target_fa):
        unresolved_queries = [query_dict[g["query_gene"]] for g in gene_results
                              if g["gene_state"] in ("NO_DETECTABLE_HOMOLOGOUS_RESIDUE", "AMBIGUOUS")]
        if unresolved_queries:
            genome_query = args.out_gene_tsv + ".global_queries.faa"
            genome_db = args.out_gene_tsv + ".global.dmnd"
            genome_hits = args.out_gene_tsv + ".global_blastx.tsv"
            with open(genome_query, "w") as genome_query_handle:
                for query_record in unresolved_queries:
                    SeqIO.write(query_record, genome_query_handle, "fasta")
            run_cmd(f"{cmd_prefix}diamond makedb --in {genome_query} -d {genome_db} --quiet")
            run_cmd(f"{cmd_prefix}diamond blastx -q {args.target_fa} -d {genome_db} -o {genome_hits} -e {args.blastx_evalue} -F 15 --subject-cover 30 --quiet")
            genome_dna_searched = True
            genome_dna_hit = os.path.getsize(genome_hits) > 0
            for path in (genome_query, genome_db + ".dmnd", genome_hits):
                if os.path.exists(path):
                    os.remove(path)
        else:
            genome_dna_searched = True

    # --------------------------------------------------------------------------
    # Quantitative OFI Metrics Calculation
    # --------------------------------------------------------------------------
    state_set = set(g["gene_state"] for g in gene_results)

    local_present_count = sum(1 for g in gene_results if g["gene_state"] == "PRESENT")
    intact_rescued_count = sum(1 for g in gene_results if g["gene_state"] == "ANNOTATION_MISS")
    pseudo_count = sum(1 for g in gene_results if "PSEUDOGENIZED" in g["gene_state"])
    local_dna_residue_count = sum(1 for hit in l3_hits.values()
                                  if hit["state"] != "NO_DETECTABLE_HOMOLOGOUS_RESIDUE")
    global_only_count = sum(1 for g in gene_results if g["gene_state"] == "PRESENT_GLOBAL_ONLY")
    location_unknown_count = sum(1 for g in gene_results if g["gene_state"] == "PRESENT_LOCATION_UNKNOWN")
    absent_count = sum(1 for g in gene_results if g["gene_state"] == "NO_DETECTABLE_HOMOLOGOUS_RESIDUE")

    total_protein_detected = local_present_count + global_only_count + location_unknown_count
    total_blastx_rescued = intact_rescued_count + pseudo_count

    # T_local: fraction of target subunits retained locally (including intact rescued ORFs)
    T_local = round((local_present_count + intact_rescued_count) / total_query_genes, 4) if total_query_genes > 0 else 0.0

    # T_global: fraction of target subunits retained across whole genome
    T_global = round((total_protein_detected + intact_rescued_count) / total_query_genes, 4) if total_query_genes > 0 else 0.0

    # R_known & I_novel: Known Replacement and Novel Insertion Occupancy Rates
    alt_kos = set()
    if getattr(args, "alternative_kos", None) and args.alternative_kos.strip():
        alt_kos = {k.strip() for k in args.alternative_kos.split(",") if k.strip()}
    elif getattr(args, "alternative_module", None) and args.alternative_module.strip():
        alt_mod = args.alternative_module.strip()
        alt_kos = MODULE_TARGET_KOS.get(alt_mod, set())
    else:
        alt_mod = ALT_MODULE_MAP.get(args.module_id, "")
        alt_kos = MODULE_TARGET_KOS.get(alt_mod, set())
    n_alt_expected = len(alt_kos) if alt_kos else 8

    # Scan intervening genes for known alternative module subunits
    n_alt_matched = 0
    matched_alt_pids = set()
    matched_target_pids = set(g["hit_id"] for g in gene_results if g["gene_state"] == "PRESENT")

    for pid in intervening_genes_list:
        if pid in matched_target_pids:
            continue
        info = gene_annotations.get(pid, {})
        ko = info.get("KO", "")
        prod = info.get("product", "").lower()
        is_alt = False
        if ko and ko in alt_kos:
            is_alt = True
        if is_alt:
            n_alt_matched += 1
            matched_alt_pids.add(pid)

    R_known = min(1.0, round(n_alt_matched / n_alt_expected, 4)) if n_alt_expected > 0 else 0.0

    # Novel insertion: intervening genes that are neither target nor alternative module
    n_intervening = len(intervening_genes_list)
    n_novel = max(0, n_intervening - len(matched_target_pids.intersection(set(intervening_genes_list))) - n_alt_matched)
    I_novel = min(1.0, round(n_novel / total_query_genes, 4)) if total_query_genes > 0 else 0.0
    if n_intervening == 0:
        R_known = 0.0
        I_novel = 0.0
    replacement_search_status = "COMPLETE" if all(pid in gene_annotations for pid in intervening_genes_list) else "INCOMPLETE"

    # Tri-State Evidence A_local and A_global (Core-7)
    local_search_complete = target_locus_searched and has_bounded_interval and same_contig == "YES" and not any(g["gene_state"] == "PRESENT_LOCATION_UNKNOWN" for g in gene_results)
    if local_present_count + intact_rescued_count + pseudo_count + local_dna_residue_count > 0:
        A_local = "0"
    elif local_search_complete and absent_count == total_query_genes:
        A_local = "1"
    else:
        A_local = "NA"  # Local negative search incomplete / not run

    if total_protein_detected + intact_rescued_count + pseudo_count + local_dna_residue_count > 0 or genome_dna_hit:
        A_global = "0"
    elif target_faa_searched and genome_dna_searched and absent_count == total_query_genes:
        A_global = "1"
    else:
        A_global = "NA"  # Global proteome search incomplete / not run

    # Flags
    annotation_issue = "YES" if intact_rescued_count > 0 else "NO"
    pseudogene_detected = "YES" if (pseudo_count > 0 or total_frameshifts > 0 or total_premature_stops > 0) else "NO"

    # Cluster State logic
    if local_present_count == total_query_genes:
        cluster_state = "COMPLETE"
    elif (local_present_count + intact_rescued_count) == total_query_genes:
        cluster_state = "COMPLETE_AFTER_RESCUE"
    elif absent_count == total_query_genes:
        cluster_state = "FULLY_ABSENT"
    elif local_present_count > 0 and absent_count > 0 and pseudo_count == 0 and intact_rescued_count == 0:
        cluster_state = "PARTIAL_ABSENT"
    elif local_present_count > 0 and intact_rescued_count > 0 and pseudo_count == 0 and absent_count == 0:
        cluster_state = "PARTIAL_INTACT"
    elif pseudo_count > 0 and absent_count == 0:
        cluster_state = "PARTIAL_PSEUDOGENIZED"
    else:
        cluster_state = "MIXED"

    cluster_row = {
        "target_id": args.target_name,
        "expected_gene_count": total_query_genes,
        "protein_detected_count": total_protein_detected,
        "location_unknown_count": location_unknown_count,
        "blastx_rescued_count": total_blastx_rescued,
        "intact_rescued_count": intact_rescued_count,
        "pseudogene_count": pseudo_count,
        "frameshift_count": total_frameshifts,
        "premature_stop_count": total_premature_stops,
        "mutation_evidence_status": ("VERSION_LIMITED" if any(g["mutation_evidence_status"] == "VERSION_LIMITED" for g in gene_results)
                                      else "ASSESSED" if any(g["mutation_evidence_status"] == "ASSESSED" for g in gene_results)
                                      else "NOT_ASSESSED"),
        "annotation_issue": annotation_issue,
        "pseudogene_detected": pseudogene_detected,
        "member_gene_states": ";".join(member_states),
        "cluster_state": cluster_state,
        "T_local": T_local,
        "T_global": T_global,
        "R_known": R_known,
        "I_novel": I_novel,
        "uncharacterized_intervening_count": n_novel,
        "novel_insertion_evidence": "NOT_ASSESSED",
        "replacement_search_status": replacement_search_status,
        "A_local": A_local,
        "A_global": A_global
    }
    cluster_row.update(classify_operon_structure(gene_results, gene_coordinates,
                                                list(query_dict), intervening_genes_list,
                                                args.max_operon_intergenic_bp))

    pd.DataFrame([cluster_row]).to_csv(args.out_cluster_tsv, sep="\t", index=False)
    print(f"Stage 2 Rescue & OFI Quantification Completed: cluster_state={cluster_state}, T_local={T_local}, T_global={T_global}, R_known={R_known}, I_novel={I_novel}, A_local={A_local}, A_global={A_global}", flush=True)

if __name__ == "__main__":
    main()
