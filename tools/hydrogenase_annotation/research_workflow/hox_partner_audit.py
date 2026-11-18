#!/usr/bin/env python3
"""Independent HoxF/HoxU partner audit (raw search and linkage provenance).

This deliberately consumes the existing full-protein search output and module
graph.  It never turns remote co-occurrence into a physical-complex claim.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, ast
from pathlib import Path

VERSION = "hox-partner-audit-v1"
ROLES = ("HoxF", "HoxU")

def _blast_loader():
    """Use the production parser when dependencies are available; otherwise
    load its literal BLAST_FIELDS assignment and apply the same named-field
    semantics (never positional guesses)."""
    try:
        from scripts.annotate import BLAST_FIELDS, load_blast
        return BLAST_FIELDS, load_blast
    except ModuleNotFoundError:
        tree=ast.parse(Path(__file__).parents[1].joinpath("scripts/annotate.py").read_text())
        node=next(n.value for n in ast.walk(tree) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=="BLAST_FIELDS" for t in n.targets))
        fields=ast.literal_eval(node.func.value).split()
        def load(path, thresholds, references=None):
            out=[]
            refs={r["reference_id"]:r for r in references or []}
            for line in Path(path).read_text().splitlines():
                h=dict(zip(fields,line.split("\t")))
                for k in ("pident","evalue","bitscore"): h[k]=float(h[k])
                for k in ("length","qstart","qend","sstart","send","qlen","slen"): h[k]=int(h[k])
                h["qcov"]=(h["qend"]-h["qstart"]+1)/h["qlen"]
                h["scov"]=(h["send"]-h["sstart"]+1)/h["slen"]
                h["coverage_scope"]=refs.get(h["sseqid"],{}).get("coverage_scope","full_reference_protein")
                h["rule_query_coverage"]=h["qcov"]
                h["passes_thresholds"]=(h["evalue"]<=thresholds["blast_evalue"] and h["pident"]>=thresholds["min_identity"] and h["rule_query_coverage"]>=thresholds["min_query_coverage"] and h["scov"]>=thresholds["min_reference_coverage"])
                out.append(h)
            return out
        return fields, load

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def rows(path: Path):
    if not path.exists(): return []
    with path.open() as f: return list(csv.DictReader(f, delimiter="\t"))

def classify(pass_hits, fragment_hits, competitors, linked, complete=True):
    """Small, testable status policy; incomplete searches are never absence."""
    if not complete: return "search_incomplete"
    if pass_hits and not linked: return "unlinked"
    if len(pass_hits) > 1: return "pairing_ambiguous"
    if pass_hits and pass_hits[0].get("query"):
        competitors=[x for x in competitors if not x.get("query") or x.get("query")==pass_hits[0].get("query")]
    if pass_hits and competitors:
        best, rival = pass_hits[0].get("bitscore", 0), competitors[0].get("bitscore", 0)
        if best and rival >= .9 * best: return "competing_family"
    if pass_hits: return "supported"
    if fragment_hits: return "fragment_only"
    if competitors: return "competing_family"
    return "not_detected"

def _float(x):
    try: return float(x)
    except (TypeError, ValueError): return 0.0

def _json(value, default):
    try: return json.loads(value) if value else default
    except (TypeError, ValueError): return default

def audit(root: Path, out: Path):
    fields, load_blast = _blast_loader()
    thresholds=json.loads((root/"config/thresholds.json").read_text())
    ref_json = root / "references/processed/references.json"
    refs = json.loads(ref_json.read_text())
    hox_refs = {x["role"]: x for x in refs if x.get("role") in ROLES}
    ref_by_id = {x.get("reference_id"): x for x in refs}
    ref_faa = root / "references/processed/references.faa"
    samples = root / "results/hadar_pilot_round49_Fes_expansion/current"
    records, inputs = [], []
    for sample in sorted(p for p in samples.iterdir() if p.is_dir()):
        anchors = [x for x in rows(sample/"catalytic_anchors.tsv") if x.get("structure_id")=="NiFe_HoxFUYH"]
        if not anchors: continue
        genes = {x.get("protein_id") or x.get("id"): x for x in rows(sample/"genes.tsv")}
        nodes = rows(sample/"module_graph_nodes.tsv")
        manifest=_json((sample/"run_manifest.json").read_text(),{}) if (sample/"run_manifest.json").exists() else {}
        search_complete=manifest.get("status")=="complete" and (sample/"blast.tsv").exists() and bool((sample/"blast.tsv").stat().st_size)
        # Existing full-protein BLAST is retained verbatim as raw evidence.
        raw_hits = load_blast(sample/"blast.tsv", thresholds, refs) if (sample/"blast.tsv").exists() else []
        inputs += [sample/"proteins.faa", sample/"genes.tsv", sample/"blast.tsv", sample/"module_graph_nodes.tsv"]
        for anchor in anchors:
            aid = anchor["protein_id"]
            ag = genes.get(aid, {}) or anchor
            for role in ROLES:
                rr = hox_refs[role]
                # module graph rows preserve attempted role resources (including
                # split/empty assignments) and are not confused with hits.
                nrows = [n for n in nodes if n.get("anchor")==aid and n.get("role")==role]
                candidates = []
                for n in nrows:
                    ids=_json(n.get("protein_ids"), [])
                    if isinstance(ids,str): ids=[ids]
                    for q in ids:
                        candidates.append({"query":q,"coverage":_float(n.get("reference_coverage")),"identity":_float(n.get("aligned_pair_identity")),
                            "reference_start":n.get("reference_start"),"reference_end":n.get("reference_end"),
                            "source":"module_graph_nodes","selected":n.get("realization"),"candidate_realizations":_json(n.get("candidate_realizations"),[])})
                # Parse all raw BLAST lines and retain every Hox role hit.  The
                # current candidate gate is identity>=30, query/ref coverage>=.65.
                for h0 in raw_hits:
                    if h0["sseqid"] != rr.get("reference_id"): continue
                    h={"query":h0["qseqid"],"reference_id":h0["sseqid"],"identity":h0["pident"],"alignment_length":h0["length"],
                       "q_start":h0["qstart"],"q_end":h0["qend"],"s_start":h0["sstart"],"s_end":h0["send"],
                       "evalue":h0["evalue"],"bitscore":h0["bitscore"],"query_length":h0["qlen"],"reference_length":h0["slen"],
                       "query_coverage":h0["qcov"],"reference_coverage":h0["scov"],"passes_thresholds":h0["passes_thresholds"],"source":"blast.tsv"}
                    candidates.append(h)
                candidates.sort(key=lambda x: (_float(x.get("bitscore")),_float(x.get("coverage"))), reverse=True)
                for h in candidates:
                    g = genes.get(h.get("query"), {})
                    if g:
                        h["query_coordinate"] = {k:g.get(k) for k in ("contig","start","end","strand","bin_id")}
                pass_hits=[x for x in candidates if x.get("passes_thresholds")]
                fragment=[x for x in candidates if _float(x.get("identity"))>=thresholds["min_identity"] and x not in pass_hits]
                # Competitor records are explicit alternative roles/families only;
                # no competitor is invented from remote presence.
                competitors=[]
                # Preserve same-query alternative families from the complete
                # search, rather than inferring competition from co-location.
                target_queries={x.get("query") for x in pass_hits}
                for h0 in raw_hits:
                    if h0["qseqid"] not in target_queries or not h0["passes_thresholds"]: continue
                    alt=ref_by_id.get(h0["sseqid"])
                    if not alt or alt.get("role") in (role,"HoxF","HoxU"): continue
                    competitors.append({"query":h0["qseqid"],"reference_id":h0["sseqid"],"role":alt.get("role"),"identity":h0["pident"],"bitscore":h0["bitscore"],"query_coverage":h0["qcov"],"reference_coverage":h0["scov"],"source":"blast.tsv"})
                competitors.sort(key=lambda x:_float(x.get("bitscore")), reverse=True)
                linked = bool(pass_hits)
                local = False
                for h in pass_hits:
                    g=genes.get(h.get("query"),{})
                    if not g or g.get("contig") != ag.get("contig"): linked=False
                    elif g.get("strand") == ag.get("strand") and abs(int(g.get("start",0))-int(ag.get("end",0))) <= thresholds["local_max_span_bp"]:
                        local=True
                status=classify(pass_hits,fragment,competitors,linked,search_complete)
                records.append({"sample":sample.name,"anchor":aid,"role":role,"status":status,"search_status":"complete" if search_complete else "incomplete",
                    "anchor_coordinate":{"contig":ag.get("contig"),"start":ag.get("start"),"end":ag.get("end"),"strand":ag.get("strand"),"bin_id":ag.get("bin_id")},
                    "reference":{"reference_id":rr.get("reference_id"),"sequence_sha256":rr.get("sequence_sha256"),"length":rr.get("source_protein_length")},
                    "all_candidate_hits":candidates,"passing_hits":pass_hits,"fragment_hits":fragment,"competing_family_hits":competitors,
                    "resource_occupancy":{"module_rows":nrows,"query_proteins":sorted({x.get("query") for x in candidates if x.get("query")}),"fusion_or_split":len(nrows)>1},
                    "linkage":{"same_contig":linked,"local_same_strand_span":local,"local_max_span_bp":thresholds["local_max_span_bp"],"physical_complex_not_inferred":True,"association_mode":"genome","mode_ownership":"genome association only; no local-cluster ownership or physical-complex ownership"},
                    "rule_version":VERSION})
    out.mkdir(parents=True,exist_ok=True)
    input_hashes={str(p):sha256(p) for p in inputs if p.exists()}
    payload={"audit_version":VERSION,"records":records,"inputs":input_hashes,"references_sha256":sha256(ref_json),"reference_fasta_sha256":sha256(ref_faa) if ref_faa.exists() else None,
             "scope":"4 Hox anchors x HoxF/HoxU; existing full-protein BLAST plus module graph raw evidence; no remote co-occurrence claim"}
    (out/"audit.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    with (out/"audit.tsv").open("w") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow(["sample","anchor","role","status","passing_hits","fragment_hits","competitors","search_status"])
        for r in records: w.writerow([r["sample"],r["anchor"],r["role"],r["status"],len(r["passing_hits"]),len(r["fragment_hits"]),len(r["competing_family_hits"]),r["search_status"]])
    return payload

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--root",type=Path,default=Path(__file__).parents[1]); ap.add_argument("--out",type=Path,default=Path("experiments/round60_hox_audit")); a=ap.parse_args(); p=a.out if a.out.is_absolute() else a.root/a.out; audit(a.root,p)
if __name__=="__main__": main()
