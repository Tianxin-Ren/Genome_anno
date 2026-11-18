"""Exact sequence identity to directly assayed hydrogenase reference proteins."""
import hashlib, io, json, math
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
RULE_VERSION = "assay_sequence_identity_v1"
STATUS = "exact_assayed_core_sequence_supported_uncalibrated"
def _sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def _digest(x): return hashlib.sha256(json.dumps(x, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
def _seqsha(x): return hashlib.sha256(str(x).encode()).hexdigest()

def _model(root):
    root=Path(root); p=root/"research_workflow/assay_sequence_identity_model.json"; model=json.loads(p.read_text())
    specs=model.get("assay_records",[])
    if model.get("rule_version")!=RULE_VERSION or not specs or len({s.get("assay_accession") for s in specs})!=len(specs): raise ValueError("assay model mismatch")
    return model

def _genes(annotation):
    rows=[]
    with (Path(annotation)/"genes.tsv").open() as fh:
        header=fh.readline().rstrip("\n").split("\t")
        if "id" not in header: raise ValueError("genes header missing id")
        for line in fh:
            if line.strip():
                v=line.rstrip("\n").split("\t")
                if len(v)!=len(header): raise ValueError("genes column mismatch")
                rows.append(dict(zip(header,v)))
    ids=[r.get("id") for r in rows]
    if not ids or any(not x for x in ids) or len(ids)!=len(set(ids)): raise ValueError("genes IDs are not unique")
    return {r["id"]:r for r in rows}

def _actual_features(annotation,features):
    p=Path(annotation)/"confirmation/features.json"; actual=json.loads(p.read_text()); supplied=list(features.values()) if isinstance(features,dict) else list(features); byid={x.get("protein_id"):x for x in actual}
    sids=[x.get("protein_id") for x in supplied]
    if len(sids)!=len(set(sids)) or any(x not in byid for x in sids) or any(_digest(x)!=_digest(byid[x.get("protein_id")]) for x in supplied): raise ValueError("features do not match actual confirmation source")
    # The caller may intentionally provide only a core-feature subset.  Bind
    # every supplied row to the canonical confirmation row, but evaluate only
    # that subset (never silently broaden its claim to every feature).
    return actual,p,[byid[x] for x in sids]

def _domain_ok(feature,spec,protein_length):
    if feature.get("family")!=spec["family"]: return None
    for d in feature.get("domain_evidence",[]):
        cs=d.get("components",[])
        if d.get("accession")!=spec["domain_accession"] or d.get("supported") is not True or d.get("evidence_kind")!="single_domain" or len(cs)!=1: continue
        c=cs[0]
        nums=[d.get("model_coverage"),d.get("query_core_overlap"),c.get("sequence_bitscore"),c.get("sequence_ga"),c.get("domain_bitscore"),c.get("domain_ga"),c.get("query_start"),c.get("query_end"),c.get("query_length"),c.get("model_length")]
        if any(not isinstance(x,(int,float)) or not math.isfinite(x) for x in nums): continue
        if any(not 0 <= d.get(k, -1) <= 1 for k in ("model_coverage", "query_core_overlap")): continue
        if d.get("model_coverage",0)<spec["domain_min_coverage"] or d.get("query_core_overlap",0)<spec["domain_min_query_overlap"] or c.get("protein_id")!=feature.get("protein_id") or c.get("model_sha256")!=spec["domain_model_sha256"] or c.get("accession")!=spec["component_accession"] or c.get("profile")!=spec["component_profile"] or c.get("model_length")!=spec["component_model_length"] or c.get("query_length")!=protein_length or c.get("sequence_ga")!=spec["component_ga"] or c.get("domain_ga")!=spec["component_ga"]: continue
        if c.get("sequence_bitscore",0)<spec["component_ga"] or c.get("domain_bitscore",0)<spec["component_ga"] or c.get("query_start",0)<1 or c.get("query_end",0)<c.get("query_start",0) or c.get("query_end",0)>protein_length: continue
        return d
    return None

def _conflict(seq,spec,root):
    digest=_seqsha(seq); bank=json.loads((Path(root)/"references/processed/catalytic_confirmation/prototypes.json").read_text())
    return any(x.get("sequence_sha256")==digest and (x.get("label")=="negative" or (x.get("family") and x.get("family")!=spec["family"])) for x in bank)

def _evaluate(annotation,features,root):
    root=Path(root); model=_model(root); actual,feature_path,selected_features=_actual_features(annotation,features); genes=_genes(annotation)
    with (Path(annotation)/"proteins.faa").open() as fh: recs=list(SeqIO.parse(fh,"fasta"))
    if len({r.id for r in recs})!=len(recs): raise ValueError("protein IDs are not unique")
    proteins={r.id:str(r.seq).upper() for r in recs}
    if set(proteins)!=set(genes): raise ValueError("genes/proteins are not one-to-one")
    specs=model["assay_records"]; out={}
    for f in selected_features:
        pid=f["protein_id"]; seq=proteins[pid]
        for spec in specs:
            if _seqsha(seq)!=spec["assay_sequence_sha256"] or _conflict(seq,spec,root): continue
            domain=_domain_ok(f,spec,len(seq))
            if not domain: continue
            raw=genes[pid].get("partial"); out[pid]={"rule_version":RULE_VERSION,"status":STATUS,"protein_id":pid,"family":spec["family"],"protein_role":spec["protein_role"],"assay_accession":spec["assay_accession"],"assay_sequence_sha256":spec["assay_sequence_sha256"],"protein_sha256":_seqsha(seq),"protein_length":len(seq),"source_partial":str(raw).lower()=="true" if str(raw).lower() in ("true","false") else None,"raw_source_partial":raw,"source_annotation":str(Path(annotation).resolve()),"assay_source_record":spec["source_record"],"assay_source_sha256":spec["source_record_sha256"],"domain_accession":spec["domain_accession"],"domain_evidence":domain,"feature_sha256":_digest(f),"tag_status":"unknown_not_bound","mature_boundary_status":"unknown_not_bound","scope":spec["scope"]}; out[pid]["record_digest"]=_digest(out[pid]); break
    files=[Path(annotation)/"genes.tsv",Path(annotation)/"proteins.faa",feature_path,Path(__file__),root/"research_workflow/assay_sequence_identity_model.json",root/model["reference_bank"],root/model["threshold_policy"],root/model["thresholds"],root/"references/processed/catalytic_confirmation/prototypes.json"]
    return out,{str(p.resolve()):_sha(p) for p in files}

def prepare(annotation,features,output,root=ROOT):
    annotation,output,root=Path(annotation),Path(output),Path(root)
    if output.exists(): raise FileExistsError("output directory already exists")
    output.mkdir(parents=True)
    try:
        before_model=_model(root); files=[annotation/"genes.tsv",annotation/"proteins.faa",annotation/"confirmation/features.json",Path(__file__),root/"research_workflow/assay_sequence_identity_model.json",root/before_model["reference_bank"],root/before_model["threshold_policy"],root/before_model["thresholds"],root/"references/processed/catalytic_confirmation/prototypes.json"]
        before={str(p.resolve()):_sha(p) for p in files}; records,inputs=_evaluate(annotation,features,root); after={p:_sha(p) for p in before}
        if before!=after: raise RuntimeError("assay input changed during evaluation")
        dp=output/"assay_sequence_identity.json"; dp.write_text(json.dumps({"schema":"assay_sequence_identity_v1","records":list(records.values())},indent=2,sort_keys=True)+"\n"); manifest={"schema":"assay_sequence_identity_manifest_v1","status":"complete","rule_version":RULE_VERSION,"source_annotation":str(annotation.resolve()),"input_hashes":inputs,"output_sha256":_sha(dp),"model_sha256":_sha(root/"research_workflow/assay_sequence_identity_model.json")}; mp=output/"assay_sequence_identity_manifest.json"; mp.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n"); return records,{**inputs,str(dp.resolve()):_sha(dp),str(mp.resolve()):_sha(mp)}
    except Exception as e:
        (output/"assay_sequence_identity_manifest.json").write_text(json.dumps({"schema":"assay_sequence_identity_manifest_v1","status":"failed","error":str(e)},indent=2)+"\n"); raise

def verify_saved(output,annotation,features,root=ROOT):
    output,annotation,root=Path(output),Path(annotation),Path(root); m=json.loads((output/"assay_sequence_identity_manifest.json").read_text()); dp=output/"assay_sequence_identity.json"
    if m.get("status")!="complete" or m.get("source_annotation")!=str(annotation.resolve()) or m.get("output_sha256")!=_sha(dp): raise ValueError("assay evidence owner/status/output mismatch")
    current,inputs=_evaluate(annotation,features,root); data=json.loads(dp.read_text())
    if data.get("records")!=list(current.values()): raise ValueError("assay evidence does not recompute")
    if m.get("input_hashes")!=inputs: raise ValueError("assay input envelope changed")
    return current

def record_supported(record,feature,root=ROOT):
    if not isinstance(record,dict) or record.get("status")!=STATUS or record.get("rule_version")!=RULE_VERSION or record.get("record_digest")!=_digest({k:v for k,v in record.items() if k!="record_digest"}) or record.get("feature_sha256")!=_digest(feature): return False
    try:
        ann=Path(record["source_annotation"]); actual=json.loads((ann/"confirmation/features.json").read_text()); current,_=_evaluate(ann,actual,root); return current.get(record.get("protein_id"))==record
    except Exception: return False
