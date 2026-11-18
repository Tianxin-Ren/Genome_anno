"""Frozen prospective end-to-end evaluation of the hydrogenase workflow."""
from __future__ import annotations
import argparse, csv, hashlib, json, os, re, subprocess, sys, tempfile
from collections import Counter
from pathlib import Path
from Bio import SeqIO
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from evaluate_identity_corpus import identity_metrics
from validate_fes_prediction import wilson
from annotate import extract_genbank

PROTECTED = ROOT / "references/validation/round35_archaeal_sources/construct_correspondence.json"
EVALUATOR = Path(__file__).resolve()
SUPPORTED_TYPES = {"genome": "--genome", "genbank": "--genbank", "proteins": "--proteins"}
SUPPORTED_MODES = {"genome", "assembly", "bin", "transcript"}
IDENTITY_STATES = {"calibrated_high", "sequence_supported_uncalibrated_core", "sequence_supported_uncalibrated_noncanonical"}
HASH_RE = re.compile(r"^[0-9a-f]{64}$")
UNIT_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
FUNCTION_LABELS = {
    "candidate_H2_F420_reduction", "candidate_H2_THF_derivative_hydrogenation_potential",
    "candidate_H2_H4MPT_hydrogenation_potential", "candidate_H2_Hdr_coupling",
    "candidate_H2_pyridine_nucleotide_coupling", "candidate_H2_ferredoxin_compatible_core",
    "candidate_H2_quinone_coupling", "candidate_H2_formate_coupling",
    "candidate_genome_H2_formate_pathway", "candidate_H2_formate_pathway_conditional_on_recoding",
    "candidate_CO_linked_H2_pathway", "candidate_membrane_associated_H2_redox_system",
    "candidate_reference_reaction_potential", "reference_sequence_matched_reaction_potential",
}
CORE_ONLY_LABELS = {"candidate_H2_ferredoxin_compatible_core"}
MEMBRANE_ONLY_LABELS = {"candidate_membrane_associated_H2_redox_system"}
CONDITIONAL_LABELS = {"candidate_H2_formate_pathway_conditional_on_recoding"}
UNRESOLVED_LABELS = {"unresolved", "unassigned", "hydrogenase_family_supported_function_unresolved", "hydrogenase_candidate_identity_unresolved", "candidate_function_rule_missing", "H2_core_supported_partner_unresolved", "H2_core_supported_substrate_unresolved"}

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def canonical(sequence): return "".join(str(sequence).split()).upper().rstrip("*")
def sequence_hash(sequence): return hashlib.sha256(canonical(sequence).encode()).hexdigest()
def resolve_path(value):
    p = Path(value); return (p if p.is_absolute() else ROOT / p).resolve()
def _targets(value): return [] if value is None else [value] if isinstance(value, str) else list(value)
def is_specific_function(label):
    return str(label) in FUNCTION_LABELS - MEMBRANE_ONLY_LABELS - CONDITIONAL_LABELS
def function_label_diagnostic(label):
    if not label or label in UNRESOLVED_LABELS or "unresolved" in str(label).lower(): return "unresolved_function_label"
    if label in CORE_ONLY_LABELS: return "core_only_function_potential"
    if label in MEMBRANE_ONLY_LABELS: return "membrane_only_function_potential"
    if label in CONDITIONAL_LABELS: return "conditional_recoding_function_potential"
    if label in FUNCTION_LABELS: return "recognized_function_label"
    return "unrecognized_function_label"
def _metrics(truth, predicted): return identity_metrics(truth, set(predicted))
def metric(tp, fp, fn, tn):
    rows = ([{"unit_id": f"p{i}", "label": "positive"} for i in range(tp + fn)] +
            [{"unit_id": f"n{i}", "label": "negative"} for i in range(fp + tn)])
    predicted = {f"p{i}" for i in range(tp)} | {f"n{i}" for i in range(fp)}
    return identity_metrics(rows, predicted)
def _truth_rows(rows, positive): return [{"unit_id": r["unit_id"], "label": "positive" if positive(r) else "negative"} for r in rows]

def _protected_hashes(path):
    result = set()
    for row in json.loads(Path(path).read_text()):
        if isinstance(row, dict): result.update(row[k] for k in ("native_sequence_sha256", "construct_sequence_sha256") if row.get(k))
    return result
def _safe_id(value):
    return bool(value) and bool(UNIT_ID_RE.fullmatch(str(value))) and str(value) not in {".", ".."}

def _input_proteins(unit):
    inp, path = unit["input"], resolve_path(unit["input"]["path"])
    if inp["type"] == "proteins": return [str(r.seq) for r in SeqIO.parse(path, "fasta")]
    if inp["type"] == "genbank": return [g["sequence"] for g in extract_genbank(path)[0]]
    with tempfile.TemporaryDirectory(prefix="hydrogenase_benchmark_preview-") as td:
        td = Path(td); faa = td / "prodigal.faa"; prodigal = ROOT / ".env/bin/prodigal"
        if not prodigal.is_file(): raise ValueError("fixed Prodigal executable is unavailable for pre-query guard")
        with (td / "prodigal.log").open("w") as log:
            subprocess.run([str(prodigal), "-i", str(path), "-a", str(faa), "-o", str(td/"prodigal.gff"), "-f", "gff", "-p", "meta", "-g", "11"], check=True, stdout=subprocess.DEVNULL, stderr=log)
        return [str(r.seq) for r in SeqIO.parse(faa, "fasta")]
def _assert_input_guard(unit, protected):
    if any(sequence_hash(s) in protected for s in _input_proteins(unit)): raise ValueError("benchmark input contains a protected validation sequence: " + unit["unit_id"])

def validate_manifest(manifest_path):
    manifest_path = Path(manifest_path).resolve(); manifest = json.loads(manifest_path.read_text())
    if manifest.get("state") != "frozen_before_prediction": raise ValueError("benchmark manifest must be frozen_before_prediction")
    model = manifest.get("model_snapshot") or {}; model_path = resolve_path(model.get("path", ""))
    if not model_path.is_file() or not HASH_RE.fullmatch(str(model.get("sha256", ""))) or sha(model_path) != model["sha256"]: raise ValueError("model snapshot is missing or changed")
    methods = json.loads(model_path.read_text()).get("method_state", {})
    if not methods: raise ValueError("model snapshot has no method_state")
    for name, expected in methods.items():
        path = resolve_path(name)
        if not path.is_file() or sha(path) != expected: raise ValueError("frozen method changed: " + name)
    guard = manifest.get("protected_guard") or {}; protected_path = resolve_path(guard.get("metadata_path", str(PROTECTED)))
    if protected_path != PROTECTED.resolve() or not protected_path.is_file() or sha(protected_path) != guard.get("metadata_sha256"): raise ValueError("protected metadata guard missing or changed")
    if guard.get("sequence_content_read", True) is not False or guard.get("protected_hash_count") != len(_protected_hashes(protected_path)): raise ValueError("protected sequence guard is invalid")
    units = manifest.get("units")
    if not units: raise ValueError("manifest has no benchmark units")
    seen_units, seen_sources, protected = set(), set(), _protected_hashes(protected_path)
    for unit in units:
        uid, source = unit.get("unit_id"), unit.get("biological_source_id")
        if not _safe_id(uid) or uid in seen_units or not source or source in seen_sources: raise ValueError("unit_id and biological_source_id must be unique and path-safe")
        seen_units.add(uid); seen_sources.add(source)
        admission = unit.get("admission") or {}; ap = resolve_path(admission.get("path", ""))
        if not ap.is_file() or sha(ap) != admission.get("sha256"): raise ValueError("admission artifact missing or changed: " + str(uid))
        inp = unit.get("input") or {}; typ, mode = inp.get("type"), inp.get("mode")
        if typ not in SUPPORTED_TYPES or mode not in SUPPORTED_MODES: raise ValueError("unsupported input type or mode: " + str(uid))
        if mode == "bin" and not inp.get("bin_map"): raise ValueError("bin mode requires a bound bin_map: " + str(uid))
        ip = resolve_path(inp.get("path", ""))
        if not ip.is_file() or sha(ip) != inp.get("sha256"): raise ValueError("input artifact missing or changed: " + str(uid))
        if inp.get("bin_map"):
            bm, bp = inp["bin_map"], resolve_path(inp["bin_map"].get("path", ""))
            if mode != "bin" or not bp.is_file() or sha(bp) != bm.get("sha256"): raise ValueError("bin_map is missing or not bound to bin mode: " + str(uid))
        if typ == "proteins" and inp.get("gff"):
            gp = resolve_path(inp["gff"].get("path", ""))
            if not gp.is_file() or sha(gp) != inp["gff"].get("sha256"): raise ValueError("protein GFF is missing or changed: " + str(uid))
        truth = unit.get("truth") or {}
        if not isinstance(truth.get("hydrogenase"), bool): raise ValueError("truth.hydrogenase must be boolean: " + str(uid))
        if not isinstance(truth.get("function_evaluable"), bool): raise ValueError("truth.function_evaluable must be boolean: " + str(uid))
        targets = _targets(truth.get("target_AA_sha256"))
        if not targets or any(not isinstance(h, str) or not HASH_RE.fullmatch(h) for h in targets): raise ValueError("truth.target_AA_sha256 must contain SHA256 protein hashes: " + str(uid))
        if truth.get("function_evaluable") and not _targets(truth.get("compatible_function_predictions")): raise ValueError("function_evaluable truth requires compatible function labels: " + str(uid))
        _assert_input_guard(unit, protected)
    return manifest, model_path, methods

def _catalog(path):
    with Path(path).open() as handle: return {row["candidate_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

def native_locus_targets(unit, genes_path):
    """Map a native locus without assuming an identical predicted start codon."""
    locus = unit['truth'].get('target_locus')
    if not locus: return []
    inp = unit['input']
    if inp['type'] not in ('genome','genbank'): raise ValueError('native-locus mapping requires DNA input')
    records = SeqIO.parse(resolve_path(inp['path']), 'fasta' if inp['type']=='genome' else 'genbank')
    contigs = {r.id: len(r.seq) for r in records if sequence_hash(r.seq)==locus['contig_sequence_sha256']}
    start,end,strand = (locus[k] for k in ('start','end','strand'))
    if (not contigs or not all(type(x) is int for x in (start,end,strand)) or
            strand not in (-1,1) or not 1<=start<=end<=min(contigs.values())):
        raise ValueError('native target locus does not belong to the bound DNA input')
    with Path(genes_path).open() as handle: genes=list(csv.DictReader(handle,delimiter='\t'))
    targets=[]
    for gene in genes:
        a,b,s=(int(gene[k]) for k in ('start','end','strand'))
        contained=(start<=a<=b<=end or a<=start<=end<=b)
        same_frame=((a-start) if strand==1 else (b-end))%3==0
        if gene['contig'] in contigs and s==strand and contained and same_frame: targets.append(gene['id'])
    return targets

def score_unit(unit, output_dir):
    output_dir = Path(output_dir); proteins, catalog_path = output_dir / "annotation/proteins.faa", output_dir / "candidate_function_catalog.tsv"
    if not proteins.is_file() or not catalog_path.is_file(): raise FileNotFoundError("workflow output is incomplete: " + str(output_dir))
    sequences = {}
    with proteins.open() as handle:
        for rec in SeqIO.parse(handle, "fasta"): sequences.setdefault(sequence_hash(rec.seq), []).append(rec.id)
    truth = unit["truth"]; target_hashes = _targets(truth["target_AA_sha256"]); exact_ids = [pid for h in target_hashes for pid in sequences.get(h, [])]
    target_ids=native_locus_targets(unit,output_dir/'annotation/genes.tsv') if truth.get('target_locus') else exact_ids
    exact_ids=[pid for pid in exact_ids if pid in target_ids]
    mapping='exact_protein_sequence' if exact_ids else 'native_locus_same_frame_containment' if target_ids else 'unmatched'
    catalog = _catalog(catalog_path); selected = [catalog[pid] for pid in target_ids if pid in catalog]
    retained=bool(selected);identity=any(r.get('identity_status') in IDENTITY_STATES for r in selected)
    labels=sorted({r.get('function_prediction','') for r in selected});specific={label for label in labels if is_specific_function(label)}
    label=labels[0] if len(labels)==1 else json.dumps(labels) if labels else ''
    assigned=bool(specific);compatible=set(_targets(truth.get('compatible_function_predictions')))
    function_evaluable=bool(truth.get('function_evaluable')) and truth.get('hydrogenase') is True and bool(compatible)
    # Do not select a favourable CDS when duplicate or fragmented calls make several claims.
    function_correct=bool(assigned and specific<=compatible)
    diagnostic = ("target_AA_and_native_locus_not_found" if truth.get('target_locus') else 'target_AA_not_found_in_called_proteome') if not target_ids else "target_not_retained_by_catalog" if not retained else 'native_gene_boundary_differs' if not exact_ids else ""
    return {"unit_id": unit["unit_id"], "biological_source_id": unit["biological_source_id"], "independence_stratum": unit.get("independence_stratum", "unspecified"), "input_mode": unit["input"]["mode"], "truth_hydrogenase": truth["hydrogenase"], "function_evaluable": function_evaluable, "target_protein_id": ';'.join(target_ids), "target_hashes_found": len(exact_ids), 'target_locus_matches':len(target_ids) if not exact_ids else 0, 'target_mapping_method':mapping, 'multiple_target_calls':len(target_ids)>1, "catalog_candidate_count": len(catalog), "unscored_catalog_candidates": len(set(catalog)-set(target_ids)), "retained": retained, "identity_supported": identity, "calibrated_high": any(r.get('identity_status')=='calibrated_high' for r in selected), "function_assigned": assigned, "function_correct": function_correct, "predicted_function": label, 'all_predicted_functions':labels, "function_label_diagnostic": function_label_diagnostic(label) if len(labels)<=1 else 'multiple_calls_all_claims_scored', "compatible_functions": sorted(compatible), "diagnostic": diagnostic}

def score_metrics(rows):
    discovery = _metrics(_truth_rows(rows, lambda r: r["truth_hydrogenase"]), {r["unit_id"] for r in rows if r["retained"]}); truth_identity = _truth_rows(rows, lambda r: r["truth_hydrogenase"]); identity = _metrics(truth_identity, {r["unit_id"] for r in rows if r["identity_supported"]}); calibrated = _metrics(truth_identity, {r["unit_id"] for r in rows if r["calibrated_high"]})
    function_rows = [r for r in rows if r["function_evaluable"]]; function_truth = _truth_rows(function_rows + [r for r in rows if not r["truth_hydrogenase"]], lambda r: r["function_evaluable"]); assigned = _metrics(function_truth, {r["unit_id"] for r in rows if r["function_assigned"]}); wrong = sum(r["function_assigned"] and not r["function_correct"] and r["function_evaluable"] for r in rows); correct = sum(r["function_correct"] for r in function_rows); fn = sum(not r["function_correct"] for r in function_rows); function_fp = assigned["FP"] + wrong; function = dict(assigned); function.update(TP=correct, FP=function_fp, FN=fn, precision=correct/(correct+function_fp) if correct+function_fp else None, recall=correct/(correct+fn) if correct+fn else None, precision_wilson95=wilson(correct, correct+function_fp), recall_wilson95=wilson(correct, correct+fn), wrong_function_FP_and_FN=wrong)
    denom = len(function_rows); coverage = sum(r["function_assigned"] for r in function_rows) / denom if denom else None; recall = correct / denom if denom else None; ft = [{"unit_id": r["unit_id"], "label": "positive"} for r in function_rows]
    negative_rows = [{"unit_id": r["unit_id"], "label": "negative"} for r in rows if not r["truth_hydrogenase"]]
    negative = _metrics(negative_rows, {r["unit_id"] for r in rows if r["function_assigned"] and not r["truth_hydrogenase"]})
    positives = [r for r in rows if r['truth_hydrogenase']]
    all_assigned = sum(r['function_assigned'] for r in positives)
    return {"discovery": discovery, "identity": identity, "identity_calibrated_high": calibrated, "function": function, "function_evaluable_positive": function, "function_assignment": assigned, "negative_function": negative,
            "positive_truth_count": len(positives), "function_evaluable_positive_count": denom,
            "automatic_assignment_coverage_all_positive": all_assigned/len(positives) if positives else None,
            "automatic_assignment_coverage_all_positive_wilson95": wilson(all_assigned,len(positives)),
            "function_assignment_coverage": coverage, "function_assignment_coverage_wilson95": identity_metrics(ft, {r["unit_id"] for r in function_rows if r["function_assigned"]})["recall_wilson95"] if denom else None, "function_end_to_end_recall": recall, "function_end_to_end_recall_wilson95": identity_metrics(ft, {r["unit_id"] for r in function_rows if r["function_correct"]})["recall_wilson95"] if denom else None, "partial_truth_unscored_catalog_candidates": sum(r["unscored_catalog_candidates"] for r in rows)}

def _run_unit(unit, out, methods, frozen_context):
    inp, out_unit = unit["input"], out / unit["unit_id"]
    if out_unit.exists(): raise FileExistsError("benchmark output already exists: " + str(out_unit))
    args = [str(ROOT / "run_hydrogenase.py"), SUPPORTED_TYPES[inp["type"]], str(resolve_path(inp["path"])), "--mode", inp["mode"], "--out", str(out_unit), "--threads", "2"]
    if inp["type"] == "proteins" and inp.get("gff"): args.extend(["--gff", str(resolve_path(inp["gff"]["path"]))])
    if inp.get("bin_map"): args.extend(["--bin-map", str(resolve_path(inp["bin_map"]["path"]))])
    before = {n: sha(resolve_path(n)) for n in methods}
    if before != frozen_context['methods']: raise RuntimeError('frozen method changed before run')
    bindings = {"manifest": sha(frozen_context["manifest_path"]), "evaluator": sha(EVALUATOR), "admission": sha(resolve_path(unit["admission"]["path"])), "input": sha(resolve_path(inp["path"]))}
    if inp.get("gff"): bindings["gff"] = sha(resolve_path(inp["gff"]["path"]))
    if inp.get("bin_map"): bindings["bin_map"] = sha(resolve_path(inp["bin_map"]["path"]))
    if bindings != frozen_context["bindings"][unit["unit_id"]]: raise RuntimeError("benchmark binding changed before run")
    env = dict(os.environ); env["PATH"] = str(ROOT / ".env/bin") + os.pathsep + env.get("PATH", "")
    proc = subprocess.run([str(ROOT / ".env/bin/python"), *args], cwd=ROOT, env=env, capture_output=True, text=True); after = {n: sha(resolve_path(n)) for n in methods}
    after_bindings = {"manifest": sha(frozen_context["manifest_path"]), "evaluator": sha(EVALUATOR), "admission": sha(resolve_path(unit["admission"]["path"])), "input": sha(resolve_path(inp["path"]))}
    if inp.get("gff"): after_bindings["gff"] = sha(resolve_path(inp["gff"]["path"]))
    if inp.get("bin_map"): after_bindings["bin_map"] = sha(resolve_path(inp["bin_map"]["path"]))
    if before != after or before != frozen_context["methods"] or after_bindings != bindings: raise RuntimeError("frozen method or benchmark binding changed during run")
    (out / (unit["unit_id"] + ".benchmark_command.json")).write_text(json.dumps({"command": [str(ROOT / ".env/bin/python"), *args], "returncode": proc.returncode, "stdout": proc.stdout, "stderr": proc.stderr}, ensure_ascii=False, indent=2) + "\n")
    if proc.returncode: raise RuntimeError("workflow failed for " + unit["unit_id"])
    return out_unit

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--manifest", type=Path, required=True); parser.add_argument("--out", type=Path, required=True); args = parser.parse_args(argv); manifest_path, out = args.manifest.resolve(), args.out.resolve()
    if out.exists(): raise SystemExit("Use a fresh benchmark output directory")
    manifest, model_path, methods = validate_manifest(manifest_path)
    frozen_context = {"methods": {n: sha(resolve_path(n)) for n in methods}, "manifest": sha(manifest_path), "manifest_path": manifest_path, "evaluator": sha(EVALUATOR), "bindings": {}}
    for unit in manifest["units"]:
        inp = unit["input"]; b = {"manifest": frozen_context["manifest"], "evaluator": frozen_context["evaluator"], "admission": sha(resolve_path(unit["admission"]["path"])), "input": sha(resolve_path(inp["path"]))}
        if inp.get("gff"): b["gff"] = sha(resolve_path(inp["gff"]["path"]))
        if inp.get("bin_map"): b["bin_map"] = sha(resolve_path(inp["bin_map"]["path"]))
        frozen_context["bindings"][unit["unit_id"]] = b
    out.mkdir(parents=True); rows = []
    for unit in manifest["units"]: rows.append(score_unit(unit, _run_unit(unit, out, methods, frozen_context)))
    grouped = {"all": rows}
    for key in ("input_mode", "independence_stratum"):
        for value in sorted({r[key] for r in rows}): grouped[f"{key}:{value}"] = [r for r in rows if r[key] == value]
    metrics = {name: score_metrics(selected) for name, selected in grouped.items()}
    with (out / "per_unit.tsv").open("w", newline="") as handle:
        fields = list(rows[0]); writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in row.items()} for row in rows)
    report = {"status": "complete", "manifest": str(manifest_path), "manifest_sha256": frozen_context["manifest"], "evaluator_sha256": frozen_context["evaluator"], "model_snapshot": str(model_path), "model_snapshot_sha256": sha(model_path), "units": len(rows), "source_count": len({r["biological_source_id"] for r in rows}), "independence_strata": dict(Counter(r["independence_stratum"] for r in rows)), "metrics": metrics, "artifacts": {"per_unit.tsv": sha(out / "per_unit.tsv")}, "limitations": ["Only reviewed target proteins are scored; other genomic candidates are unscored partial truth, not genome-wide negatives.", "Function metrics require nonempty compatible labels and remain computational predictions."]}
    (out / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n"); print(json.dumps({"status": report["status"], "units": len(rows), "metrics": metrics}, ensure_ascii=False))
if __name__ == "__main__": main()
