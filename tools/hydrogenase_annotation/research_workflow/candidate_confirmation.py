"""Confirm high-recall candidates with the existing identity evidence stack.

This is a research-only, uncalibrated extension.  It never edits a core
annotation and never emits ``calibrated_high`` or a reaction/function claim.
"""
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from candidate_screen import sha as file_sha, table as read_table

FAMILY_MAP = {
    "NiFe": "NiFe", "FeFe": "FeFe", "Fe": "Fe",
    "NiFe_related": "NiFe", "NiFe_or_Nuo_related": "NiFe",
    "FeFe_related": "FeFe", "Hmd_related": "Fe",
}
SUPPORTED = {"NiFe", "FeFe", "Fe"}


def _write_table(path, rows, fields=None):
    rows = list(rows)
    fields = fields or sorted({k for row in rows for k in row})
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({k: row.get(k, "") for k in fields} for row in rows)


def _json_field(value):
    if isinstance(value, str):
        try:
            return json.loads(value)
        except (TypeError, json.JSONDecodeError):
            return value
    return value


def _candidate_family(row):
    values = _json_field(row.get("family_candidates", []))
    if isinstance(values, str):
        values = [values]
    normalized = {FAMILY_MAP.get(value, "") for value in values}
    normalized.discard("")
    return sorted(normalized), sorted(values)


def _copy_subset(annotation, work, ids):
    """Copy source metadata, then replace query-owned files with the subset."""
    work.mkdir(parents=True)
    for name in ["family_candidates.tsv", "rejected.tsv"]:
        source = annotation / name
        if source.exists():
            shutil.copy2(source, work / name)
    genes = [r for r in read_table(annotation / "genes.tsv") if r.get("id") in ids]
    _write_table(work / "genes.tsv", genes, list(read_table(annotation / "genes.tsv")[0]))
    source_faa = annotation / "proteins.faa"
    from Bio import SeqIO
    with source_faa.open() as handle:
        records = [r for r in SeqIO.parse(handle, "fasta") if r.id in ids]
    SeqIO.write(records, work / "proteins.faa", "fasta")
    (work / "candidate_subset_manifest.json").write_text(json.dumps({
        "status": "candidate_subset_prepared",
        "source_annotation": str(annotation.resolve()),
        "candidate_proteins": sorted(ids),
        "genes_sha256": file_sha(work / "genes.tsv"),
        "proteins_sha256": file_sha(work / "proteins.faa"),
        "scope": "Subset evidence input only; not a complete core annotation run."
    }, indent=2) + "\n")
    return genes


def _verify_screen(screen_dir):
    manifest_path = screen_dir / "candidate_screen_manifest.json"
    if not manifest_path.exists():
        raise ValueError("Missing candidate_screen_manifest.json")
    manifest = json.loads(manifest_path.read_text())
    for name, digest in manifest.get("output_sha256", {}).items():
        if file_sha(screen_dir / name) != digest:
            raise ValueError("Candidate screen output changed: " + name)
    annotation = Path(manifest["annotation"])
    if not annotation.exists():
        raise ValueError("Candidate screen annotation source is unavailable")
    if manifest.get("status") != "complete":
        raise ValueError("Candidate screen is not complete")
    for path, digest in manifest.get("input_sha256", {}).items():
        if not Path(path).is_file() or file_sha(Path(path)) != digest:
            raise ValueError("Candidate screen input changed: " + path)
    if not manifest.get("input_sha256"):
        raise ValueError("Candidate screen lacks bound input hashes")
    return manifest, annotation


def confirm(annotation_dir, screen_dir, out_dir, threads=2):
    screen_dir = Path(screen_dir).resolve(); out_dir = Path(out_dir).resolve()
    manifest, source_annotation = _verify_screen(screen_dir)
    source_inputs = {p: file_sha(Path(p)) for p in manifest["input_sha256"]}
    annotation_dir = Path(annotation_dir).resolve()
    if annotation_dir != source_annotation.resolve():
        raise ValueError("annotation_dir does not match the candidate screen source")
    if out_dir.exists():
        raise FileExistsError("Use a fresh candidate-confirmation output directory")
    screened = read_table(screen_dir / "candidate_screen.tsv")
    core = {r["protein_id"] for r in read_table(annotation_dir / "catalytic_anchors.tsv")}
    extras = [r for r in screened if r["protein_id"] not in core]
    out_dir.mkdir(parents=True)
    decisions = []
    evidence = {}
    backbone_inputs = {}
    fe_inputs = {}
    profile_inputs = {}
    experimental_backbone_inputs = {}
    assay_inputs = {}
    coding_rows = {}
    coding_inputs = {}
    nife_region_inputs = {}
    nife_region_proofs = {}
    for row in extras:
        pid = row["protein_id"]
        families, raw_families = _candidate_family(row)
        base = dict(candidate_id=pid, raw_family_candidates=raw_families,
                    normalized_families=families, source_sets=_json_field(row.get("source_sets", [])),
                    screen_status=row.get("discovery_status", ""), rejected_reason=row.get("rejected_reason", ""),
                    identity_status="unconfirmed", function_status="unassigned",
                    sequence_status="unconfirmed", blocking_reasons=[], identity_reasons=[])
        if len(families) != 1:
            base["blocking_reasons"] = ["multi_family_ambiguity" if len(families) > 1 else "family_unresolved"]
            base["identity_reasons"] = list(base["blocking_reasons"])
            decisions.append(base); continue
        family = families[0]; base["family"] = family
        if family not in SUPPORTED:
            base["blocking_reasons"] = ["family_outside_identity_calibration_scope"]
            base["identity_reasons"] = list(base["blocking_reasons"])
            decisions.append(base); continue
        if row.get("rejected_reason") and "non_catalytic" in row["rejected_reason"]:
            base["blocking_reasons"] = ["core_rejected_non_catalytic_homolog"]
            base["identity_reasons"] = list(base["blocking_reasons"])
            decisions.append(base); continue
        decisions.append(base)

    runnable = [d for d in decisions if d.get("family") in SUPPORTED and not d["blocking_reasons"]]
    if runnable:
        ids = {d["candidate_id"] for d in runnable}
        from coding_completeness import prepare as prepare_coding, load_evidence as load_coding
        coding_dir = out_dir / 'coding_completeness'
        _, coding_inputs = prepare_coding(annotation_dir, coding_dir, candidate_ids=ids, root=ROOT)
        _, coding_rows = load_coding(coding_dir, annotation=annotation_dir, root=ROOT)
        coding_inputs.update({str(p): file_sha(p) for p in [coding_dir/'coding_completeness.json', coding_dir/'coding_completeness_manifest.json']})
        work = out_dir / "candidate_annotation"
        subset_genes=_copy_subset(annotation_dir, work, ids)
        from active_centers import annotate_active_centers
        (work / "active_centers").mkdir()
        annotate_active_centers(work / "proteins.faa", work / "active_centers", threads=threads,
                                root=ROOT, anchors=ids)
        # Keep the exact filenames expected by extract_features.
        for name in ["active_center_blast.tsv", "active_center_summary.tsv", "active_center_template_evidence.tsv"]:
            source = work / "active_centers" / name
            if source.exists(): shutil.copy2(source, work / name)
        from noncanonical_center_evidence import annotate_evidence
        fe_ids = {d["candidate_id"] for d in runnable if d.get("family") == "FeFe"}
        if fe_ids:
            partial_ids={g['id'] for g in subset_genes if g.get('partial') in [True,'True']}
            annotate_evidence(work / "proteins.faa", work / "noncanonical_centers", fe_ids,partial_ids,
                              threads=threads, root=ROOT)
        else:
            (work / "noncanonical_centers").mkdir()
        anchors = [dict(protein_id=d["candidate_id"], catalytic_family=d["family"],
                        catalytic_identity_confidence="not_high_confidence") for d in runnable]
        from catalytic_confirmation import extract_features
        features = extract_features(work, work / "confirmation", ROOT, threads=threads, anchors=anchors)
        by_feature = {f["protein_id"]: f for f in features}
        from identity_scope import reference_backbones
        backbones,backbone_inputs = reference_backbones(work, by_feature, ROOT)
        from canonical_hmd import identity_evidence
        fe_evidence,fe_inputs=identity_evidence(work,by_feature,subset_genes,ROOT,coding_completeness=(coding_rows,coding_inputs))
        from explain_candidates import sequence_status, load_lineage_review, identity_reasons
        review = load_lineage_review(ROOT)
        rules = json.loads((ROOT / "config/catalytic_confirmation_thresholds.json").read_text())["families"]
        from identity_scope import NIFE_CORE_REGION_STATUS, nife_core_regions
        nife_region_proofs, nife_region_inputs = nife_core_regions(work, by_feature, ROOT)
        from fefe_profile_ligands import prepare as prepare_profile_ligands
        profile_evidence,profile_inputs=prepare_profile_ligands(work,by_feature,out_dir/'fefe_profile_ligands',rules.get('FeFe',{}),ROOT)
        from fefe_experimental_backbone import prepare as prepare_experimental_backbone
        experimental_backbone_evidence,experimental_backbone_inputs=prepare_experimental_backbone(work,by_feature,out_dir/'fefe_experimental_backbone',rules.get('FeFe',{}),ROOT)
        assay_evidence = {}
        if fe_ids:
            from assay_sequence_identity import prepare as prepare_assay
            assay_evidence,assay_inputs=prepare_assay(work,by_feature,out_dir/'assay_sequence_identity',ROOT)
        alternatives = {r["protein_id"]: r for r in read_table(work / "noncanonical_centers/alternatives.tsv")} if (work / "noncanonical_centers/alternatives.tsv").exists() else {}
        from Bio import SeqIO
        with (work/'proteins.faa').open() as handle:
            protein_sha256={r.id:hashlib.sha256(str(r.seq).upper().encode()).hexdigest() for r in SeqIO.parse(handle,'fasta')}
        genes_by_id = {r["id"]: dict(r,protein_sha256=protein_sha256[r['id']]) for r in read_table(work / "genes.tsv")}
        for decision in runnable:
            feature = by_feature.get(decision["candidate_id"])
            alternative = alternatives.get(decision["candidate_id"], {})
            alternative = dict(alternative, query_partial=genes_by_id.get(decision["candidate_id"], {}).get("partial", "False") in [True, "True"])
            if decision["candidate_id"] in nife_region_proofs:
                alternative["NiFe_core_region"] = nife_region_proofs[decision["candidate_id"]]
            if decision['candidate_id'] in backbones:alternative['NiFe_reference_backbone']=backbones[decision['candidate_id']]
            if decision['candidate_id'] in fe_evidence:alternative['Fe_reference_identity']=fe_evidence[decision['candidate_id']]
            if decision['candidate_id'] in profile_evidence:alternative['FeFe_profile_ligands']=profile_evidence[decision['candidate_id']]
            if decision['candidate_id'] in experimental_backbone_evidence:alternative['FeFe_experimental_backbone']=experimental_backbone_evidence[decision['candidate_id']]
            if decision['candidate_id'] in assay_evidence:alternative['Assay_sequence_identity']=assay_evidence[decision['candidate_id']]
            alternative['Coding_completeness'] = coding_rows.get(decision['candidate_id'])
            status = sequence_status({}, feature, alternative, rules.get(decision["family"]), review,
                                     coding_completeness=alternative['Coding_completeness'], gene=genes_by_id[decision['candidate_id']])
            # No additional-candidate status may inherit high-confidence calibration.
            if status == "calibrated_high" or status not in {"sequence_supported_uncalibrated_core", "sequence_supported_uncalibrated_noncanonical", NIFE_CORE_REGION_STATUS, "sequence_support_center_model_unresolved", "hydrogenase_like_lineage_unconfirmed"}:
                status = "unconfirmed"
            if status == NIFE_CORE_REGION_STATUS:
                decision["NiFe_core_region"] = alternative["NiFe_core_region"]
            decision["sequence_status"] = status
            decision["identity_status"] = status
            decision["features"] = feature
            decision["identity_reasons"] = identity_reasons(feature, rules.get(decision["family"])) if feature else ["identity_features_missing"]
            decision["partial_cds"] = genes_by_id.get(decision["candidate_id"], {}).get("partial", "")
            decision['NiFe_reference_backbone'] = backbones.get(decision['candidate_id'])
            decision['Fe_reference_identity'] = fe_evidence.get(decision['candidate_id'])
            decision['FeFe_profile_ligands'] = profile_evidence.get(decision['candidate_id'])
            decision['FeFe_experimental_backbone'] = experimental_backbone_evidence.get(decision['candidate_id'])
            decision['Assay_sequence_identity'] = assay_evidence.get(decision['candidate_id'])
            decision['Coding_completeness'] = coding_rows.get(decision['candidate_id'])
            evidence[decision["candidate_id"]] = {"annotation": str(work), "features": str(work / "confirmation/features.json"),
                "active_centers": str(work / "active_center_summary.tsv")}
            if status == NIFE_CORE_REGION_STATUS:
                evidence[decision["candidate_id"]]["NiFe_core_region"] = "identity_scope.nife_core_regions"
            if (work / "noncanonical_centers/alternatives.tsv").is_file():
                evidence[decision["candidate_id"]]["noncanonical"] = str(work / "noncanonical_centers/alternatives.tsv")
    for decision in decisions:
        if decision["sequence_status"] == "unconfirmed" and not decision["blocking_reasons"]:
            decision["blocking_reasons"] = ["sequence_core_evidence_insufficient"]
    if source_inputs != {p: file_sha(Path(p)) for p in source_inputs}:
        raise ValueError("Candidate screen source inputs changed during confirmation")
    if any(file_sha(Path(p))!=h for p,h in backbone_inputs.items()):raise ValueError('NiFe backbone input changed')
    if any(file_sha(Path(p))!=h for p,h in fe_inputs.items()):raise ValueError('Fe identity input changed')
    if any(file_sha(Path(p))!=h for p,h in profile_inputs.items()):raise ValueError('FeFe profile ligand input changed')
    if any(file_sha(Path(p))!=h for p,h in experimental_backbone_inputs.items()):raise ValueError('FeFe experimental backbone input changed')
    if any(file_sha(Path(p))!=h for p,h in assay_inputs.items()):raise ValueError('Assayed sequence identity input changed')
    if any(file_sha(Path(p))!=h for p,h in coding_inputs.items()):raise ValueError('Coding completeness input or proof changed')
    if any(file_sha(Path(p))!=h for p,h in nife_region_inputs.items()):raise ValueError('NiFe core-region input or proof changed')
    code_paths = [Path(__file__), ROOT / "research_workflow/explain_candidates.py", ROOT / "research_workflow/candidate_screen.py",
                  ROOT / "research_workflow/identity_scope.py", ROOT / "scripts/active_centers.py",
                  ROOT / "scripts/noncanonical_center_evidence.py", ROOT / "scripts/catalytic_confirmation.py", ROOT / "research_workflow/fefe_profile_ligands.py", ROOT / "research_workflow/fefe_experimental_backbone.py",
                  ROOT / 'research_workflow/coding_completeness.py', ROOT / 'research_workflow/canonical_hmd.py']
    model_paths = [ROOT / "references/processed/catalytic_confirmation/manifest.json", ROOT / "references/processed/active_centers/manifest.json",
                   ROOT / "config/catalytic_confirmation_thresholds.json", ROOT / "config/catalytic_confirmation.json",
                   ROOT / "config/noncanonical_center_evidence.json", ROOT / "research_workflow/identity_scope_source_review.json", ROOT / "research_workflow/fefe_profile_ligands_model.json", ROOT / "research_workflow/fefe_experimental_backbone_model.json"]
    extra_core_region_count = sum(d.get("sequence_status") == "sequence_supported_uncalibrated_core_region" for d in decisions)
    report = dict(status="complete", annotation=str(annotation_dir), screen=str(screen_dir), candidate_count=len(extras),
                  extra_candidate_ids=sorted(d["candidate_id"] for d in decisions), decisions=decisions,
                  evidence=evidence, source_screen_manifest_sha256=file_sha(screen_dir / "candidate_screen_manifest.json"),
                  source_screen_outputs={n: file_sha(screen_dir / n) for n in manifest["output_sha256"]},
                  source_screen_inputs=source_inputs,
                  backbone_input_sha256=backbone_inputs,
                  Fe_input_sha256=fe_inputs,
                  FeFe_profile_input_sha256=profile_inputs,
                  FeFe_experimental_backbone_input_sha256=experimental_backbone_inputs,
                  Assay_sequence_identity_input_sha256=assay_inputs,
                  coding_completeness_input_sha256=coding_inputs,
                  NiFe_core_region_input_sha256=nife_region_inputs,
                  code_sha256={str(p): file_sha(p) for p in code_paths},
                  model_sha256={str(p): file_sha(p) for p in model_paths},
                  scope="Research-only identity evidence; no calibrated_high, function, reaction or physiology assignment.")
    (out_dir / "candidate_confirmation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    (out_dir / "candidate_confirmation_manifest.json").write_text(json.dumps({
        "status": "complete", "report_sha256": file_sha(out_dir / "candidate_confirmation.json"),
        "source_screen_manifest_sha256": report["source_screen_manifest_sha256"],
        "candidate_count": len(extras),
        "extra_confirmed_count": sum(d["sequence_status"] in {"sequence_supported_uncalibrated_core", "sequence_supported_uncalibrated_noncanonical"} for d in decisions),
        "extra_core_region_count": extra_core_region_count
    }, ensure_ascii=False, indent=2) + "\n")
    return report


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    print(json.dumps(confirm(args.annotation, args.screen, args.out, args.threads), ensure_ascii=False))
