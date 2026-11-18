#!/usr/bin/env python3
"""Evaluate development/calibration identity data, keeping final validation reserved."""
import argparse
import collections
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from annotate import ROOT, write_tsv
from run_acceptance_checks import state_hashes
from validate_fes_prediction import wilson
from discovery_evidence import stage_fingerprint


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def identity_metrics(truth, predicted):
    tp = sum(r["label"] == "positive" and r["unit_id"] in predicted for r in truth)
    fn = sum(r["label"] == "positive" and r["unit_id"] not in predicted for r in truth)
    fp = sum(r["label"] == "negative" and r["unit_id"] in predicted for r in truth)
    tn = sum(r["label"] == "negative" and r["unit_id"] not in predicted for r in truth)
    return dict(TP=tp, FP=fp, FN=fn, TN=tn, recall=tp/(tp+fn) if tp+fn else None,
                precision=tp/(tp+fp) if tp+fp else None, recall_wilson95=wilson(tp, tp+fn),
                precision_wilson95=wilson(tp, tp+fp), negative_hit_rate=fp/(fp+tn) if fp+tn else None)


def family(row):
    values = set()
    for category in row["categories"]:
        for name in ["NiFe", "FeFe", "Fe"]:
            if category == name or category.startswith("["+name+"]"):
                values.add(name)
    return ";".join(sorted(values))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--split", choices=["development", "calibration", "final_validation"], default="development")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument('--training-audit',type=Path,help='Optional source-only core homology audit frozen before evaluating this partition')
    parser.add_argument('--discovery-lock',type=Path,help='Required for final validation: fixed candidate-discovery state and completed calibration')
    args = parser.parse_args()
    if args.split=='final_validation' and (not args.training_audit or not args.discovery_lock):
        parser.error('Final candidate-discovery validation requires --training-audit and --discovery-lock')
    corpus, out = args.corpus.resolve(), args.out.resolve()
    if (out/"summary.json").exists():
        parser.error("Use a fresh evaluation directory")
    frozen = json.loads((corpus/"freeze.json").read_text())
    if frozen["status"] != "frozen_before_first_prediction":
        raise ValueError("Expected frozen partitions")
    for filename, expected in frozen["artifacts"].items():
        if sha(corpus/filename) != expected:
            raise ValueError(f"Corpus artifact changed: {filename}")
    from identity_truth import load_units
    all_rows = load_units(corpus)
    rows = [r for r in all_rows if r["split"] == args.split]
    profiles=json.loads((ROOT/'references/processed/profile_manifest.json').read_text())
    training_ids={uid for profile in profiles for uid in profile.get('training_corpus_unit_ids',[])}
    confirmation_prototypes=json.loads((ROOT/'references/processed/catalytic_confirmation/prototypes.json').read_text())
    confirmation_training_ids={s['unit_id'] for r in confirmation_prototypes for s in r['sources'] if s.get('unit_id')}
    confirmation_training_ids.update(s['identity_review']['source_unit_id'] for r in confirmation_prototypes for s in r['sources'] if s.get('identity_review'))
    audit=None;independent_ids=None;audited_remote_ids=set();audit_sha_before=None
    if args.training_audit:
        audit=json.loads((args.training_audit/'summary.json').read_text())
        if audit['status']!='source_audit_frozen_before_partition_prediction' or audit['corpus_freeze_sha256']!=sha(corpus/'freeze.json'):
            raise ValueError('Training core audit is not frozen for this corpus')
        if audit['profile_manifest_sha256']!=sha(ROOT/'references/processed/profile_manifest.json'):
            raise ValueError('Training core audit is stale for current profiles')
        audit_sha_before=sha(args.training_audit/'summary.json')
        for name,expected in audit['artifacts'].items():
            if sha(args.training_audit/name)!=expected:raise ValueError('Training core audit artifact changed: '+name)
        independent_ids=set(audit['independent_unit_ids'])
        audited_remote_ids={r['unit_id'] for r in table(args.training_audit/'per_unit_audit.tsv') if r['remote_Group4_after_core_audit']=='True'}
    query = corpus/(args.split+".faa")
    if not rows or {r.id for r in SeqIO.parse(query, "fasta")} != {r["unit_id"] for r in rows}:
        raise ValueError("Empty or inconsistent partition FASTA")
    if args.split!='final_validation':
        assert not {r["unit_id"] for r in rows} & {r["unit_id"] for r in all_rows if r["split"] == "final_validation"}
    out.mkdir(parents=True, exist_ok=True)
    before = state_hashes()
    if args.split=='final_validation':
        lock=json.loads(args.discovery_lock.read_text())
        calibration_path=Path(lock['calibration_summary']);calibration=json.loads(calibration_path.read_text())
        if (lock['stage']!='candidate_discovery_only' or lock['model_state_hashes']!=before or lock['discovery_stage_fingerprint']!=stage_fingerprint(ROOT) or
            lock['corpus_freeze_sha256']!=sha(corpus/'freeze.json') or lock['training_audit_summary_sha256']!=audit_sha_before or
            lock['calibration_summary_sha256']!=sha(calibration_path) or calibration['split']!='calibration' or calibration['model_state_hashes']!=before or
            lock['acceptance_targets_sha256']!=sha(ROOT/'config/acceptance_targets.json')):
            raise ValueError('Candidate discovery lock/calibration is missing, changed or stale')
    command = [sys.executable, str(ROOT/"scripts/annotate.py"), "--proteins", str(query), "--mode", "transcript", "--threads", "4", "--out", str(out/"annotation")]
    with (out/"annotation.log").open("w") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
    anchors = {r["protein_id"]: r for r in table(out/"annotation/catalytic_anchors.tsv")}
    hmm = table(out/"annotation/family_candidates.tsv")
    hits = collections.defaultdict(list)
    for line in (out/"annotation/family.domtbl").read_text().splitlines():
        if line.startswith("#"):
            continue
        c = line.split()
        hits[c[0]].append(dict(profile=c[3], evalue=float(c[12]), score=float(c[13]),
                              coverage=(int(c[16])-int(c[15])+1)/int(c[5]), model_start=int(c[15]), model_end=int(c[16])))
    rejected = {r["protein_id"]: r for r in table(out/"annotation/rejected.tsv")}
    sets = dict(hmm_discovery={r["protein_id"] for r in hmm}, combined_candidate=set(anchors),
                supported_catalytic={k for k, r in anchors.items() if r["catalytic_status"] == "supported"})
    sets['calibrated_high_identity']={k for k,r in anchors.items() if r.get('catalytic_identity_confidence')=='calibrated_high'}
    hmm_families=collections.defaultdict(set)
    for hit in hmm:hmm_families[hit['protein_id']].add(hit['family'])
    sets['hmm_correct_family_discovery']={r['unit_id'] for r in rows if r['unit_id'] in sets['hmm_discovery'] and (r['label']=='negative' or family(r) in hmm_families[r['unit_id']])}
    sets['group4_profile_discovery']={r['protein_id'] for r in hmm if r['profile']=='NiFe_Group4'}
    predictions = []
    thresholds = json.loads((ROOT/"config/thresholds.json").read_text())
    for row in rows:
        uid = row["unit_id"]
        domains = sorted(hits[uid], key=lambda h: h["score"], reverse=True)
        reasons = []
        if row["label"] == "positive" and uid not in sets["hmm_discovery"]:
            if not domains:
                reasons.append("no_reported_HMM_domain")
            elif all(h["evalue"] > thresholds["hmm_domain_evalue"] for h in domains):
                reasons.append("HMM_evalue_above_discovery_threshold")
            else:
                reasons.append("significant_HMM_domains_below_model_coverage_threshold")
        if uid in rejected:
            reasons.append(rejected[uid]["reason"]+":"+rejected[uid]["reference"])
        anchor = anchors.get(uid, {})
        predictions.append(dict(unit_id=uid, truth=row["label"], frozen_truth=row.get('frozen_label',row['label']), truth_family=family(row), categories=";".join(row["categories"]),
            evidence_tiers=";".join(row["evidence_tiers"]), accessions=";".join(sorted({s["accession"] for s in row["sources"]})),
            genus_names=";".join(row["genus_names"]), homolog_component=row["homolog_component"],
            used_for_discovery_training=uid in training_ids,independent_after_core_audit=uid in independent_ids if independent_ids is not None else '',
            used_for_confirmation_training=uid in confirmation_training_ids,
            max_comparable_training_identity=row["max_comparable_training_identity"], remote_group4=row["remote_group4"],
            hmm_discovery=uid in sets["hmm_discovery"], combined_candidate=uid in sets["combined_candidate"],
            hmm_families=';'.join(sorted(hmm_families[uid])),hmm_correct_family_discovery=uid in sets['hmm_correct_family_discovery'],
            group4_profile_discovery=uid in sets['group4_profile_discovery'],
            supported_catalytic=uid in sets["supported_catalytic"], predicted_family=anchor.get("catalytic_family", ""),
            calibrated_high_identity=uid in sets['calibrated_high_identity'],identity_confidence=anchor.get('catalytic_identity_confidence','not_output'),
            catalytic_status=anchor.get("catalytic_status", "not_output"), structure_status=anchor.get("structure_status", "not_output"),
            function_status=anchor.get("function_status", "not_output"), function_ids=anchor.get("potential_function_ids", ""),
            diagnostic_reasons=";".join(reasons), hmm_domains=json.dumps(domains)))
    metrics = {}
    # Deterministic component representatives are not selected for successful predictions.
    representatives = {}
    for row in sorted(rows, key=lambda r: r["unit_id"]):
        representatives.setdefault((row["homolog_component"], row["label"], family(row)), row)
    strata = {"all": rows, "component_representatives": list(representatives.values())}
    strata['not_used_for_discovery_training']=[r for r in rows if r['unit_id'] not in training_ids]
    strata['not_used_for_confirmation_training']=[r for r in rows if r['unit_id'] not in confirmation_training_ids]
    strata['not_used_for_either_identity_training']=[r for r in rows if r['unit_id'] not in training_ids|confirmation_training_ids]
    if independent_ids is not None:
        strata['independent_after_core_audit']=[r for r in rows if r['unit_id'] in independent_ids]
        for name in ['NiFe','FeFe','Fe']:
            strata['independent_after_core_audit_family_'+name]=[r for r in rows if r['unit_id'] in independent_ids and r['label']=='positive' and family(r)==name]
        strata['remote_Group4_after_core_audit']=[r for r in rows if r['unit_id'] in audited_remote_ids]
    for name in ["NiFe", "FeFe", "Fe"]:
        strata["positive_family_"+name] = [r for r in rows if r["label"] == "positive" and family(r) == name]
    strata["remote_Group4"] = [r for r in rows if r["label"] == "positive" and r["remote_group4"]]
    for category in sorted({c for r in rows if r["label"] == "negative" for c in r["categories"]}):
        strata["negative_"+category] = [r for r in rows if r["label"] == "negative" and category in r["categories"]]
    for tier in sorted({t for r in rows for t in r["evidence_tiers"]}):
        strata["evidence_"+tier] = [r for r in rows if tier in r["evidence_tiers"]]
    for stratum, selected in strata.items():
        metrics[stratum] = {stage: identity_metrics(selected, predicted) for stage, predicted in sets.items()}
    after = state_hashes()
    if before != after:
        raise ValueError("Code, policy or models changed during evaluation")
    if args.training_audit and sha(args.training_audit/'summary.json')!=audit_sha_before:raise ValueError('Core audit changed during evaluation')
    write_tsv(out/"per_unit.tsv", predictions)
    write_tsv(out/"metrics.tsv", [dict(stratum=s, stage=stage, **m) for s, stages in metrics.items() for stage, m in stages.items()])
    report = dict(status={"development":"development_diagnostic_complete","calibration":"calibration_diagnostic_complete","final_validation":"final_discovery_validation_complete"}[args.split],
        split=args.split, corpus=str(corpus), corpus_freeze_sha256=sha(corpus/"freeze.json"), command=command,
        model_state_hashes=before, after_state_hashes=after, units=len(rows), labels=dict(collections.Counter(r["label"] for r in rows)),
        discovery_stage_fingerprint=stage_fingerprint(ROOT),
        homolog_components=len({r["homolog_component"] for r in rows}), genera=len({g for r in rows for g in r["genus_taxids"]}),
        training_unit_ids=sorted(training_ids),confirmation_training_unit_ids=sorted(confirmation_training_ids),
        core_audit_independence_scope='Candidate-discovery training only; confirmation-training overlap is separately explicit.',
        training_audit=str(args.training_audit.resolve()) if args.training_audit else None,
        training_audit_summary_sha256=sha(args.training_audit/'summary.json') if args.training_audit else None,
        discovery_lock=str(args.discovery_lock.resolve()) if args.discovery_lock else None,
        discovery_lock_sha256=sha(args.discovery_lock) if args.discovery_lock else None,
        independent_positive_families={name:len(strata.get('independent_after_core_audit_family_'+name,[])) for name in ['NiFe','FeFe','Fe']},
        audited_remote_Group4_metadata=dict(units=len(strata.get('remote_Group4_after_core_audit',[])),
            branches=sorted({c for r in strata.get('remote_Group4_after_core_audit',[]) for c in r['categories'] if c.startswith('[NiFe] Group 4')}),
            genera=sorted({g for r in strata.get('remote_Group4_after_core_audit',[]) for g in r['genus_names']})),
        metrics=metrics, wrong_family_supported=[r["unit_id"] for r in predictions if r["truth"] == "positive" and r["supported_catalytic"] and r["truth_family"] != r["predicted_family"]],
        negative_function_outputs=[r["unit_id"] for r in predictions if r["truth"] == "negative" and r["function_ids"]],
        artifacts={"per_unit.tsv": sha(out/"per_unit.tsv"), "metrics.tsv": sha(out/"metrics.tsv"), "run_manifest.json": sha(out/"annotation/run_manifest.json")},
        limitations=["Final_validation applies only to the locked candidate-discovery stage; other statuses are exploratory and do not validate complex/function performance." if args.split=='final_validation' else "These are development/calibration results, not final blind-test acceptance.",
            "HydDB author-curated identity and UniProt reviewed identity are not per-protein experimental activity.",
            "Candidate false positives are measured separately from supported catalytic identity; supported is an uncalibrated current status, not a newly validated high-confidence label.",
            "Confirmation prototypes include development identities; in-training confirmation metrics are not independent accuracy. Use leave-component/genus-out development features and separately locked calibration/final confirmation evidence.",
            "Transcript mode prevents artificial operons; it cannot evaluate complete-complex or molecular-function sensitivity.",
            "Sequence-level Wilson intervals are descriptive because homologs remain correlated; deterministic component-representative metrics are also reported.",
            "Final candidate-discovery results must not be reused to tune that stage." if args.split=='final_validation' else "Final validation sequences were not annotated."])
    (out/"summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({k: report[k] for k in ["status", "units", "labels", "homolog_components", "genera", "negative_function_outputs"]}))
    print(json.dumps(report["metrics"]["all"]))


if __name__ == "__main__":
    main()
