"""Report maturation search observations without inferring cofactor assembly."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
VERSION = "maturation_observation_v3"


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def audit(annotation, output):
    """Read bound annotation products and write a separate observation report."""
    annotation, output = Path(annotation).resolve(), Path(output)
    manifest_path = annotation / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    bindings = {"run_manifest.json": _sha(manifest_path)}
    missing = []
    tables = {}
    for name in ("genes.tsv", "domain_hits.tsv", "maturation.tsv"):
        path = annotation / name
        expected = manifest.get("output_sha256", {}).get(name)
        if not path.is_file() or not expected:
            missing.append(name)
            tables[name] = []
            continue
        bindings[name] = _sha(path)
        if bindings[name] != expected:
            raise ValueError(f"Maturation observation input changed: {name}")
        tables[name] = _table(path)
    fasta = annotation / "proteins.faa"
    proteins = {}
    if not fasta.is_file() or not manifest.get("normalized_proteins_sha256"):
        missing.append("proteins.faa")
    else:
        bindings["proteins.faa"] = _sha(fasta)
        if bindings["proteins.faa"] != manifest["normalized_proteins_sha256"]:
            raise ValueError("Maturation observation input changed: proteins.faa")
        proteins = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    genes = tables["genes.tsv"]
    ids = [row["id"] for row in genes]
    if not missing and (len(ids) != len(set(ids)) or set(ids) != proteins.keys()):
        raise ValueError("Maturation genes and proteins have inconsistent ownership")
    if any(row.get("protein_id") not in proteins for name in ("domain_hits.tsv", "maturation.tsv")
           for row in tables[name]) and not missing:
        raise ValueError("Maturation evidence belongs to an unknown protein")

    paths = {"policy": ROOT / "config/maturation_domain_policy.json",
             "audit_summary": ROOT / "references/processed/maturation_domain_audit/summary.json",
             "model": ROOT / "references/processed/maturation_domain_audit/maturation_domains.hmm"}
    sources = {key: {"path": str(path.relative_to(ROOT)), "sha256": _sha(path)} for key, path in paths.items()}
    policy = json.loads(paths["policy"].read_text())
    summary = json.loads(paths["audit_summary"].read_text())
    if (summary.get("policy_sha256") != sources["policy"]["sha256"]
            or summary.get("combined_model_sha256") != sources["model"]["sha256"]):
        raise ValueError("Maturation model does not match its audit summary")
    for name, actual in (("maturation_domain_policy_sha256", sources["policy"]["sha256"]),
                         ("maturation_domain_model_sha256", sources["model"]["sha256"])):
        if name not in manifest:
            missing.append(name)
        elif manifest[name] != actual:
            raise ValueError(f"Maturation annotation used a different model: {name}")

    log = annotation / "maturation_domains/hmmsearch.log"
    text = log.read_text() if log.is_file() else ""
    target = re.search(r"^# target sequence database:\s*(.+)$", text, re.MULTILINE)
    model = re.search(r"^# query HMM file:\s*(.+)$", text, re.MULTILINE)
    target_bound = bool(target and Path(target[1].strip()).resolve() == fasta)
    model_bound = bool(model and Path(model[1].strip()).is_file()
                       and _sha(Path(model[1].strip())) == sources["model"]["sha256"])
    complete = bool(text.rstrip().endswith("\n[ok]") and target_bound and model_bound
                    and not missing and manifest.get("status") == "complete")

    mode = manifest.get("input_mode", "unknown")
    groups = defaultdict(list)
    for gene in genes:
        if mode == "genome":
            scope = "input_genome"
        elif mode == "bin" and gene.get("bin_id"):
            scope = "bin:" + gene["bin_id"]
        else:
            scope = "contig:" + gene.get("contig", "unavailable")
        groups[scope].append(gene)
    domains = defaultdict(set)
    for row in tables["domain_hits.tsv"]:
        if row.get("passes") == "True":
            domains[row["protein_id"]].add(row["accession"].split(".")[0])
    matches = defaultdict(set)
    for row in tables["maturation.tsv"]:
        matches[row["maturation_gene"]].add(row["protein_id"])
    rows = []
    for scope, members in sorted(groups.items()):
        member_ids = {gene["id"] for gene in members}
        unavailable = any(gene.get("coordinate_evidence") != "genomic_CDS" for gene in members)
        open_boundary = any(gene.get("partial") == "True" or gene.get("topology") != "circular" for gene in members)
        observability = ("coordinates_unavailable" if unavailable else
                         "open_contig_boundaries" if open_boundary else "provided_sequences_only")
        for role, rule in policy["rules"].items():
            hits = sorted(matches[role] & member_ids)
            architecture = sorted(pid for pid in member_ids
                                  if set(rule.get("all_of", ())) <= domains[pid]
                                  and (not rule.get("any_of") or set(rule["any_of"]) & domains[pid]))
            supported = sorted(set(hits) & set(architecture))
            rows.append(dict(scope_id=scope, maturation_role=role,
                             search_status=("not_evaluated" if not complete else "detected" if hits else "searched_not_detected"),
                             evidence_status=("profile_architecture_supported" if supported else "homology_only" if hits else "none"),
                             genomic_observability=observability, searched_protein_count=len(member_ids),
                             matched_protein_ids=hits, profile_architecture_protein_ids=supported,
                             profile_only_protein_ids=sorted(set(architecture) - set(hits)),
                             specificity_scope=rule["scope"]))
    receipt = dict(schema=VERSION, status="complete" if complete else "not_evaluated",
                   annotation=str(annotation), input_mode=mode, missing_inputs=missing,
                   input_sha256=bindings, source_bindings=sources,
                   search_proof=dict(complete=complete, target_bound=target_bound, model_bound=model_bound,
                                     log_sha256=_sha(log) if log.is_file() else None), rows=rows,
                   family_rule_scope={"NiFe": dict(status="modeled_factors", roles=[r for r in policy["rules"] if r.startswith("hyp")]),
                                      "FeFe": dict(status="modeled_factors", roles=[r for r in policy["rules"] if r.startswith("hyd")]),
                                      "Fe_Hmd": dict(status="not_applicable", roles=[],
                                                     reason="Covered NiFe/FeFe models do not assess FeGP maturation.")},
                   scope_note="Search is limited to supplied predicted proteins. No biological absence, mandatory maturation checklist, core ownership, assembly capacity or cofactor occupancy is inferred.")
    output.mkdir(parents=True, exist_ok=True)
    (output / "maturation_observation.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.annotation, args.out)
    print(json.dumps(dict(status=result["status"], rows=len(result["rows"]))))
