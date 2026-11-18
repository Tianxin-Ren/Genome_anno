#!/usr/bin/env python3
"""Research-only Fe-S refinement runner.

This packages the reviewed 34-site core plus the provisional 6BTM sites.  It
uses the production mapper and its site-window rescue unchanged; outputs are
structural Fe-S evidence only and are never hydrogenase/function calls.
"""
import argparse, hashlib, json, shutil, subprocess, sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts' / 'fes_sites.py'
MODEL = ROOT / 'research_workflow' / 'fes_refinement_model.json'
MODEL_DIR = ROOT / 'references' / 'processed' / 'fes_refinement'
CORE_DIR = ROOT / 'references' / 'processed' / 'fes'
sys.path.insert(0, str(ROOT / 'scripts'))
from annotate import BLAST_FIELDS, load_blast
from fes_sites import SITE_FIELDS, map_hits, site_evidence, rescue_with_site_windows, full_reference_coordinates, write_tsv


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _check_model():
    spec = json.loads(MODEL.read_text())
    core = json.loads((CORE_DIR / 'sites.json').read_text())
    sites = json.loads((MODEL_DIR / 'sites.json').read_text())
    if len(core) != spec['core_site_count'] or len(sites) != spec['site_count']:
        raise ValueError('unexpected Fe-S model site count')
    if sites[:len(core)] != core:
        raise ValueError('refinement model changed a core Fe-S record')
    if not all(s.get('pdb_id') != '6BTM' or s.get('cofactor_occupancy', s.get('query_cofactor_occupancy', 'unknown')) == 'unknown' for s in sites):
        raise ValueError('6BTM occupancy must remain unknown')
    return spec, sites


def refine(annotation_dir, out_dir, threads=2):
    """Run structural Fe-S refinement on an annotation directory.

    ``annotation_dir`` must contain ``proteins.faa``.  The returned manifest
    records all input/model hashes for reproducibility.
    """
    spec, sites = _check_model()
    annotation_dir, out_dir = Path(annotation_dir), Path(out_dir)
    query = annotation_dir / 'proteins.faa'
    if not query.exists():
        raise FileNotFoundError(query)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_sha = sha(query)
    db = out_dir / 'refinement_seed'
    raw = out_dir / 'fes_blast.tsv'
    log = out_dir / 'fes_blast.log'
    subprocess.run(['makeblastdb', '-blastdb_version', '4', '-in', str(MODEL_DIR/'templates.faa'), '-dbtype', 'prot', '-out', str(db)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with log.open('w') as handle:
        subprocess.run(['blastp', '-query', str(query), '-db', str(db), '-evalue', str(spec['parameters']['max_evalue']), '-max_target_seqs', '1000', '-max_hsps', '10', '-num_threads', str(threads), '-outfmt', '6 ' + ' '.join(BLAST_FIELDS), '-out', str(raw)], check=True, stdout=subprocess.DEVNULL, stderr=handle)
    thresholds = dict(blast_evalue=spec['parameters']['max_evalue'], min_identity=0, min_query_coverage=0, min_reference_coverage=0)
    rows, evidence = map_hits(load_blast(raw, thresholds), sites, spec['parameters'])
    for row in evidence:
        row['alignment_source'] = 'full_reference_search'

    windows = json.loads((MODEL_DIR/'site_window_index.json').read_text())
    window_db = out_dir / 'refinement_site_seed'
    window_raw = out_dir / 'fes_site_window_blast.tsv'
    subprocess.run(['makeblastdb', '-blastdb_version', '4', '-in', str(MODEL_DIR/'site_windows.faa'), '-dbtype', 'prot', '-out', str(window_db)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with (out_dir/'fes_site_window_blast.log').open('w') as handle:
        subprocess.run(['blastp', '-query', str(query), '-db', str(window_db), '-evalue', str(spec['parameters']['max_evalue']), '-max_target_seqs', '1000', '-max_hsps', '10', '-num_threads', str(threads), '-outfmt', '6 ' + ' '.join(BLAST_FIELDS), '-out', str(window_raw)], check=True, stdout=subprocess.DEVNULL, stderr=handle)
    by_site = {s['site_id']: s for s in sites}
    local = []
    for hit in load_blast(window_raw, thresholds):
        site = by_site[windows[hit['sseqid']]]
        local.append(dict(site_evidence(full_reference_coordinates(hit, site), site, spec['parameters']), alignment_source='reviewed_site_window_search'))
    rows = rescue_with_site_windows(rows, local)
    evidence.extend(local)
    # These rows are deliberately structural and unlinked; no H2 identity is inferred.
    for row in rows:
        row['complex_anchors'] = []
        row['association_scope'] = 'unlinked_FeS_homolog_not_hydrogenase_identity'
        row['cofactor_occupancy'] = 'unknown'
    def serial(row):
        return {k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v for k, v in row.items()}
    write_tsv(out_dir/'fes_sites.tsv', [serial(r) for r in rows], SITE_FIELDS)
    write_tsv(out_dir/'fes_template_evidence.tsv', [serial(r) for r in evidence], list(evidence[0]) if evidence else ['protein_id', 'template_site_id', 'status'])
    ligands = [dict(fes_site_id=r['fes_site_id'], protein_id=r['protein_id'], template_site_id=r['template_site_id'], **ligand) for r in rows for ligand in r['ligand_residues']]
    write_tsv(out_dir/'fes_ligands.tsv', ligands, ['fes_site_id','protein_id','template_site_id','reference_position','expected_residue','query_position','observed_residue','status','coverage_reason'])
    manifest = {'status':'research_refinement_complete', 'input_annotation_dir':str(annotation_dir), 'query_fasta_sha256':input_sha, 'model':str(MODEL.relative_to(ROOT)), 'model_sha256':sha(MODEL), 'core_manifest_sha256':spec['core_manifest_sha256'], 'fes_sites_script_sha256':spec['fes_sites_script_sha256'], 'parameters':spec['parameters'], 'site_count':len(sites), 'functional_scope':spec['functional_scope'], 'outputs':{p.name:sha(p) for p in [out_dir/'fes_sites.tsv',out_dir/'fes_ligands.tsv',out_dir/'fes_template_evidence.tsv']}}
    (out_dir/'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('annotation_dir')
    ap.add_argument('out_dir')
    ap.add_argument('--threads', type=int, default=2)
    args = ap.parse_args()
    print(json.dumps(refine(args.annotation_dir, args.out_dir, args.threads), indent=2))


if __name__ == '__main__':
    main()
