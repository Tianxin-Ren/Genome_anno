"""Audit an experimentally observed remote partner and evaluate the unchanged predictor."""
import argparse
import copy
import csv
import hashlib
import json
import subprocess
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT, annotate, extract_genbank
from catalytic_confirmation import sha, table
from complex_extensions import load_contracts
from prepare_natural_complexes import anonymous_genes
from run_acceptance_checks import state_hashes


def exposure_audit(out, native, policy):
    selected = {g['id']: g for g in native if g['id'] in policy['observed_members'].values()}
    fasta = out / 'observed_members.faa'
    fasta.write_text(''.join('>' + k + '\n' + g['sequence'] + '\n' for k, g in sorted(selected.items())))
    files = [ROOT / 'references/processed' / name for name in
             ['NiFe.faa', 'FeFe.faa', 'Fe.faa', 'NiFe_Group4.faa', 'NiFe_Group4.core.faa',
              'FeFe_diverse_core.faa', 'FeFe_diverse_core.core.faa', 'references.faa',
              'complex_extensions/references.faa', 'catalytic_confirmation/prototypes.faa']]
    targets = {}
    for path in files:
        for record in SeqIO.parse(path, 'fasta'):
            targets['target_' + str(len(targets) + 1)] = dict(source_file=str(path.relative_to(ROOT)),
                source_id=record.id, sequence=str(record.seq))
    subject = out / 'exposed_sequences.faa'
    subject.write_text(''.join('>' + k + '\n' + r['sequence'] + '\n' for k, r in targets.items()))
    columns = 'qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore'.split()
    command = [str(ROOT / '.env/bin/blastp'), '-query', str(fasta), '-subject', str(subject),
               '-evalue', '1e-5', '-max_hsps', '1', '-outfmt', '6 ' + ' '.join(columns),
               '-out', str(out / 'source_overlap.tsv')]
    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with (out / 'source_overlap.tsv').open() as handle:
        hits = list(csv.DictReader(handle, fieldnames=columns, delimiter='\t'))
    units = json.loads((ROOT / 'references/validation/round09_identity_corpus/units.json').read_text())
    rows = []
    for locus, gene in sorted(selected.items()):
        sequence = gene['sequence']; digest = hashlib.sha256(sequence.encode()).hexdigest()
        related = []
        for hit in hits:
            if hit['qseqid'] != locus:
                continue
            qcov = (int(hit['qend']) - int(hit['qstart']) + 1) / int(hit['qlen'])
            scov = (int(hit['send']) - int(hit['sstart']) + 1) / int(hit['slen'])
            shorter = qcov if int(hit['qlen']) <= int(hit['slen']) else scov
            if float(hit['pident']) >= policy['homolog_identity_percent'] and shorter >= policy['homolog_shorter_sequence_coverage']:
                target = targets[hit['sseqid']]
                related.append(dict(source_file=target['source_file'], source_id=target['source_id'],
                    identity_percent=float(hit['pident']), query_coverage=qcov, subject_coverage=scov))
        rows.append(dict(locus=locus, protein_accession=gene['protein_accession'], sequence_sha256=digest,
            exact_exposed_targets=[{k: t[k] for k in ['source_file', 'source_id']} for t in targets.values() if t['sequence'] == sequence],
            homolog_overlaps=related, identity_corpus_exact=[dict(unit_id=u['unit_id'], split=u['split'],
                homolog_component=u['homolog_component']) for u in units if u['sequence_sha256'] == digest]))
    report = dict(status='source_exposure_audited_before_prediction', independent_final_eligible=False,
                  rows=rows, command=command, input_files={str(p.relative_to(ROOT)): sha(p) for p in files},
                  artifacts={p.name: sha(p) for p in [fasta, subject, out / 'source_overlap.tsv']},
                  limitation='Homolog overlap rejects final eligibility; this assay source remains development-only. No taxonomy-independent or atomic-coordination performance claim.')
    (out / 'source_exposure.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(); out = args.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    path = ROOT / 'config/round30_distributed_function.json'; policy = json.loads(path.read_text())
    for name, digest in policy['source_files'].items():
        if sha(ROOT / name) != digest:
            raise ValueError('Experimental source changed: ' + name)
    before = state_hashes(); native = extract_genbank(ROOT / policy['genome'])[0]
    by_id = {g['id']: g for g in native}
    for locus, accession in policy['protein_accessions'].items():
        if by_id[locus]['protein_accession'] != accession:
            raise ValueError('Exact native protein accession mismatch: ' + locus)
    audit = exposure_audit(out, native, policy)
    genes, mapping = anonymous_genes(native)
    observed = {mapping[n] for n in policy['observed_members'].values()}
    hyd = {mapping[n] for n in policy['hydrogenase_loci']}; anchor_id = mapping[policy['anchor']]
    freeze = dict(status='development_truth_frozen_before_native_prediction', protocol_sha256=sha(path),
                  source_exposure_sha256=sha(out / 'source_exposure.json'), source_to_query=mapping,
                  proteins=len(genes), observed_purified_members=sorted(observed), independent_final_eligible=False,
                  before_hashes=before, limits=policy['limits'])
    (out / 'freeze.json').write_text(json.dumps(freeze, indent=2) + '\n')
    contracts = {c['contract_id']: c['reaction_key'] for c in load_contracts()}; results = []
    execution = dict(status='running', before_hashes=before, runs=results)
    save = lambda: (out / 'execution.json').write_text(json.dumps(execution, indent=2) + '\n')
    save()
    for mode in policy['modes']:
        altered = copy.deepcopy(genes)
        if mode in ['bin', 'cross_bin']:
            for g in altered:
                g['bin_id'] = 'bin_one'
                if mode == 'cross_bin' and g['id'] == mapping[policy['observed_members']['FdhA']]:
                    g.update(bin_id='bin_two', contig='other_bin_contig')
        run = out / mode; run.mkdir()
        annotations = annotate(altered, run, threads=2, mode='bin' if mode == 'cross_bin' else mode,
                               input_files=[ROOT / policy['genome']])
        anchor = next((a for a in annotations if a['protein_id'] == anchor_id), {})
        members = {r['protein_id'] for r in table(run / 'complex_members.tsv') if r['anchor'] == anchor_id}
        assigned = {contracts[f] for f in anchor.get('potential_function_ids', '').split(';') if f}
        positive, negative = set(policy['positive_reactions']), set(policy['negative_reactions'])
        correct = {k for k in assigned & positive if (hyd if k == 'H2_NAD_Fd_coupled_reduction' else observed) <= members}
        hypotheses = []
        for h in table(run / 'complex_hypotheses.tsv'):
            if h['anchor'] != anchor_id:
                continue
            hs = {k: h[k] for k in ['structure_id', 'composition_complete', 'architecture_status',
                  'anchor_competition_guard', 'decision', 'missing_structure', 'required_member_ids']}
            hs['observed_member_coverage'] = len(set(json.loads(h['required_member_ids'])) & observed)
            hypotheses.append(hs)
        row = dict(mode=mode, anchor=anchor, observed_member_coverage=len(members & observed),
                   missing_observed_members=sorted(observed - members), inferred_extra_members=sorted(members - observed),
                   hydrogenase_members_covered=hyd <= members, distributed_partner_covered=observed <= members,
                   assigned_reactions=sorted(assigned), reaction=dict(TP=len(correct), FP=len(assigned & negative),
                   FN=len(positive - correct), precision=len(correct)/(len(correct)+len(assigned & negative))
                   if correct or assigned & negative else None,
                   automatic_assignment_coverage=len(assigned & positive)/len(positive)),
                   candidate_hypotheses=hypotheses, ungraded_reactions=sorted(assigned - positive - negative),
                   run_manifest_sha256=sha(run / 'run_manifest.json'))
        results.append(row); save(); print(mode, row['observed_member_coverage'], row['reaction'], flush=True)
    after = state_hashes()
    if before != after:
        raise ValueError('Source/model state changed while evaluating the case')
    for name, digest in audit['input_files'].items():
        if sha(ROOT / name) != digest:
            raise ValueError('Source exposure audit input changed: ' + name)
    boundary = next(r for r in results if r['mode'] == 'cross_bin')
    cross_bin_ok = not (set(boundary['assigned_reactions']) & (set(policy['positive_reactions']) - {'H2_NAD_Fd_coupled_reduction'}))
    summary = dict(status='development_distributed_function_evaluated', independent_validation=False,
                   native_proteins=len(genes), source_exposure_sha256=sha(out / 'source_exposure.json'),
                   unchanged_predictor=before == after, rows=results, cross_bin_boundary_passed=cross_bin_ok,
                   policy_sha256=sha(path), limits=policy['limits'])
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    execution.update(status='completed', after_hashes=after, summary_sha256=sha(out / 'summary.json')); save()
    if not cross_bin_ok or any(r['reaction']['FP'] for r in results):
        raise SystemExit('Assay-scoped negative function or cross-bin function was assigned')


if __name__ == '__main__':
    main()
