#!/usr/bin/env python3
"""Audit confirmation prototypes against held-out full proteins and author cores before prediction."""
import argparse
import collections
import json
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT,write_tsv
from catalytic_confirmation import sha,table


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args();out=args.out.resolve()
    if (out/'summary.json').exists():parser.error('Use a fresh source audit directory')
    out.mkdir(parents=True,exist_ok=True);policy=json.loads((ROOT/'config/catalytic_confirmation.json').read_text());corpus=ROOT/policy['corpus'];old=ROOT/policy['training_audit']
    freeze=json.loads((corpus/'freeze.json').read_text());previous=json.loads((old/'summary.json').read_text());model=ROOT/'references/processed/catalytic_confirmation'
    if previous['corpus_freeze_sha256']!=sha(corpus/'freeze.json'):raise ValueError('Confirmation source audit is stale')
    for name,digest in freeze['artifacts'].items():
        if sha(corpus/name)!=digest:raise ValueError('Confirmation corpus changed')
    for name,digest in previous['artifacts'].items():
        if sha(old/name)!=digest:raise ValueError('Original core audit changed')
    from identity_truth import load_units
    ip=json.loads((ROOT/'config/identity_corpus_policy.json').read_text());units=load_units(corpus);byid={r['unit_id']:r for r in units};prototypes=json.loads((model/'prototypes.json').read_text())
    candidates=[r for r in units if r['split'] in ['calibration','final_validation']];candidate_ids={r['unit_id'] for r in candidates}
    before=dict(corpus_freeze_sha256=sha(corpus/'freeze.json'),confirmation_model_manifest_sha256=sha(model/'manifest.json'),old_core_audit_sha256=sha(old/'summary.json'),
        source_policy='Preserve prior core exclusions; exclude prototype genus overlaps and >=70% identity over >=70% shorter protein/core. Propagate new overlap to existing refined homolog components. No predictions or success labels select independence.')
    (out/'policy_before_alignment.json').write_text(json.dumps(before,indent=2)+'\n')
    training_genera={g for r in prototypes for g in r['genus_taxids']};refined={r['unit_id']:r['refined_component'] for r in table(old/'per_unit_audit.tsv')}
    fragments=[r for r in json.loads((old/'author_fragments.json').read_text()) if candidate_ids&set(r['unit_ids'])]
    query_records=[SeqRecord(Seq(r['sequence'].replace('U','C')),id=r['unit_id'],description='full_source') for r in candidates]
    query_records.extend(SeqRecord(Seq(r['sequence'].replace('U','C')),id=r['fragment_id'],description='author_core_source') for r in fragments)
    links={r['unit_id']:[r['unit_id']] for r in candidates};links.update({r['fragment_id']:list(candidate_ids&set(r['unit_ids'])) for r in fragments})
    SeqIO.write(query_records,out/'source_queries.faa','fasta');SeqIO.write([SeqRecord(Seq(r['sequence'].replace('U','C')),id=r['prototype_id'],description='leakage_normalization_only') for r in prototypes],out/'prototype_sources.faa','fasta')
    with (out/'alignment.log').open('w') as log:
        subprocess.run(['makeblastdb','-in',str(out/'prototype_sources.faa'),'-dbtype','prot','-out',str(out/'prototype_db')],stdout=log,stderr=log,check=True)
        subprocess.run(['blastp','-query',str(out/'source_queries.faa'),'-db',str(out/'prototype_db'),'-evalue','1e-5','-max_target_seqs','10000','-max_hsps','1','-num_threads','4',
            '-outfmt','6 qseqid sseqid pident qlen slen qstart qend sstart send bitscore evalue','-out',str(out/'source_vs_prototypes.tsv')],stderr=log,check=True)
    overlaps=collections.defaultdict(list);maximum={}
    for line in (out/'source_vs_prototypes.tsv').read_text().splitlines():
        q,s,identity,qlen,slen,qs,qe,ss,se,score,evalue=line.split('\t');identity=float(identity);qcov=(int(qe)-int(qs)+1)/int(qlen);scov=(int(se)-int(ss)+1)/int(slen)
        shorter=qcov if int(qlen)<=int(slen) else scov
        for uid in links[q]:
            if min(qcov,scov)>=.7:maximum[uid]=max(maximum.get(uid,0),identity)
            if identity>=ip['homolog_identity_percent'] and shorter>=ip['homolog_shorter_sequence_coverage']:
                overlaps[refined[uid]].append(dict(query_id=q,prototype_id=s,identity=identity,shorter_coverage=shorter))
    prior=set(previous['independent_unit_ids']);rows=[];accepted=[]
    for r in candidates:
        uid=r['unit_id'];reasons=[]
        if uid not in prior:reasons.append('prior_core_audit_excluded')
        if set(r['genus_taxids'])&training_genera:reasons.append('prototype_genus_overlap')
        if refined[uid] in overlaps:reasons.append('prototype_full_or_author_core_homology_overlap')
        if not reasons:accepted.append(uid)
        rows.append(dict(unit_id=uid,split=r['split'],label=r['label'],family_categories=';'.join(r['categories']),refined_component=refined[uid],
            eligible=not reasons,reasons=';'.join(reasons),genus_taxids=';'.join(r['genus_taxids']),max_comparable_prototype_identity=maximum.get(uid),overlap_evidence=json.dumps(overlaps.get(refined[uid],[]))))
    write_tsv(out/'per_unit_audit.tsv',rows)
    counts={s:dict(total=sum(r['split']==s for r in rows),independent=sum(r['split']==s and r['eligible'] for r in rows),
        positive=sum(r['split']==s and r['eligible'] and r['label']=='positive' for r in rows),negative=sum(r['split']==s and r['eligible'] and r['label']=='negative' for r in rows)) for s in ['calibration','final_validation']}
    summary=dict(status='source_audit_frozen_before_confirmation_prediction',**before,independent_unit_ids=sorted(accepted),partition_counts=counts,
        limitations='Own-library homology/genus audit under the registered operational thresholds. Public Pfam training membership is not independently reconstructable; no claim of no evolutionary relatedness or per-protein experimental activity.',
        artifacts={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.suffix in ['.faa','.json','.tsv']})
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(counts,indent=2))


if __name__=='__main__':main()
