#!/usr/bin/env python3
"""Refine identity independence with actual training spans and author domain fragments."""
import argparse
import collections
import csv
import hashlib
import json
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT,write_tsv
from fetch_validation_metadata import accession


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def component_eligibility(splits,training_exposed):
    if training_exposed or 'reference_overlap' in splits:return 'training_or_prior_exposure'
    if 'source_review' in splits:return 'source_review'
    if 'split_bridge' in splits or len(splits)>1:return 'core_split_bridge'
    return 'independent'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--corpus',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();corpus=args.corpus.resolve();out=args.out.resolve()
    if (out/'summary.json').exists():parser.error('Use a fresh source audit directory')
    out.mkdir(parents=True,exist_ok=True)
    frozen=json.loads((corpus/'freeze.json').read_text())
    for name in ['units.json','sequence_variant_audit.tsv']:
        if sha(corpus/name)!=frozen['artifacts'][name]:raise ValueError('Frozen identity data changed: '+name)
    units=json.loads((corpus/'units.json').read_text());byid={r['unit_id']:r for r in units}
    policy=json.loads((ROOT/'config/identity_corpus_policy.json').read_text())
    profiles_path=ROOT/'references/processed/profile_manifest.json';profiles=json.loads(profiles_path.read_text())
    training_ids={uid for p in profiles for uid in p.get('training_corpus_unit_ids',[])}
    assert all(byid[uid]['split']=='development' for uid in training_ids)
    # Source-only policy is recorded before any new alignment or partition prediction.
    source_policy=dict(original_freeze_sha256=sha(corpus/'freeze.json'),profile_manifest_sha256=sha(profiles_path),
        definition='Same preregistered >=70% identity and >=70% shorter-sequence coverage; author catalytic fragments also link biological units.',
        rules=['Keep all original units and predictions; mark extra overlap before calibration/final prediction.',
            'Union original full-protein components with author FeFe fragment components.',
            'Propagate prior exposure, training overlap, source review and cross-split bridges to the whole refined component.',
            'No reassignment to another accepted partition and no selection based on prediction success.',
            'Report all frozen-partition metrics separately from the source-audited independent subset.'])
    (out/'policy_before_alignment.json').write_text(json.dumps(source_policy,indent=2)+'\n')
    training={}
    def add_training(sequence,source):
        digest=hashlib.sha256(sequence.encode()).hexdigest()
        row=training.setdefault(digest,dict(training_id='training_'+digest[:20],sequence=sequence,sources=[]))
        row['sources'].append(source)
    for name in ['references.faa','reference_source_proteins.faa']:
        with (ROOT/'references/processed'/name).open() as handle:
            for rec in SeqIO.parse(handle,'fasta'):add_training(str(rec.seq),dict(kind='reference_full',source_file=name,reference_id=rec.id))
    for profile in profiles:
        if profile.get('retained_alignment_columns_1based'):
            path=ROOT/'references/processed'/(profile['profile']+'.aligned.faa')
            keep=[i-1 for i in profile['retained_alignment_columns_1based']]
            with path.open() as handle:
                for rec in SeqIO.parse(handle,'fasta'):
                    seq=str(rec.seq);columns=[i for i in keep if seq[i]!='-']
                    left,right=min(columns),max(columns)
                    start=sum(a!='-' for a in seq[:left])+1;end=sum(a!='-' for a in seq[:right+1])
                    span=seq[left:right+1].replace('-','')
                    assert len(span)==end-start+1
                    add_training(span,dict(kind='actual_training_core_contiguous_span',profile=profile['profile'],reference_id=rec.id,
                        protein_start=start,protein_end=end,alignment_sha256=sha(path)))
        for uid in profile.get('training_corpus_unit_ids',[]):
            add_training(byid[uid]['sequence'],dict(kind='development_training_full',unit_id=uid))
    SeqIO.write([SeqRecord(Seq(r['sequence'].replace('U','C')),id=r['training_id'],description='') for r in training.values()],out/'training_spans.faa','fasta')
    (out/'training_spans.json').write_text(json.dumps(list(training.values()),indent=2)+'\n')
    SeqIO.write([SeqRecord(Seq(r['sequence'].replace('U','C')),id=r['unit_id'],description='') for r in units],out/'source_units.faa','fasta')
    aliases=collections.defaultdict(set)
    for unit in units:
        for source in unit['sources']:
            if source['source'].startswith('HydDB'):
                aliases[source['accession']].add(unit['unit_id'])
    raw_core=ROOT/'references/raw/validation_corpus/HydDB_all_hydrogenases.faa'
    if sha(raw_core)!=frozen['sources'][str(raw_core.relative_to(ROOT))]:raise ValueError('Author domain fragment source changed')
    fragments={}
    with raw_core.open() as handle:
        for rec in SeqIO.parse(handle,'fasta'):
            if '[FeFe]' not in rec.description:continue
            ids=aliases[accession(rec.id)]
            if not ids:raise ValueError('Author FeFe fragment has no frozen biological-unit link')
            seq=str(rec.seq).upper();digest=hashlib.sha256(seq.encode()).hexdigest()
            entry=fragments.setdefault(digest,dict(fragment_id='fragment_'+digest[:20],sequence=seq,unit_ids=set(),source_headers=[]))
            entry['unit_ids'].update(ids);entry['source_headers'].append(rec.description)
    for row in fragments.values():row['unit_ids']=sorted(row['unit_ids'])
    SeqIO.write([SeqRecord(Seq(r['sequence'].replace('U','C')),id=r['fragment_id'],description='') for r in fragments.values()],out/'author_fragments.faa','fasta')
    (out/'author_fragments.json').write_text(json.dumps(list(fragments.values()),indent=2)+'\n')
    commands=[]
    def run(cmd):
        commands.append(cmd);(out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
        with (out/'blast.log').open('a') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    for name in ['training_spans','author_fragments']:
        run(['makeblastdb','-in',str(out/(name+'.faa')),'-dbtype','prot','-out',str(out/name)])
    fmt='6 qseqid sseqid pident qlen slen qstart qend sstart send length bitscore evalue'
    for name,query,db,count in [('units_vs_training','source_units','training_spans',len(training)),
        ('fragments_vs_training','author_fragments','training_spans',len(training)),
        ('fragment_pairs','author_fragments','author_fragments',len(fragments))]:
        print(name,flush=True)
        run(['blastp','-query',str(out/(query+'.faa')),'-db',str(out/db),'-evalue','1e-5','-max_hsps','1',
             '-max_target_seqs',str(count),'-num_threads','4','-outfmt',fmt,'-out',str(out/(name+'.tsv'))])
    parent={uid:uid for uid in byid}
    def find(a):
        while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
        return a
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    for row in units:union(row['unit_id'],row['homolog_component'])
    fragment_byid={r['fragment_id']:r for r in fragments.values()}
    for fragment in fragments.values():
        for uid in fragment['unit_ids']:union(uid,fragment['unit_ids'][0])
    def alignments(name):
        with (out/(name+'.tsv')).open() as handle:
            for line in handle:
                q,s,identity,qlen,slen,qs,qe,ss,se,*_=line.split('\t')
                qcov=(int(qe)-int(qs)+1)/int(qlen);scov=(int(se)-int(ss)+1)/int(slen)
                yield q,s,float(identity),qcov,scov,qcov if int(qlen)<=int(slen) else scov
    for q,s,identity,qcov,scov,shorter in alignments('fragment_pairs'):
        if identity>=policy['homolog_identity_percent'] and shorter>=policy['homolog_shorter_sequence_coverage']:
            union(fragment_byid[q]['unit_ids'][0],fragment_byid[s]['unit_ids'][0])
    exposed=set(training_ids);maximum={}
    for name in ['units_vs_training','fragments_vs_training']:
        for q,s,identity,qcov,scov,shorter in alignments(name):
            ids=[q] if name=='units_vs_training' else fragment_byid[q]['unit_ids']
            if identity>=policy['homolog_identity_percent'] and shorter>=policy['homolog_shorter_sequence_coverage']:exposed.update(ids)
            if name=='units_vs_training' and min(qcov,scov)>=policy['remote_group4_both_sequence_coverage']:
                maximum[q]=max(maximum.get(q,0),identity)
    groups=collections.defaultdict(list)
    for row in units:groups[find(row['unit_id'])].append(row)
    audit=[];independent=[]
    for cid,members in sorted(groups.items()):
        splits={r['split'] for r in members};status=component_eligibility(splits,bool({r['unit_id'] for r in members}&exposed))
        for row in members:
            uid=row['unit_id'];eligible=status=='independent' and row['split'] in ['development','calibration','final_validation']
            if eligible:independent.append(uid)
            remote=bool(eligible and any(c.startswith('[NiFe] Group 4') for c in row['categories']) and uid in maximum and maximum[uid]<=policy['remote_group4_max_training_identity_percent'])
            audit.append(dict(unit_id=uid,original_split=row['split'],refined_component=cid,eligibility=status,independent_for_identity=eligible,
                directly_matches_current_training=uid in exposed,component_original_splits=';'.join(sorted(splits)),max_comparable_current_training_identity=maximum.get(uid),
                remote_Group4_after_core_audit=remote,genus_names=';'.join(row['genus_names']),categories=';'.join(row['categories'])))
    # Verify the accepted source-only refinements share neither components nor genera.
    accepted={s:[r for r in audit if r['original_split']==s and r['independent_for_identity']] for s in ['development','calibration','final_validation']}
    for a in accepted:
        for b in accepted:
            if a!=b:
                assert not {r['refined_component'] for r in accepted[a]} & {r['refined_component'] for r in accepted[b]}
                assert not {g for r in accepted[a] for g in byid[r['unit_id']]['genus_taxids']} & {g for r in accepted[b] for g in byid[r['unit_id']]['genus_taxids']}
    write_tsv(out/'per_unit_audit.tsv',audit)
    counts={s:dict(independent_units=len(rs),positive=sum(byid[r['unit_id']]['label']=='positive' for r in rs),
        negative=sum(byid[r['unit_id']]['label']=='negative' for r in rs),remote_Group4=sum(r['remote_Group4_after_core_audit'] for r in rs),
        components=len({r['refined_component'] for r in rs})) for s,rs in accepted.items()}
    summary=dict(status='source_audit_frozen_before_partition_prediction',corpus=str(corpus),corpus_freeze_sha256=sha(corpus/'freeze.json'),
        profile_manifest_sha256=sha(profiles_path),training_unit_ids=sorted(training_ids),independent_unit_ids=sorted(independent),partition_counts=counts,
        source_only=True,no_HMM_prediction_used=True,author_fragment_sequences=len(fragments),training_span_sequences=len(training),
        author_fragment_source_sha256=sha(raw_core),artifacts={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.suffix in ['.faa','.json','.tsv']},
        limitations=['An operational source-only domain/full-protein independence audit, not functional truth.',
            'Original partitions and all records retained; this audit can reduce eligible final sample size without any prediction-based filtering.',
            'Nonexact author fragments link only for conservative leakage exclusion, never as replacement sequence or activity truth.'])
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(counts,indent=2),flush=True)


if __name__=='__main__':main()
