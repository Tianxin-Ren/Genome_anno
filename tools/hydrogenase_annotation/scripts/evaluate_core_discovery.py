#!/usr/bin/env python3
"""Compare existing-seed FeFe discovery masks on development data only."""
import argparse
import collections
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT, write_tsv
from profile_models import build_profile
from evaluate_profile_design import search
from evaluate_identity_corpus import identity_metrics
from discovery_evidence import load_domains,select_domains


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--corpus',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();out=args.out.resolve();corpus=args.corpus.resolve()
    if (out/'summary.json').exists():parser.error('Use a fresh diagnostic directory')
    out.mkdir(parents=True,exist_ok=True)
    frozen=json.loads((corpus/'freeze.json').read_text())
    for filename in ['units.json','development.faa']:
        if hashlib.sha256((corpus/filename).read_bytes()).hexdigest()!=frozen['artifacts'][filename]:
            raise ValueError('Frozen development data changed')
    units=json.loads((corpus/'units.json').read_text())
    truth=[r for r in units if r['split']=='development' and (r['label']=='negative' or any(c.startswith('[FeFe]') for c in r['categories']))]
    query=corpus/'development.faa'
    ids={s.id for s in SeqIO.parse(ROOT/'references/processed/FeFe.faa','fasta')}
    records=[s for s in SeqIO.parse(ROOT/'references/processed/references.faa','fasta') if s.id in ids]
    anchor=next(s for s in records if s.id=='FeFe_CpI_monomer__CpI')
    core=dict(reference_id=anchor.id,sequence_sha256=hashlib.sha256(str(anchor.seq).encode()).hexdigest(),start=285,end=518,
        required_positions={'300':'C','355':'C','499':'C','503':'C'},source='https://www.rcsb.org/structure/1FEH',
        source_site_id='1FEH__C_581_SF4',source_sha256=hashlib.sha256((ROOT/'references/raw/structures/1FEH.cif').read_bytes()).hexdigest(),
        scope='H-cluster ligand span plus 15-residue flanks; discovery core, not complete H-domain or activity proof')
    (out/'design_before_prediction.json').write_text(json.dumps(dict(training_ids=sorted(ids),core=core,
        alternatives=['existing_legacy','occupancy_only_0.6','structure_anchored_ligand_span_occupancy_0.6'],
        unchanged_thresholds=dict(domain_evalue=1e-5,model_coverage=.55),data_split='development_only'),indent=2)+'\n')
    models=[('legacy',ROOT/'references/processed/FeFe.hmm',None)]
    for name,region in [('FeFe_occupancy',None),('FeFe_Hcluster_core',core)]:
        model,manifest=build_profile(records,out/name,name,0.6,region)
        models.append((name,model,manifest))
    rows=[];metrics={};manifests={}
    policy=json.loads((ROOT/'config/discovery_evidence.json').read_text())
    for name,model,manifest in models:
        domains=search(model,query,out,name,reference_sequences=policy['evalue_reference_sequences'])
        predicted={d['protein_id'] for d in domains if d['passes']}
        metrics[name]=identity_metrics(truth,predicted)
        metrics[name]['missed_positive_categories']=dict(collections.Counter(c for r in truth if r['label']=='positive' and r['unit_id'] not in predicted for c in r['categories']))
        selected=select_domains(load_domains(out/(name+'.domtbl'),{'FeFe' if name=='legacy' else name:'FeFe'}),
            dict(hmm_domain_evalue=1e-5,hmm_coverage=.55),policy)
        chained={r['protein_id'] for r in selected}
        metrics[name]['with_collinear_domains']=identity_metrics(truth,chained)
        metrics[name]['with_collinear_domains']['missed_positive_categories']=dict(collections.Counter(c for r in truth if r['label']=='positive' and r['unit_id'] not in chained for c in r['categories']))
        write_tsv(out/(name+'_selected.tsv'),selected)
        rows.extend(dict(experiment=name,**d) for d in domains)
        manifests[name]=manifest
    write_tsv(out/'domains.tsv',rows)
    report=dict(status='development_design_comparison',metrics=metrics,model_manifests=manifests,
        corpus_freeze_sha256=hashlib.sha256((corpus/'freeze.json').read_bytes()).hexdigest(),
        limitation='Development-only model selection. No external corpus sequence is added to training; negatives are candidate-stage hits, not final function calls.')
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(metrics,indent=2))


if __name__=='__main__':main()
