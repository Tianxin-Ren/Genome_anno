#!/usr/bin/env python3
"""Add local TMbed topology evidence without changing catalytic/complex/function labels."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT,write_tsv


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def runtime_versions(python):
    code="import importlib.metadata as m,json,platform; print(json.dumps(dict(python=platform.python_version(),distributions={d.metadata['Name']:d.version for d in m.distributions()}),sort_keys=True))"
    return json.loads(subprocess.check_output([str(python),'-c',code],text=True))


def verify_runtime(policy,model,root=ROOT):
    if runtime_versions(Path(root)/policy['python'])!=model['runtime_versions']:raise ValueError('Localization Python/package versions changed; revalidate the external model')


def table(path):
    with Path(path).open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def verify_source_file(run,manifest,name):
    expected=manifest.get('normalized_proteins_sha256') if name=='proteins.faa' else manifest['output_sha256'].get(name)
    if expected!=sha(Path(run)/name):raise ValueError('Source annotation changed: '+str(Path(run)/name))


def parse_topology(path,expected):
    lines=Path(path).read_text().splitlines()
    if len(lines)%3:raise ValueError('TMbed output must contain complete three-line records')
    result={}
    for i in range(0,len(lines),3):
        header,sequence,labels=lines[i:i+3];qid=header[1:].split()[0] if header.startswith('>') else ''
        if not qid or qid in result or qid not in expected:raise ValueError('Unknown/repeated TMbed query ownership')
        if sequence!=expected[qid] or len(labels)!=len(sequence) or set(labels)-set('BbHhSio'):
            raise ValueError('TMbed sequence, coordinates or output alphabet mismatch')
        result[qid]=labels
    if set(result)!=set(expected):raise ValueError('TMbed omitted selected queries')
    return result


def segments(labels):
    result=[];start=0
    for stop in range(1,len(labels)+1):
        if stop<len(labels) and labels[stop]==labels[start]:continue
        code=labels[start]
        if code in 'BbHhS':
            result.append(dict(start=start+1,end=stop,kind='alpha_helix' if code in 'Hh' else 'beta_strand' if code in 'Bb' else 'signal_peptide',
                orientation='predicted_inside_to_outside' if code in 'BH' else 'predicted_outside_to_inside' if code in 'bh' else 'not_assigned'))
        start=stop
    return result


def summarize(labels,partial):
    spans=segments(labels);tm=[s for s in spans if s['kind']!='signal_peptide'];signal=[s for s in spans if s['kind']=='signal_peptide']
    return dict(localization_status='topology_predicted_compartment_unresolved',transmembrane_segments=len(tm),
        signal_peptide_predicted=bool(signal),N_terminal_absence_interpretation='not_interpretable_partial_or_unknown_CDS' if partial else 'prediction_only_start_not_experimentally_verified',
        signal_pathway='unresolved_Tat_vs_Sec_or_other',cell_compartment='unresolved',
        membrane_binding_without_TM='not_assessed',confidence='provisional_external_topology_model_not_hydrogenase_calibrated',
        limitations='No TM/SP prediction does not prove cytoplasmic location; co-export of a catalytic partner, peripheral binding, lipid anchors and actual translocation are not established.')


def project_node(labels,positions):
    positions=set(positions)
    if any(not isinstance(p,int) or isinstance(p,bool) or not 1<=p<=len(labels) for p in positions):raise ValueError('Module topology coordinates outside query')
    return dict(mapped_residues=len(positions),mapped_TM_residues=sum(labels[p-1] in 'BbHh' for p in positions),
        mapped_SP_residues=sum(labels[p-1]=='S' for p in positions))


def verify_completed(directory,root=ROOT):
    directory=Path(directory);root=Path(root);report=json.loads((directory/'localization_manifest.json').read_text())
    if report['status']!='complete' or not report['inputs_and_model_unchanged']:raise ValueError('Incomplete localization execution')
    before=report['before'];policy_path=root/'config/localization.json';policy=json.loads(policy_path.read_text());model_path=root/policy['model_manifest']
    if sha(policy_path)!=before['policy_sha256'] or sha(model_path)!=before['model_manifest_sha256'] or sha(root/'scripts/localization.py')!=before['script_sha256']:
        raise ValueError('Localization predictor changed after execution')
    model=json.loads(model_path.read_text());verify_runtime(policy,model,root)
    paths={**{Path(n):d for n,d in before['inputs'].items()},**{directory/n:d for n,d in report['artifacts'].items()},
        **{root/n:d for n,d in model['files'].items()}}
    if any(sha(p)!=digest for p,digest in paths.items()):raise ValueError('Localization input, prediction or model changed after execution')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs',nargs='+',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    runs=sorted({p.resolve() for p in args.runs})
    if len(runs)!=len(args.runs):raise ValueError('Duplicate source run')
    policy_path=ROOT/'config/localization.json';policy=json.loads(policy_path.read_text())
    model_path=ROOT/policy['model_manifest'];model=json.loads(model_path.read_text())
    if model['status']!='downloaded_and_verified':raise ValueError('Localization model not ready')
    verify_runtime(policy,model)
    for name,digest in model['files'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Localization model/code changed: '+name)
    out=args.out.resolve()
    if out.exists():raise FileExistsError('Use a fresh localization output directory')
    owners=[];queries={};nodes=[];inputs={};before={}
    for run in runs:
        manifest_path=run/'run_manifest.json';manifest=json.loads(manifest_path.read_text())
        if manifest['status']!='complete':raise ValueError('Source hydrogenase run is incomplete')
        files=['run_manifest.json','proteins.faa','genes.tsv','catalytic_anchors.tsv','module_graph_nodes.tsv']
        for name in files:
            path=run/name;inputs[str(path)]=sha(path)
            if name!='run_manifest.json':verify_source_file(run,manifest,name)
        proteins=list(SeqIO.parse(run/'proteins.faa','fasta'));byid={p.id:str(p.seq).upper() for p in proteins}
        if len(byid)!=len(proteins):raise ValueError('Repeated protein ID in source run')
        genes={g['id']:g for g in table(run/'genes.tsv')};selected={a['protein_id'] for a in table(run/'catalytic_anchors.tsv')}
        for node in table(run/'module_graph_nodes.tsv'):
            choices=json.loads(node['candidate_realizations']);resources={}
            for c in choices:
                for pid,positions in c['query_resources'].items():resources.setdefault(pid,set()).update(positions)
            for pid,positions in json.loads(node['query_resources']).items():resources.setdefault(pid,set()).update(positions)
            selected.update(resources)
            for pid,positions in resources.items():nodes.append(dict(source_run=str(run),protein_id=pid,anchor=node['anchor'],node_id=node['node_id'],
                role=node['role'],module_category=node['module_category'],selected_for_node=pid in json.loads(node['query_resources']),positions=positions))
        if selected-set(byid):raise ValueError('Module/anchor protein missing from source FASTA')
        for pid in sorted(selected):
            sequence=byid[pid];qid=hashlib.sha256(sequence.encode()).hexdigest();queries[qid]=sequence
            owners.append(dict(source_run=str(run),protein_id=pid,query_id=qid,partial=genes.get(pid,{}).get('partial')!='False'))
    out.mkdir(parents=True,exist_ok=False)
    fasta=out/'queries.faa';SeqIO.write([SeqRecord(Seq(s),id=k,description='local_topology_query') for k,s in sorted(queries.items())],fasta,'fasta')
    before=dict(policy_sha256=sha(policy_path),model_manifest_sha256=sha(model_path),script_sha256=sha(Path(__file__)),inputs=inputs)
    command=[str(ROOT/policy['python']),'-m','tmbed','predict','-f',str(fasta),'-p',str(out/'raw.pred'),'--out-format','4',
        '--no-use-gpu','--threads',str(policy['threads']),'--batch-size',str(policy['batch_residues']),'--model-dir',str(ROOT/policy['encoder_directory'])]
    report=dict(status='running',scope='Auxiliary topology evidence only; upstream identity/structure/function and bin ownership are unchanged.',
        before=before,command=command,selected_protein_owners=len(owners),unique_sequences=len(queries),sources=len(runs),model_source=policy['source'])
    record=out/'localization_manifest.json';record.write_text(json.dumps(report,indent=2)+'\n');start=time.monotonic()
    if queries:
        env=dict(os.environ);env.pop('HF_HOME',None);env.update(HF_HUB_OFFLINE='1',TRANSFORMERS_OFFLINE='1',HF_HUB_DISABLE_TELEMETRY='1',HF_HUB_DISABLE_XET='1')
        with (out/'command.log').open('w') as log:process=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env)
        report['exit_code']=process.returncode
        if process.returncode:
            report.update(status='prediction_failed',elapsed_seconds=time.monotonic()-start);record.write_text(json.dumps(report,indent=2)+'\n')
            raise RuntimeError('Local TMbed prediction failed; see command.log')
        try:predictions=parse_topology(out/'raw.pred',queries)
        except ValueError as error:
            report.update(status='invalid_prediction_output',error=str(error),elapsed_seconds=time.monotonic()-start)
            record.write_text(json.dumps(report,indent=2)+'\n');raise
    else:predictions={}
    rows=[];spans=[];module_rows=[]
    by_owner={(o['source_run'],o['protein_id']):o for o in owners}
    for owner in owners:
        labels=predictions[owner['query_id']];rows.append(dict(**owner,**summarize(labels,owner['partial'])))
        spans.extend(dict(source_run=owner['source_run'],protein_id=owner['protein_id'],query_id=owner['query_id'],**s) for s in segments(labels))
    for node in nodes:
        owner=by_owner[node['source_run'],node['protein_id']];positions=node.pop('positions')
        module_rows.append(dict(**node,query_id=owner['query_id'],**project_node(predictions[owner['query_id']],positions),
            evidence_scope='Predicted residue labels projected only onto this module candidate; no physical interaction or function promotion.'))
    write_tsv(out/'localization.tsv',rows,['source_run','protein_id','query_id','partial',*summarize('',True)])
    write_tsv(out/'localization_segments.tsv',spans,['source_run','protein_id','query_id','start','end','kind','orientation'])
    write_tsv(out/'module_localization.tsv',module_rows,['source_run','protein_id','anchor','node_id','role','module_category','selected_for_node','query_id','mapped_residues','mapped_TM_residues','mapped_SP_residues','evidence_scope'])
    unchanged=all(sha(Path(p))==digest for p,digest in inputs.items()) and sha(policy_path)==before['policy_sha256'] and sha(model_path)==before['model_manifest_sha256'] and sha(Path(__file__))==before['script_sha256']
    unchanged &= runtime_versions(ROOT/policy['python'])==model['runtime_versions']
    unchanged &= all(sha(ROOT/name)==digest for name,digest in model['files'].items())
    report.update(status='complete' if unchanged else 'source_or_model_changed',elapsed_seconds=time.monotonic()-start,inputs_and_model_unchanged=unchanged,
        artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file() and p!=record})
    record.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['status','sources','selected_protein_owners','unique_sequences','elapsed_seconds']},indent=2))
    if not unchanged:raise SystemExit(1)


if __name__=='__main__':main()
