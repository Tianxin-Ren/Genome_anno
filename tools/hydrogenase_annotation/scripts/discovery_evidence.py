"""Retain individual HMM domains and bounded collinear candidate evidence."""
import collections
import json
import math
import hashlib
import subprocess
from pathlib import Path
from tool_paths import resolve_tool


FIELDS='protein_id family profile independent_evalue domain_bitscore profile_coverage query_start query_end model_start model_end evidence_kind evalue_scope bitscore_scope model_segments query_segments component_domains evalue_reference_sequences'.split()


def _hmmsearch_version_sha(root):
    out = subprocess.check_output([str(resolve_tool(root, 'hmmsearch')), '-h'], text=True, stderr=subprocess.STDOUT)
    line = next((l.strip() for l in out.splitlines() if l.strip().startswith('# HMMER ')), '')
    return hashlib.sha256(line.encode('utf-8')).hexdigest()


def stage_fingerprint(root):
    files=['scripts/discovery_evidence.py','config/discovery_profiles.json','config/discovery_evidence.json',
           'references/processed/catalytic_families.hmm','references/processed/profile_manifest.json']
    thresholds=json.loads((root/'config/thresholds.json').read_text())
    return dict(files={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in files},
        thresholds={k:thresholds[k] for k in ['hmm_domain_evalue','hmm_coverage']},
        hmmsearch_binary_sha256=_hmmsearch_version_sha(root))


def search_discovery(fasta,out,root,threads,thresholds):
    dom=out/'family.domtbl'
    tbl=out/'family.tbl'
    policy_path=root/'config/discovery_evidence.json';policy=json.loads(policy_path.read_text())
    size=policy['evalue_reference_sequences']
    if not isinstance(size,int) or isinstance(size,bool) or size<=0:raise ValueError('Invalid discovery E-value reference size')
    with (out/'hmmsearch.log').open('w') as log:
        subprocess.run([str(resolve_tool(root, 'hmmsearch')),'--cpu',str(threads),'-Z',str(size),'--domtblout',str(dom),'--tblout',str(tbl),str(root/'references/processed/catalytic_families.hmm'),str(fasta)],stdout=log,check=True)
    manifest_path=root/'references/processed/profile_manifest.json'
    families={r['profile']:r['family'] for r in json.loads(manifest_path.read_text())}
    hits=select_domains(load_domains(dom,families),thresholds,policy)
    return hits,manifest_path,policy_path


def load_targets(path, profile_families):
    """Read HMMER full-sequence hits, independently of domain thresholds."""
    path = Path(path)
    targets = []
    with path.open() as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split()
            # Standard hmmsearch tblout has four name columns plus fourteen
            # full/domain statistic columns; a missing tail is not a valid row.
            if len(fields) < 18:
                raise ValueError('Invalid HMM full-sequence row')
            profile = fields[2]
            if profile not in profile_families:
                raise ValueError('Invalid HMM full-sequence row or unknown profile')
            try:
                evalue, score, bias = (float(fields[i]) for i in (4, 5, 6))
            except (TypeError, ValueError):
                raise ValueError('Invalid HMM full-sequence score')
            if (not math.isfinite(evalue) or evalue < 0 or
                    not math.isfinite(score) or not math.isfinite(bias) or bias < 0):
                raise ValueError('Invalid HMM full-sequence score')
            targets.append(dict(protein_id=fields[0], family=profile_families[profile],
                profile=profile, sequence_evalue=evalue, sequence_bitscore=score,
                sequence_bias=bias))
    return targets


def load_domains(path,profile_families):
    domains=[]
    with path.open() as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():continue
            c=line.split()
            if len(c)<22 or c[3] not in profile_families:
                raise ValueError('Invalid HMM domain row or unknown profile')
            row=dict(protein_id=c[0],family=profile_families[c[3]],profile=c[3],query_length=int(c[2]),model_length=int(c[5]),
                independent_evalue=float(c[12]),domain_bitscore=float(c[13]),model_start=int(c[15]),model_end=int(c[16]),query_start=int(c[17]),query_end=int(c[18]))
            if (not 1<=row['model_start']<=row['model_end']<=row['model_length'] or
                not 1<=row['query_start']<=row['query_end']<=row['query_length'] or
                not math.isfinite(row['independent_evalue']) or row['independent_evalue']<0 or not math.isfinite(row['domain_bitscore'])):
                raise ValueError('Invalid HMM coordinates or scores')
            domains.append(row)
    return domains


def chain_compatible(a,b,policy):
    if (a['protein_id'],a['profile'],a['model_length'],a['query_length'])!=(b['protein_id'],b['profile'],b['model_length'],b['query_length']):
        return False
    if b['model_start']<=a['model_end'] or b['query_start']<=a['query_end']:
        return False
    query_gap=b['query_start']-a['query_end']-1
    model_gap=b['model_start']-a['model_end']-1
    return query_gap-model_gap <= a['model_length']*policy['max_extra_query_gap_model_fraction']


def best_chain(domains,policy):
    ordered=sorted(domains,key=lambda d:(d['query_start'],d['model_start'],d['query_end'],d['model_end']))
    chains=[]
    def score(path):
        return sum(d['model_end']-d['model_start']+1 for d in path),sum(d['domain_bitscore'] for d in path),-len(path)
    for i,d in enumerate(ordered):
        candidates=[[d]]+[chains[j]+[d] for j in range(i) if chain_compatible(ordered[j],d,policy)]
        chains.append(max(candidates,key=score))
    return max(chains,key=score) if chains else []


def evidence_row(components,policy):
    first,last=components[0],components[-1]
    chained=len(components)>1
    return dict(protein_id=first['protein_id'],family=first['family'],profile=first['profile'],
        independent_evalue=max(d['independent_evalue'] for d in components),domain_bitscore=sum(d['domain_bitscore'] for d in components),
        profile_coverage=sum(d['model_end']-d['model_start']+1 for d in components)/first['model_length'],
        query_start=first['query_start'],query_end=last['query_end'],model_start=first['model_start'],model_end=last['model_end'],
        evidence_kind='collinear_domain_chain' if chained else 'single_domain',
        evalue_scope='max_component_domain_E_not_combined_pvalue' if chained else 'single_domain_E',
        bitscore_scope='sum_component_bitscores' if chained else 'single_domain_bitscore',
        model_segments=json.dumps([[d['model_start'],d['model_end']] for d in components]),
        query_segments=json.dumps([[d['query_start'],d['query_end']] for d in components]),component_domains=json.dumps(components),
        evalue_reference_sequences=policy['evalue_reference_sequences'])


def select_domains(domains,thresholds,policy):
    groups=collections.defaultdict(list)
    for d in domains:
        if d['independent_evalue']<=thresholds['hmm_domain_evalue']:
            groups[(d['protein_id'],d['profile'])].append(d)
    selected=[]
    for key,values in sorted(groups.items()):
        single=[evidence_row([d],policy) for d in values if (d['model_end']-d['model_start']+1)/d['model_length']>=thresholds['hmm_coverage']]
        if single:
            selected.extend(single)
            continue
        chain=best_chain(values,policy)
        if len(chain)>1:
            row=evidence_row(chain,policy)
            if row['profile_coverage']>=thresholds['hmm_coverage']:
                selected.append(row)
    return selected
