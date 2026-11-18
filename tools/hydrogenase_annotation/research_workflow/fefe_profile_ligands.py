"""Uncalibrated FeFe ligand-pattern mapping through the official PF02906 model."""
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = Path(__file__).with_name('fefe_profile_ligands_model.json')
RULE_VERSION = 'FeFe-PF02906-profile-ligand-pattern-v1'
AMBIGUOUS = set('XBZJOU')


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _load_model(root):
    model = json.loads(MODEL_PATH.read_text())
    if model['rule_version'] != RULE_VERSION:
        raise ValueError('FeFe profile mapping rule changed')
    hmm, refs = Path(root)/model['hmm_path'], Path(root)/model['reference_model_path']
    ref = next(r for r in json.loads(refs.read_text()) if r['reference_id'] == model['reference_id'])
    if model['reference_positions'] != [300,355,499,503] or model['expected_match_columns'] != [77,130,272,276]:
        raise ValueError('FeFe reference ligand coordinates changed')
    if any(ref['ligands'].get(str(p)) != 'C' or ref['sequence'][p-1] != 'C' for p in model['reference_positions']):
        raise ValueError('FeFe reference ligand residues changed')
    return model, hmm, ref, refs


def _fasta(path):
    with Path(path).open() as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
    result = {r.id:str(r.seq).upper() for r in records}
    if len(records) != len(result):
        raise ValueError('Duplicate FASTA identifiers')
    if any(not s or any(not c.isalpha() for c in s) for s in result.values()):
        raise ValueError('Invalid FeFe profile input sequence')
    return result


def _stockholm(text):
    if not text.startswith('# STOCKHOLM 1.0') or text.rstrip().splitlines()[-1] != '//':
        raise ValueError('Incomplete Stockholm alignment')
    seqs, pp, rf = {}, {}, ''
    for line in text.splitlines():
        if line.startswith('#=GC RF'):
            rf += line.split(None,2)[2]
        elif line.startswith('#=GR '):
            x = line.split(None,3)
            if len(x) == 4 and x[2] == 'PP': pp[x[1]] = pp.get(x[1],'') + x[3]
        elif line and not line.startswith(('#','//')):
            x = line.split()
            if len(x) != 2: raise ValueError('Invalid Stockholm sequence row')
            seqs[x[0]] = seqs.get(x[0],'') + x[1]
    return seqs, pp, rf


def _query_domain(feature, model, min_domain_coverage=None, protein_id=None, query_length=None):
    minimum = model['min_domain_coverage'] if min_domain_coverage is None else float(min_domain_coverage)
    accepted = []
    for d in feature.get('domain_evidence',[]):
        if d.get('accession') != 'PF02906' or d.get('supported') is not True or d.get('evidence_kind') != 'single_domain': continue
        if len(d.get('components',[])) != 1: continue
        c = d['components'][0]
        if c.get('protein_id') != protein_id or c.get('model_sha256') != model['hmm_sha256']: continue
        if str(c.get('accession','')).split('.')[0] != 'PF02906': continue
        keys = ['sequence_bitscore','domain_bitscore','sequence_ga','domain_ga','query_length','query_start','query_end','model_length','model_start','model_end','model_coverage']
        if any(k not in c for k in keys) or any(k not in d for k in ['model_coverage','query_core_overlap']): continue
        values = [float(c[k]) for k in keys] + [float(d['model_coverage']),float(d['query_core_overlap']),minimum]
        if not all(math.isfinite(x) for x in values): raise ValueError('Non-finite PF02906 evidence')
        if any(float(c[k]) != int(c[k]) for k in ['query_length','query_start','query_end','model_length','model_start','model_end']): continue
        if c['query_length'] != query_length or int(c['model_length']) != 280: continue
        lo, hi = int(c['query_start']), int(c['query_end'])
        if not 1 <= lo <= hi <= query_length: continue
        if not 1 <= int(c['model_start']) <= int(c['model_end']) <= 280: continue
        if float(c['sequence_ga']) != model['official_sequence_ga'] or float(c['domain_ga']) != model['official_domain_ga']: continue
        if float(c['sequence_bitscore']) < float(c['sequence_ga']) or float(c['domain_bitscore']) < float(c['domain_ga']): continue
        if float(d['model_coverage']) < minimum or float(c['model_coverage']) < minimum or float(d['query_core_overlap']) < model['min_query_core_overlap']: continue
        accepted.append((lo,hi,d))
    return accepted


def _site_map(row, pp_row, columns):
    result, pos = {}, 0
    for col, residue in enumerate(row,1):
        if residue not in '-.': pos += 1
        if col in columns:
            result[col] = dict(query_position=None if residue in '-.' else pos, residue=residue.upper(), pp=pp_row[col-1])
    return result


def _eligible(queries, features, model, parameters):
    return {q:ds[0] for q in sorted(features) if q in queries and features[q].get('family') == 'FeFe'
            and len(ds := _query_domain(features[q],model,parameters.get('min_domain_coverage'),q,len(queries[q]))) == 1}


def _decode(text, queries, features, eligible, opaque, model, ref):
    seqs, pps, rf = _stockholm(text)
    expected = {'ref_1FEH':ref['sequence'].upper(), **{opaque[q]:queries[q].replace('U','X') for q in eligible}}
    if set(seqs) != set(expected) or set(pps) != set(expected): raise ValueError('Stockholm identifiers changed')
    if any(len(row) != len(rf) or len(pps[q]) != len(rf) for q,row in seqs.items()): raise ValueError('Stockholm row lengths changed')
    if any(''.join(c for c in seqs[q] if c not in '-.').upper() != s for q,s in expected.items()): raise ValueError('Stockholm sequence reconstruction failed')
    matches = [i for i,c in enumerate(rf,1) if c.lower() == 'x']
    if len(matches) != 280: raise ValueError('PF02906 match count changed')
    cols = {p:matches[k-1] for p,k in zip(model['reference_positions'],model['expected_match_columns'])}
    refmap = _site_map(seqs['ref_1FEH'],pps['ref_1FEH'],set(cols.values()))
    if any(refmap[col] != dict(query_position=p,residue='C',pp='*') for p,col in cols.items()): raise ValueError('PF02906 reference ligand mapping failed')
    records = {}
    for q,(lo,hi,domain) in eligible.items():
        mapping = _site_map(seqs[opaque[q]],pps[opaque[q]],set(cols.values()))
        sites, blockers = [], []
        for rp,k in zip(model['reference_positions'],model['expected_match_columns']):
            site = dict(reference_position=rp,match_column=k,**mapping[cols[rp]])
            site['original_residue'] = queries[q][site['query_position']-1] if site['query_position'] is not None else ''
            state = 'not_covered' if site['query_position'] is None else 'ambiguous' if site['residue'] in AMBIGUOUS else 'present' if site['residue'] == 'C' else 'substituted'
            site.update(residue_status=state,posterior_status='high' if site['pp'] == '*' else 'low_or_missing')
            site['within_domain'] = site['query_position'] is not None and lo <= site['query_position'] <= hi
            site['status'] = 'supported' if state == 'present' and site['pp'] == '*' and site['within_domain'] else 'unresolved'
            if state != 'present': blockers.append(state)
            if site['pp'] != '*': blockers.append('low_profile_posterior')
            if not site['within_domain']: blockers.append('domain_interval_not_owned')
            sites.append(site)
        positions = [s['query_position'] for s in sites if s['query_position'] is not None]
        if len(positions) != len(set(positions)) or positions != sorted(positions): blockers.append('nonunique_or_unordered_positions')
        records[q] = dict(protein_id=q,profile_ligand_status='unresolved' if blockers else 'profile_mapped_ligand_pattern_supported',
            query_sha256=hashlib.sha256(queries[q].encode()).hexdigest(),query_length=len(queries[q]),feature_sha256=_digest(features[q]),
            domain_ownership=dict(accession='PF02906',query_start=lo,query_end=hi,domain_sha256=_digest(domain)),
            sites=sites,blockers=sorted(set(blockers)),rule_version=RULE_VERSION,hmm_sha256=model['hmm_sha256'],reference_model_sha256=model['reference_model_sha256'],
            normalization='U_to_X_same_length_only',scope=model['scope'])
    return records


def prepare(annotation, features, output, parameters=None, root=ROOT):
    annotation, output, root = Path(annotation).resolve(),Path(output).resolve(),Path(root).resolve()
    parameters = parameters or {}
    model,hmm,ref,refs = _load_model(root)
    proteins = annotation/'proteins.faa'
    queries = _fasta(proteins)
    eligible = _eligible(queries,features,model,parameters)
    if not eligible: return {},{}
    sources = [proteins,hmm,refs,MODEL_PATH,Path(__file__)]
    feature_file = annotation/'confirmation/features.json'
    if feature_file.exists():
        saved = {f['protein_id']:f for f in json.loads(feature_file.read_text())}
        if any(saved.get(q) != features[q] for q in eligible): raise ValueError('FeFe feature/source mismatch')
        sources.append(feature_file)
    inputs = {str(p.resolve()):_sha(p) for p in sources}
    output.mkdir(parents=True,exist_ok=False)
    opaque = {q:f'q_{i:08d}' for i,q in enumerate(eligible)}
    inp = output/'PF02906.input.faa'
    inp.write_text('>ref_1FEH\n'+ref['sequence']+'\n'+''.join('>'+opaque[q]+'\n'+queries[q].replace('U','X')+'\n' for q in eligible))
    (output/'PF02906.opaque_ids.json').write_text(json.dumps(opaque,sort_keys=True,indent=2)+'\n')
    hmmalign = root/'.env/bin/hmmalign'
    if not hmmalign.is_file():
        hmmalign = Path(shutil.which('hmmalign') or hmmalign).resolve()
    command = [str(hmmalign),'--outformat','stockholm',str(hmm),str(inp)]
    binary_hash = _sha(Path(command[0]))
    completed = subprocess.run(command,check=True,text=True,capture_output=True)
    (output/'PF02906.alignment.sto').write_text(completed.stdout)
    result = _decode(completed.stdout,queries,features,eligible,opaque,model,ref)
    from fefe_group_b_evidence import prepare as corroborate_assayed_lineage
    inputs.update(corroborate_assayed_lineage(annotation,output,result,features,queries,parameters,root))
    (output/'PF02906.query_mapping.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    outputs = {p.name:_sha(p) for p in output.iterdir() if p.is_file()}
    manifest = dict(status='complete',rule_version=RULE_VERSION,annotation=str(annotation),parameters=parameters,
        feature_sha256=_digest({q:features[q] for q in eligible}),input_sha256=inputs,output_sha256=outputs,
        command=command,hmmalign_binary_sha256=binary_hash,scope=model['scope'])
    (output/'PF02906.manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
    if any(_sha(Path(p)) != h for p,h in inputs.items()) or _sha(Path(command[0])) != binary_hash: raise ValueError('FeFe mapping input changed during execution')
    return result,{**inputs,**{str(output/n):h for n,h in outputs.items()},str(output/'PF02906.manifest.json'):_sha(output/'PF02906.manifest.json')}


def verify_saved(output, annotation=None, features=None, parameters=None, root=ROOT):
    output = Path(output).resolve()
    manifest = json.loads((output/'PF02906.manifest.json').read_text())
    if manifest.get('status') != 'complete' or manifest.get('rule_version') != RULE_VERSION: raise ValueError('Incomplete/stale FeFe mapping')
    if annotation is None or features is None or parameters is None: raise ValueError('FeFe verification requires actual input context')
    annotation = Path(annotation).resolve()
    if manifest['annotation'] != str(annotation) or manifest['parameters'] != parameters: raise ValueError('FeFe input scope/parameters changed')
    for p,h in manifest['input_sha256'].items():
        if _sha(Path(p)) != h: raise ValueError('FeFe mapping source changed')
    for n,h in manifest['output_sha256'].items():
        if Path(n).name != n or _sha(output/n) != h: raise ValueError('FeFe mapping output changed')
    model,hmm,ref,refs = _load_model(root)
    queries = _fasta(annotation/'proteins.faa')
    eligible = _eligible(queries,features,model,parameters)
    if manifest['feature_sha256'] != _digest({q:features[q] for q in eligible}): raise ValueError('FeFe mapping features changed')
    opaque = json.loads((output/'PF02906.opaque_ids.json').read_text())
    if set(opaque) != set(eligible) or len(set(opaque.values())) != len(opaque) or 'ref_1FEH' in opaque.values(): raise ValueError('FeFe opaque query ownership changed')
    expected = _decode((output/'PF02906.alignment.sto').read_text(),queries,features,eligible,opaque,model,ref)
    if (output/'group_B_assay_blast.tsv').is_file():
        from fefe_group_b_evidence import prepare as corroborate_assayed_lineage
        corroborate_assayed_lineage(annotation,output,expected,features,queries,parameters,root,replay=True)
    saved = json.loads((output/'PF02906.query_mapping.json').read_text())
    if saved != expected: raise ValueError('FeFe saved mapping contradicts actual alignment')
    return saved,manifest
