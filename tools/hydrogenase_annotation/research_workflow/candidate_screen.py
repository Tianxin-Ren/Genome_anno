"""Retain potential catalytic homologs without promoting identity or function."""
import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from Bio import SeqIO
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'research_workflow'))
from identity_scope import (NIFE_CORE_REGION_RULE, NIFE_CORE_REGION_STATUS,
                            nife_core_regions, nife_core_region_support, load_lineage_review)
if __package__ == 'research_workflow':
    from scripts.query_domain_organization import verify_completed
    from scripts.discovery_evidence import load_domains, load_targets, select_domains, search_discovery
    from scripts.annotate import BLAST_FIELDS, load_blast
    from research_workflow.dna_rescue import load_candidates as load_dna_candidates
    from research_workflow import preferred_prediction as _preferred_policy
    PREDICTION_FIELDS, predict = _preferred_policy.PREDICTION_FIELDS, _preferred_policy.predict
elif __package__:
    from ..scripts.query_domain_organization import verify_completed
    from ..scripts.discovery_evidence import load_domains, load_targets, select_domains, search_discovery
    from ..scripts.annotate import BLAST_FIELDS, load_blast
    from .dna_rescue import load_candidates as load_dna_candidates
    from . import preferred_prediction as _preferred_policy
    PREDICTION_FIELDS, predict = _preferred_policy.PREDICTION_FIELDS, _preferred_policy.predict
else:
    from query_domain_organization import verify_completed
    from discovery_evidence import load_domains, load_targets, select_domains, search_discovery
    from annotate import BLAST_FIELDS, load_blast
    from dna_rescue import load_candidates as load_dna_candidates
    import preferred_prediction as _preferred_policy
    PREDICTION_FIELDS, predict = _preferred_policy.PREDICTION_FIELDS, _preferred_policy.predict

PROFILES={'PF00374':'NiFe_related','PF00346':'NiFe_or_Nuo_related','PF02906':'FeFe_related','PF03201':'Hmd_related','PF22616':'Hmd_related'}
FIELDS='protein_id source_sets discovery_status family_candidates GA_profiles coordinates core_identity core_function rejected_reason potential_only domain_evidence raw_signal_evidence'.split()
# Preferred predictions are deliberately first: the legacy fields below remain
# byte-for-byte compatible in meaning and continue to carry the strict labels.
PREFERRED_CATALOG_FIELDS=list(PREDICTION_FIELDS)
CATALOG_FIELDS=['candidate_id', *PREFERRED_CATALOG_FIELDS, *'candidate_sources screen_status screen_families candidate_kind dna_locus_evidence identity_status function_prediction function_confidence supported_attributes rejected_reason identity_blockers evidence coarse_function_prediction mechanistic_subtype reference_reaction_status reference_contract_ids reference_reaction_keys reference_assay_reactions reference_assay_capabilities NAD_NADP_selectivity in_vivo_direction query_activity electron_role_context electron_role_context_scope'.split()]
from reference_reactions import TUPLE_POTENTIAL_STATUSES
from function_rules import FIELDS as FUNCTION_RULE_FIELDS
CATALOG_FIELDS += FUNCTION_RULE_FIELDS
REFERENCE_LABELS={
    'exact_reference_sequence_potential':'reference_sequence_matched_reaction_potential',
    'exact_reference_coding_locus_potential':'reference_coding_locus_matched_reaction_potential',
}


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def table(path):
    if not path.exists(): return []
    with path.open() as handle: return list(csv.DictReader(handle,delimiter='\t'))
def write(path,rows,fields):
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');writer.writeheader()
        writer.writerows({k:json.dumps(r.get(k),ensure_ascii=False) if isinstance(r.get(k),(list,dict)) else r.get(k,'') for k in fields} for r in rows)


def _preferred_prediction(row, evidence=None):
    """Project the shared preferred-prediction policy onto a catalog row."""
    evidence = evidence if isinstance(evidence, dict) else {}
    # Functional annotations are the validated producer.  Preserve their
    # bound preferred projection verbatim; only legacy tables use policy
    # fallback here.
    source = evidence.get('function') if isinstance(evidence.get('function'), dict) else evidence
    if isinstance(source, dict) and all(k in source and source.get(k) not in ('', None) for k in PREFERRED_CATALOG_FIELDS):
        return {k: source[k] for k in PREFERRED_CATALOG_FIELDS}
    result = predict(row, evidence=evidence)
    if not isinstance(result, dict):
        raise ValueError('preferred prediction policy returned a non-object')
    return {k: result.get(k, '') for k in PREFERRED_CATALOG_FIELDS}


def _preferred_value(value):
    """Normalize TSV JSON cells for provenance comparisons."""
    if isinstance(value, str):
        try:
            return json.loads(value)
        except (TypeError, ValueError):
            return value
    return value


def _validate_preferred_manifest(function_dir, fmanifest, function_rows):
    """Validate the bound preferred-prediction sidecar when present.

    R135-era manifests have no preferred keys and intentionally retain the
    legacy fallback path.  Once any preferred key is emitted, all provenance
    bindings are required and the TSV, sidecar, and manifest must agree.
    """
    keys = ('preferred_prediction_output_sha256', 'preferred_prediction_policy_version',
            'preferred_prediction_policy_sha256')
    present = [k for k in keys if k in fmanifest]
    if not present:
        return False
    if len(present) != len(keys):
        raise ValueError('Preferred prediction manifest binding is incomplete')
    sidecar = function_dir / 'preferred_predictions.json'
    if not sidecar.is_file() or sha(sidecar) != fmanifest['preferred_prediction_output_sha256']:
        raise ValueError('Preferred prediction output changed or is missing')
    try:
        predictions = json.loads(sidecar.read_text())
    except (OSError, TypeError, ValueError) as error:
        raise ValueError('Preferred prediction output is invalid JSON') from error
    if not isinstance(predictions, dict):
        raise ValueError('Preferred prediction output must be keyed by candidate ID')
    rows = {r.get('candidate_id'): r for r in function_rows}
    if None in rows or set(predictions) != set(rows):
        raise ValueError('Preferred prediction candidate IDs do not match function table')
    for candidate_id, row in rows.items():
        value = predictions[candidate_id]
        if not isinstance(value, dict) or set(PREDICTION_FIELDS) - set(value):
            raise ValueError('Preferred prediction fields are incomplete')
        for field in PREDICTION_FIELDS:
            if _preferred_value(row.get(field, '')) != value[field]:
                raise ValueError('Preferred prediction differs between function table and sidecar')
    # The policy source is part of the binding.  Validate the source digest
    # directly when the module is available; this catches hand-edited output
    # even when a policy implementation does not expose a constant.
    policy_path = ROOT / 'research_workflow/preferred_prediction.py'
    if policy_path.is_file() and sha(policy_path) != fmanifest['preferred_prediction_policy_sha256']:
        raise ValueError('Preferred prediction policy source changed')
    expected_version = getattr(_preferred_policy, 'VERSION', getattr(_preferred_policy, 'POLICY_VERSION', None))
    if not isinstance(expected_version, str) or not expected_version:
        raise ValueError('Preferred prediction policy version is unavailable')
    if fmanifest['preferred_prediction_policy_version'] != expected_version:
        raise ValueError('Preferred prediction policy version differs from current policy')
    if not isinstance(fmanifest['preferred_prediction_policy_version'], str) or not fmanifest['preferred_prediction_policy_version']:
        raise ValueError('Preferred prediction policy version is invalid')
    return True


def _raw_blast_hits(annotation, run, genes):
    """Keep reported catalytic/unresolved/decoy BLAST evidence without threshold promotion."""
    blast=annotation/'blast.tsv'
    if not blast.exists(): return [], None, {}
    references_path=ROOT/'references/processed/references.json'
    expected=run.get('references_json_sha256')
    if not expected or sha(references_path)!=expected:
        raise ValueError('Reference manifest does not bind references.json for raw BLAST evidence')
    references=json.loads(references_path.read_text())
    refs={r.get('reference_id'):r for r in references}
    proteins_path=annotation/'proteins.faa'
    reference_fasta=ROOT/'references/processed/references.faa'
    if not proteins_path.is_file() or not reference_fasta.is_file() or run.get('reference_fasta_sha256') != sha(reference_fasta):
        raise ValueError('Raw BLAST evidence lacks a bound query/reference sequence library')
    with proteins_path.open() as handle: queries={r.id:str(r.seq).upper().rstrip('*') for r in SeqIO.parse(handle,'fasta')}
    with reference_fasta.open() as handle: refs_seq={r.id:str(r.seq).upper().rstrip('*') for r in SeqIO.parse(handle,'fasta')}
    for line in blast.read_text().splitlines():
        fields=line.split('\t')
        if len(fields)!=len(BLAST_FIELDS): raise ValueError('Invalid BLAST row column count')
        if fields[0] not in genes or fields[0] not in queries or fields[1] not in refs or fields[1] not in refs_seq: raise ValueError('BLAST query/reference ownership is unknown')
        try:
            numeric={k:float(fields[BLAST_FIELDS.index(k)]) for k in ('pident','evalue','bitscore')}
            integer={k:int(fields[BLAST_FIELDS.index(k)]) for k in ('length','qstart','qend','sstart','send','qlen','slen')}
        except (TypeError,ValueError): raise ValueError('Invalid BLAST score or coordinate')
        if any(not math.isfinite(v) for v in numeric.values()) or not 0 <= numeric['pident'] <= 100 or numeric['evalue']<0 or numeric['bitscore']<0:
            raise ValueError('Invalid BLAST score or coordinate')
        if integer['qlen'] != len(queries[fields[0]]) or integer['slen'] != len(refs_seq[fields[1]]) or not 1<=integer['qstart']<=integer['qend']<=integer['qlen'] or not 1<=integer['sstart']<=integer['send']<=integer['slen']:
            raise ValueError('Invalid BLAST score or coordinate')
        qseq,sseq=fields[BLAST_FIELDS.index('qseq')],fields[BLAST_FIELDS.index('sseq')]
        if integer['length']<1 or len(qseq)!=integer['length'] or len(sseq)!=integer['length']:
            raise ValueError('Invalid BLAST alignment length')
        if len(qseq.replace('-','')) != integer['qend']-integer['qstart']+1 or len(sseq.replace('-','')) != integer['send']-integer['sstart']+1:
            raise ValueError('Invalid BLAST alignment coordinates')
        if qseq.replace('-','').upper() != queries[fields[0]][integer['qstart']-1:integer['qend']].upper():
            raise ValueError('BLAST query sequence does not match bound protein')
        if sseq.replace('-','').upper() != refs_seq[fields[1]][integer['sstart']-1:integer['send']].upper():
            raise ValueError('BLAST reference sequence does not match bound library')
    # Reuse the canonical parser for coverage/typed fields, then restore the
    # kind filter so module/maturation references never enter this queue.
    parsed=load_blast(blast,run['thresholds'],references)
    return [h for h in parsed if refs[h['sseqid']].get('kind') in {'catalytic','unresolved','decoy'}], sha(blast), refs


def _raw_target_hits(annotation, genes, path=None):
    """Read full-sequence HMM hits when hmmsearch emitted a target table."""
    path=Path(path) if path is not None else annotation/'family.tbl'
    if not path.exists(): return [], None
    manifest=json.loads((ROOT/'references/processed/profile_manifest.json').read_text())
    families={r['profile']:r['family'] for r in manifest}
    rows=load_targets(path,families)
    if any(r['protein_id'] not in genes for r in rows):
        raise ValueError('Whole-target candidate lacks source gene ownership')
    return rows, sha(path)


def _recover_target_table(annotation, run, out):
    """Re-run the source proteins into a new output-only cache.

    Legacy annotations can predate the target table. Recovery is bound to the
    same model/profile/policy and records the current implementation and
    executable; it never writes into the annotation directory.
    """
    model=ROOT/'references/processed/catalytic_families.hmm'
    profile_manifest=ROOT/'references/processed/profile_manifest.json'
    policy=ROOT/'config/discovery_evidence.json'
    script=ROOT/'scripts/discovery_evidence.py'
    from tool_paths import resolve_tool
    binary=resolve_tool(ROOT, 'hmmsearch')
    expected={
        'model_sha256': sha(model),
        'discovery_profile_manifest_sha256': sha(profile_manifest),
        'discovery_evidence_policy_sha256': sha(policy),
    }
    script_sha=sha(script)
    for key,value in expected.items():
        if run.get(key) != value:
            raise ValueError('Source discovery binding does not match current '+key)
    if not binary.is_file():
        raise ValueError('Current hmmsearch binary is unavailable for family-target recovery')
    recovery=out/'family_target_recovery'
    recovery.mkdir(parents=True,exist_ok=False)
    manifest_path=recovery/'manifest.json'
    manifest=dict(status='running',source_kind='recovered_table',source_annotation=str(annotation),
        source_run_manifest_sha256=sha(annotation/'run_manifest.json'),
        source_proteins=str(annotation/'proteins.faa'),source_proteins_sha256=sha(annotation/'proteins.faa'),
        model=str(model),model_sha256=expected['model_sha256'],profile_manifest_sha256=expected['discovery_profile_manifest_sha256'],
        policy_sha256=expected['discovery_evidence_policy_sha256'],discovery_evidence_script_sha256=script_sha,
        hmmsearch_binary=str(binary),hmmsearch_binary_sha256=sha(binary),output_sha256={})
    manifest['source_family_tbl_absent_before']=not (annotation/'family.tbl').exists()
    bound_paths={
        'source_run_manifest': annotation/'run_manifest.json',
        'source_proteins': annotation/'proteins.faa',
        'model': model,
        'profile_manifest': profile_manifest,
        'policy': policy,
        'discovery_evidence_script': script,
        'hmmsearch_binary': binary,
    }
    manifest['bound_input_sha256']={name:sha(path) for name,path in bound_paths.items()}
    manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    try:
        threads=run.get('threads',2)
        if not isinstance(threads,int) or isinstance(threads,bool) or threads<1: threads=2
        thresholds=run.get('thresholds')
        if not isinstance(thresholds,dict): raise ValueError('Source run lacks discovery thresholds')
        search_discovery(annotation/'proteins.faa',recovery,ROOT,threads,thresholds)
        if (annotation/'family.tbl').exists():
            raise ValueError('Source family.tbl appeared during recovery scan')
        after={name:sha(path) for name,path in bound_paths.items()}
        if after != manifest['bound_input_sha256']:
            raise ValueError('Recovery-bound input changed during scan')
        outputs=[recovery/name for name in ('family.domtbl','family.tbl','hmmsearch.log')]
        if any(not p.is_file() for p in outputs): raise ValueError('Family-target recovery did not produce all outputs')
        manifest['status']='complete'
        manifest['output_sha256']={str(p.relative_to(recovery)):sha(p) for p in outputs}
    except Exception as exc:
        manifest['status']='failed';manifest['error']=f'{type(exc).__name__}: {exc}'
        manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
        raise
    manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    return recovery,manifest


def _load_dna_rescue_candidates(dna_rescue_dir, annotation=None, screen_dir=None, genome=None, protein_ids=()):
    """Validate DNA rescue provenance and return unresolved locus records only."""
    dna_rescue_dir=Path(dna_rescue_dir).resolve()
    manifest,candidates=load_dna_candidates(dna_rescue_dir,annotation=annotation,screen_dir=screen_dir,genome=genome,root=ROOT)
    status=manifest.get('status')
    if status != 'complete':
        if status not in {'not_applicable','not_evaluated'}:
            raise ValueError('DNA rescue is not complete: '+str(status))
        if candidates or manifest.get('candidate_count',0)!=0:
            raise ValueError('Non-complete DNA rescue has candidate rows')
        return [],manifest
    ids=[c.get('candidate_id') for c in candidates]
    if any(not isinstance(i,str) or not i.startswith('DNA_') for i in ids) or len(set(ids))!=len(ids) or set(ids)&set(protein_ids):
        raise ValueError('DNA rescue candidate ID collision or invalid ID')
    rows=[]
    for c in candidates:
        if c.get('identity_status')!='unconfirmed' or c.get('function_prediction')!='potential_homolog_function_unresolved':
            raise ValueError('DNA rescue candidate contains an upgraded identity or function label')
        try:
            contig=c['contig'];start=int(c['start']);end=int(c['end']);strand=int(c['strand']);frame=int(c['frame'])
        except (KeyError,TypeError,ValueError) as error: raise ValueError('Invalid DNA locus evidence') from error
        if not isinstance(contig,str) or not contig or start<1 or end<start or strand not in {-1,1} or abs(frame) not in {1,2,3} or (1 if frame>0 else -1)!=strand:
            raise ValueError('Invalid DNA locus evidence')
        tables=c.get('translation_tables')
        if not isinstance(tables,list) or not tables or any(not isinstance(x,int) or isinstance(x,bool) or x not in {1,4,6,9,10,11,12,13,14,15,16,21,22,23} for x in tables) or len(set(tables))!=len(tables): raise ValueError('Invalid DNA translation tables')
        if not isinstance(c.get('DNA_sha256'),str) or len(c['DNA_sha256'])!=64 or any(x not in '0123456789abcdefABCDEF' for x in c['DNA_sha256']): raise ValueError('Invalid DNA sequence hash')
        if not isinstance(c.get('translation_sha256'),str) or len(c['translation_sha256'])!=64 or any(x not in '0123456789abcdefABCDEF' for x in c['translation_sha256']): raise ValueError('Invalid translation hash')
        if not isinstance(c.get('source_hsp_ids'),list) or not c['source_hsp_ids'] or any(not isinstance(x,str) or not x for x in c['source_hsp_ids']): raise ValueError('Invalid DNA source HSP IDs')
        rows.append(dict(candidate_id=c['candidate_id'],candidate_sources=['dna_translated_homology'],screen_status='DNA_locus_unresolved',
            screen_families=c.get('families',[]),identity_status='unconfirmed',function_prediction='potential_homolog_function_unresolved',
            function_confidence='unresolved',supported_attributes=[],rejected_reason='',identity_blockers=[],candidate_kind='dna_locus',
            dna_locus_evidence=dict(scope=c.get('scope','translated_DNA_interval_not_reconstructed_CDS'),contig=contig,start=start,end=end,strand=strand,frame=frame,translation_tables=tables,
                                    reference_ids=c.get('reference_ids',[]),families=c.get('families',[]),profile_evidence=c.get('profile_evidence',[]),
                                    source_hsp_ids=c.get('source_hsp_ids',[]),translation_has_stop=bool(c.get('translation_has_stop',False)),
                                    translation_stop_positions=c.get('translation_stop_positions',[]),DNA_sha256=c.get('DNA_sha256',''),translation_sha256=c.get('translation_sha256','')),
            evidence='dna_rescue/dna_rescue_candidates.json;dna_rescue/dna_rescue_manifest.json',coarse_function_prediction='potential_homolog_function_unresolved',
            reference_reaction_status='not_assessed',reference_contract_ids=[],reference_reaction_keys=[],reference_assay_reactions=[],
            reference_assay_capabilities=[],NAD_NADP_selectivity='unresolved',in_vivo_direction='unknown',query_activity='not_tested',
            electron_role_context={},electron_role_context_scope='not_assessed'))
        rows[-1]['identity_blockers']=['DNA_interval_not_validated_CDS','catalytic_identity_not_assessed']
        rows[-1].update(_preferred_prediction(rows[-1], {'dna_rescue': c, 'screen_status': 'DNA_locus_unresolved'}))
    return rows,manifest


def catalog(screen_dir,function_dir,output,confirmation_dir=None,dna_rescue_dir=None):
    """Keep every screened homolog visible, with only its actual supported function."""
    screen_dir,function_dir,output=map(Path,[screen_dir,function_dir,output])
    report=json.loads((screen_dir/'candidate_screen_manifest.json').read_text())
    if any(sha(screen_dir/p)!=d for p,d in report['output_sha256'].items()): raise ValueError('Candidate screen changed')
    fmanifest=json.loads((function_dir/'manifest.json').read_text())
    if sha(function_dir/'functional_annotations.tsv')!=fmanifest['functional_annotations_sha256']: raise ValueError('Function table changed')
    source_manifest=str(Path(report['annotation'])/'run_manifest.json')
    if source_manifest not in fmanifest['input_sha256'] or fmanifest['input_sha256'][source_manifest]!=report['input_sha256'].get(source_manifest): raise ValueError('Catalog source annotation mismatch')
    function_rows=table(function_dir/'functional_annotations.tsv')
    preferred_bound = _validate_preferred_manifest(function_dir, fmanifest, function_rows)
    class_predictions={}
    if fmanifest.get('function_rules_manifest_sha256'):
        from function_rules import load_predictions
        if sha(function_dir/'function_rules/manifest.json')!=fmanifest['function_rules_manifest_sha256']:
            raise ValueError('Functional rule manifest changed')
        class_predictions=load_predictions(function_dir/'function_rules',report['annotation'])
        if set(class_predictions)!={r['candidate_id'] for r in function_rows}:
            raise ValueError('Functional rule candidates differ from function table')
        for row in function_rows:
            for field in FUNCTION_RULE_FIELDS:
                value=json.loads(row[field]) if field in {'function_class_ids','function_rule_alternatives'} else row[field]
                if value!=class_predictions[row['candidate_id']][field]:raise ValueError('Functional rule projection differs: '+field)
    role_contexts={}
    if fmanifest.get('module_role_inventory_sha256'):
        role_path=function_dir/'independent_roles/inventory.json'
        if not role_path.is_file() or sha(role_path)!=fmanifest['module_role_inventory_sha256']:raise ValueError('Independent module-role evidence changed or missing')
        role_inventory=json.loads(role_path.read_text());role_contexts=role_inventory['candidate_contexts']
        if not isinstance(role_contexts,dict):raise ValueError('Independent module-role candidate contexts are invalid')
        if role_inventory.get('input_sha256',{}).get(source_manifest)!=report['input_sha256'].get(source_manifest):raise ValueError('Module-role context belongs to another annotation')
    functions={r['candidate_id']:r for r in function_rows}
    if len(functions)!=len(function_rows): raise ValueError('Repeated function candidate ID')
    decisions={}
    if confirmation_dir is not None:
        confirmation_dir=Path(confirmation_dir)
        lock=json.loads((confirmation_dir/'candidate_confirmation_manifest.json').read_text())
        if sha(confirmation_dir/'candidate_confirmation.json')!=lock['report_sha256']: raise ValueError('Additional-candidate confirmation changed')
        confirmation=json.loads((confirmation_dir/'candidate_confirmation.json').read_text())
        if confirmation['source_screen_manifest_sha256']!=sha(screen_dir/'candidate_screen_manifest.json'): raise ValueError('Confirmation belongs to another screen')
        decisions={r['candidate_id']:r for r in confirmation['decisions']}
        region_inputs = confirmation.get('NiFe_core_region_input_sha256', {})
        for path, digest in region_inputs.items():
            if not Path(path).is_file() or sha(Path(path)) != digest:
                raise ValueError('NiFe core-region proof input changed or missing')
        for decision in decisions.values():
            if decision.get('sequence_status') == NIFE_CORE_REGION_STATUS:
                if not region_inputs:
                    raise ValueError('NiFe core-region proof lacks bound inputs')
                proof = decision.get('NiFe_core_region')
                if (not isinstance(proof, dict) or proof.get('status') != NIFE_CORE_REGION_STATUS or
                        proof.get('rule_version') != NIFE_CORE_REGION_RULE or
                        proof.get('protein_id') != decision.get('candidate_id')):
                    raise ValueError('NiFe core-region decision lacks a valid bound proof')
                try:
                    feature_path = confirmation['evidence'][decision['candidate_id']]['features']
                    feature_rows = json.loads(Path(feature_path).read_text())
                    feature = next(x for x in feature_rows if x.get('protein_id') == decision['candidate_id'])
                    genes = {x['id']: x for x in table(Path(report['annotation']) / 'genes.tsv')}
                    from Bio import SeqIO
                    with (Path(report['annotation']) / 'proteins.faa').open() as handle:
                        protein = next(x for x in SeqIO.parse(handle, 'fasta') if x.id == decision['candidate_id'])
                    gene = dict(genes[decision['candidate_id']], protein_sha256=hashlib.sha256(str(protein.seq).upper().encode()).hexdigest())
                    recomputed, _ = nife_core_regions(Path(feature_path).resolve().parents[1], {decision['candidate_id']: feature}, ROOT)
                    if recomputed.get(decision['candidate_id']) != proof:
                        raise ValueError('NiFe core-region proof recomputation differs')
                    parameters = json.loads((ROOT / 'config/catalytic_confirmation_thresholds.json').read_text())['families']['NiFe']
                    if not nife_core_region_support(feature, {'NiFe_core_region': proof}, parameters, load_lineage_review(ROOT), gene):
                        raise ValueError('NiFe core-region proof does not validate against current feature')
                except (KeyError, StopIteration, OSError, ValueError) as error:
                    raise ValueError('NiFe core-region proof does not validate against current feature') from error
        context=fmanifest.get('candidate_context_metadata',{})
        if context.get('confirmation') and Path(context['confirmation']).resolve()!=confirmation_dir.resolve(): raise ValueError('Function manifest confirmation context mismatch')
        if context.get('confirmation_manifest_sha256') and context['confirmation_manifest_sha256'] != sha(confirmation_dir/'candidate_confirmation_manifest.json'): raise ValueError('Function manifest confirmation hash mismatch')
        if any(r['identity_status']=='calibrated_high' or r['function_status']!='unassigned' for r in decisions.values()): raise ValueError('Additional discovery cannot promote calibrated identity or function')
    screened=table(screen_dir/'candidate_screen.tsv')
    if set(role_contexts)-{r['protein_id'] for r in screened}:raise ValueError('Module-role context absent from candidate screen')
    if set(functions)-{r['protein_id'] for r in screened}: raise ValueError('Function candidate absent from high-recall screen')
    core_ids={r['protein_id'] for r in table(screen_dir/'candidate_screen.tsv') if r.get('discovery_status')=='existing_anchor'}
    extra_function_ids=set(functions)-core_ids
    if extra_function_ids and not confirmation_dir: raise ValueError('Extra function rows require the matching candidate confirmation')
    if extra_function_ids-set(decisions): raise ValueError('Extra function row lacks confirmed candidate decision')
    if extra_function_ids:
        context=fmanifest.get('candidate_context_metadata',{})
        if not context.get('confirmation') or Path(context['confirmation']).resolve()!=confirmation_dir.resolve(): raise ValueError('Extra function confirmation context missing or mismatched')
        if context.get('confirmation_manifest_sha256')!=sha(confirmation_dir/'candidate_confirmation_manifest.json'): raise ValueError('Extra function confirmation binding missing or changed')
        if set(fmanifest.get('candidate_scoped_research_predictions',[]))!=extra_function_ids: raise ValueError('Extra function context IDs missing or mismatched')
        if any(decisions[q].get('sequence_status') not in {'sequence_supported_uncalibrated_core','sequence_supported_uncalibrated_noncanonical'} or decisions[q].get('blocking_reasons') for q in extra_function_ids): raise ValueError('Unconfirmed candidate cannot acquire a function row')
    if any(functions[q].get('identity_status') != decisions[q].get('sequence_status') for q in extra_function_ids): raise ValueError('Extra function identity bypassed confirmation sequence status')
    reference_path=function_dir/'reference_reaction_evidence.json'
    has_reference=any(f.get('reference_reaction_status','not_assessed')!='not_assessed' for f in functions.values())
    if has_reference:
        if not reference_path.is_file() or sha(reference_path)!=fmanifest.get('reference_reaction_output_sha256'):raise ValueError('Reference reaction evidence changed or missing')
        references=json.loads(reference_path.read_text())
        if set(references)-set(functions):raise ValueError('Reference reaction evidence has an unknown function candidate')
        for q,f in functions.items():
            e=references.get(q,{});caps=e.get('capabilities',[])
            if f.get('reference_reaction_status','not_assessed')!=e.get('status','not_assessed'):raise ValueError('Reference reaction status differs from evidence')
            expected={'reference_assay_capabilities':caps,'reference_contract_ids':[c['contract_id'] for c in caps],
                'reference_reaction_keys':[c['reaction_key'] for c in caps],'reference_assay_reactions':[c['molecular_reaction'] for c in caps]}
            if any(json.loads(f.get(k,'[]'))!=v for k,v in expected.items()):raise ValueError('Reference reaction detail differs from evidence')
    for f in functions.values():
        status=f.get('reference_reaction_status')
        if status not in TUPLE_POTENTIAL_STATUSES and f.get('function_prediction') in REFERENCE_LABELS.values():raise ValueError('Reference main label lacks its positive status')
    rows=[]
    for r in screened:
        f=functions.get(r['protein_id'],{});decision=decisions.get(r['protein_id'],{})
        identity=f.get('identity_status',decision.get('identity_status','not_confirmed'))
        region_status = identity == NIFE_CORE_REGION_STATUS or decision.get('sequence_status') == NIFE_CORE_REGION_STATUS
        if region_status and f and r['protein_id'] not in core_ids:
            raise ValueError('NiFe core-region candidate cannot acquire a function row')
        default_label=('NiFe_catalytic_core_region_supported_function_unresolved' if region_status else
                       'hydrogenase_family_supported_function_unresolved' if identity in ['sequence_supported_uncalibrated_core','sequence_supported_uncalibrated_noncanonical'] else 'potential_homolog_function_unresolved')
        supported = f.get('supported_attributes','[]')
        if region_status:
            parsed_supported = json.loads(supported) if isinstance(supported,str) else list(supported or [])
            if 'NiFe_catalytic_core_region_potential' not in parsed_supported:
                parsed_supported.append('NiFe_catalytic_core_region_potential')
            supported = parsed_supported
        function_prediction = default_label if region_status else f.get('function_prediction',default_label)
        coarse_prediction = default_label if region_status else f.get('coarse_function_prediction',f.get('function_prediction',default_label))
        reference_status = f.get('reference_reaction_status','not_assessed') if f else 'not_assessed'
        if region_status and f:
            if (f.get('function_prediction') != default_label or
                    f.get('coarse_function_prediction', f.get('function_prediction')) != default_label or
                    f.get('confidence', 'unresolved') != 'uncalibrated_core_region' or
                    reference_status in TUPLE_POTENTIAL_STATUSES):
                raise ValueError('Native NiFe core-region function row has an incompatible label')
            for key in ('reference_contract_ids','reference_reaction_keys','reference_assay_reactions','reference_assay_capabilities'):
                value = f.get(key, '[]')
                if json.loads(value) if isinstance(value, str) else value:
                    raise ValueError('Native NiFe core-region function row has reaction capabilities')
        rows.append(dict(candidate_id=r['protein_id'],candidate_sources=r['source_sets'],screen_status=r['discovery_status'],screen_families=r['family_candidates'],candidate_kind='protein',dna_locus_evidence={},
            identity_status=identity,function_prediction=function_prediction,
            coarse_function_prediction=coarse_prediction,
            mechanistic_subtype='' if region_status else f.get('mechanistic_subtype',''),
            reference_reaction_status=reference_status,
            reference_contract_ids='[]' if region_status else f.get('reference_contract_ids','[]'),reference_reaction_keys='[]' if region_status else f.get('reference_reaction_keys','[]'),
            reference_assay_reactions=f.get('reference_assay_reactions','[]'),reference_assay_capabilities=f.get('reference_assay_capabilities','[]'),
            function_confidence='uncalibrated_core_region' if region_status else f.get('confidence','unresolved'),supported_attributes=supported,
            NAD_NADP_selectivity=f.get('NAD_NADP_selectivity','unresolved'),in_vivo_direction=f.get('in_vivo_direction','unknown'),query_activity=f.get('query_activity','not_tested'),
            rejected_reason=r['rejected_reason'],identity_blockers=decision.get('blocking_reasons',[]),
            electron_role_context=role_contexts.get(r['protein_id'],{}),
            electron_role_context_scope=role_contexts.get(r['protein_id'],{}).get('scope','not_assessed'),
            evidence=('candidate_confirmation/candidate_confirmation.json;' if r['protein_id'] in extra_function_ids else '')+'candidates/candidate_screen.tsv'+(';functions/reference_reaction_evidence.json' if f.get('reference_reaction_status','not_assessed')!='not_assessed' else '')+(';functions/functional_annotations.tsv' if f else ';candidate_confirmation/candidate_confirmation.json' if decision else '')+(';functions/independent_roles/inventory.json' if r['protein_id'] in role_contexts else '')))
        preferred_evidence = dict(f) if isinstance(f, dict) else {}
        preferred_evidence.update({'function': f, 'screen': r, 'decision': decision,
                                   'role_context': role_contexts.get(r['protein_id'], {})})
        rows[-1].update(_preferred_prediction(rows[-1], preferred_evidence))
        if r['protein_id'] in class_predictions:rows[-1].update(class_predictions[r['protein_id']])
    if dna_rescue_dir is not None:
        dna_rows,dna_manifest=_load_dna_rescue_candidates(dna_rescue_dir,Path(report['annotation']).resolve(),screen_dir.resolve(),set(functions)|{r['protein_id'] for r in screened})
        rows.extend(dna_rows)
    if output.exists(): raise FileExistsError('Use a fresh catalog output path')
    write(output,rows,CATALOG_FIELDS)
    return rows


def dna_only_catalog(dna_rescue_dir, output, genome):
    """Write a DNA-only unresolved catalog without an annotation/function stage."""
    output=Path(output).resolve()
    if output.exists(): raise FileExistsError('Use a fresh catalog output path')
    rows,_manifest=_load_dna_rescue_candidates(dna_rescue_dir,genome=Path(genome).resolve())
    write(output,rows,CATALOG_FIELDS)
    return rows


def screen(annotation,domains,out):
    annotation,domains,out=map(lambda p:Path(p).resolve(),[annotation,domains,out])
    if out.exists(): raise FileExistsError('Use a fresh candidate-screen output directory')
    dman=verify_completed(domains)
    if str(annotation) not in dman['source_runs']: raise ValueError('Domain evidence belongs to another annotation source')
    run=json.loads((annotation/'run_manifest.json').read_text())
    resources={ROOT/'references/processed/profile_manifest.json':run['discovery_profile_manifest_sha256'],ROOT/'config/discovery_evidence.json':run['discovery_evidence_policy_sha256']}
    if any(sha(p)!=d for p,d in resources.items()): raise ValueError('Saved discovery resources changed')
    out.mkdir(parents=True)
    inputs=[annotation/n for n in ['run_manifest.json','genes.tsv','proteins.faa','catalytic_anchors.tsv','family_candidates.tsv','rejected.tsv','family.domtbl','family.tbl','blast.tsv'] if (annotation/n).exists()]
    from tool_paths import resolve_tool
    inputs += [domains/'query_domain_manifest.json',domains/'query_domain_organizations.tsv',Path(__file__),
               ROOT/'scripts/discovery_evidence.py',resolve_tool(ROOT, 'hmmsearch'),
               ROOT/'references/processed/catalytic_families.hmm',*resources,
               ROOT/'references/processed/references.json',ROOT/'references/processed/references.faa']
    target_table=domains/'query_target_signals.tsv'
    if target_table.exists(): inputs.append(target_table)
    before={str(p):sha(p) for p in inputs}
    anchors={r['protein_id']:r for r in table(annotation/'catalytic_anchors.tsv')}
    rejected={r['protein_id']:r for r in table(annotation/'rejected.tsv')}
    genes={r['id']:r for r in table(annotation/'genes.tsv')}
    with (annotation/'proteins.faa').open() as handle: protein_lengths={r.id:len(r.seq) for r in SeqIO.parse(handle,'fasta')}
    if set(protein_lengths) != set(genes): raise ValueError('Candidate source proteins and genes do not have identical ownership')
    raw_blast,raw_blast_sha,refs=_raw_blast_hits(annotation,run,genes)
    recovery_manifest=None
    target_source='saved_table'
    target_path=annotation/'family.tbl'
    domain_source='saved_table'
    domain_path=annotation/'family.domtbl'
    if not target_path.exists():
        recovery_dir,recovery_manifest=_recover_target_table(annotation,run,out)
        target_path=recovery_dir/'family.tbl'
        domain_path=recovery_dir/'family.domtbl'
        target_source='recovered_table'
        domain_source='recovered_table'
        # Bind every recovered artifact into the source-integrity check and
        # candidate-screen manifest.  The source annotation remains untouched.
        inputs += [recovery_dir/name for name in ('family.domtbl','family.tbl','hmmsearch.log','manifest.json')]
        # Keep the pre-recovery hashes for source inputs immutable. Only add
        # hashes for the newly created recovery artifacts.
        before.update({str(p):sha(p) for p in inputs if p.exists() and str(p) not in before})
    raw_targets,raw_target_sha=_raw_target_hits(annotation,genes,target_path)
    auxiliary_targets=[]
    if target_table.exists():
        for r in table(target_table):
            if r.get('source_run') != str(annotation): continue
            if r.get('protein_id') not in genes: raise ValueError('Auxiliary whole-target query lacks source gene ownership')
            accession=r.get('accession','').split('.')[0]
            if accession in PROFILES:
                auxiliary_targets.append(dict(r, family=PROFILES[accession], profile=r.get('profile',accession), accession=accession))
    family={}
    for r in table(annotation/'family_candidates.tsv'): family.setdefault(r['protein_id'],[]).append(r)
    ga={}
    for row in table(domains/'query_domain_organizations.tsv'):
        if row['source_run']!=str(annotation): continue
        evidence=[e for e in json.loads(row['profile_evidence']) if e['accession'].split('.')[0] in PROFILES and e['supported'] is True]
        for e in evidence:
            e['query_protein_id']=row['protein_id']
            for part in e['components']: part['protein_id']=row['protein_id']
        if evidence: ga[row['protein_id']]=evidence
    rows=[]
    for pid in sorted(set(anchors)|set(rejected)|set(family)|set(ga)):
        if pid not in genes: raise ValueError('Candidate lacks source gene ownership: '+pid)
        a=anchors.get(pid,{});gene=genes[pid];evidence=ga.get(pid,[])
        sources=[name for name,pool in [('existing_anchor',anchors),('family_candidate',family),('verified_GA_domain',ga),('core_rejected_retained',rejected)] if pid in pool]
        coordinate={k:gene.get(k,'') for k in ['contig','start','end','strand','bin_id','coordinate_evidence']}
        if run.get('mode')=='transcript' or coordinate['coordinate_evidence']=='unavailable_translated_transcript':
            coordinate.update(start='',end='',strand='')
        rows.append(dict(protein_id=pid,source_sets=sources,discovery_status='existing_anchor' if a else 'noncatalytic_homolog_favored_retained' if pid in rejected else 'discovery_only',
            family_candidates=sorted({r['family'] for r in family.get(pid,[])}|{PROFILES[e['accession'].split('.')[0]] for e in evidence}),
            GA_profiles=sorted({e['accession'] for e in evidence}),coordinates=coordinate,core_identity=a.get('catalytic_identity_confidence','not_assessed'),
            core_function=a.get('function_status','unassigned'),rejected_reason=rejected.get(pid,{}).get('reason',''),potential_only=not bool(a),domain_evidence=evidence))
    ledger=[]
    if domain_path.exists():
        manifest=json.loads((ROOT/'references/processed/profile_manifest.json').read_text())
        raw=load_domains(domain_path,{r['profile']:r['family'] for r in manifest})
        policy=json.loads((ROOT/'config/discovery_evidence.json').read_text());thresholds=run['thresholds']
        selected=select_domains(raw,thresholds,policy)
        def key(d): return tuple(d[k] for k in ['protein_id','profile','model_start','model_end','query_start','query_end'])
        selected_parts={key(d) for r in selected for d in json.loads(r['component_domains'])}
        for d in raw:
            if d['protein_id'] not in genes: raise ValueError('Unknown raw-domain query')
            if d['query_length'] != protein_lengths[d['protein_id']]: raise ValueError('Raw-domain query length differs from bound protein')
            cov=(d['model_end']-d['model_start']+1)/d['model_length'];reasons=[]
            if d['independent_evalue']>thresholds['hmm_domain_evalue']: reasons.append('above_saved_E_threshold')
            if cov<thresholds['hmm_coverage']: reasons.append('individual_domain_below_saved_coverage')
            selected_part=key(d) in selected_parts
            if selected_part: reasons.append('retained_as_single_domain_or_chain_component')
            elif not reasons: reasons.append('not_retained_by_saved_chain_rule')
            ledger.append(dict(d,profile_coverage=cov,retained_by_core_discovery=selected_part,filter_reason=';'.join(reasons)))
    # Keep strict candidates first; extend recall with every valid HMMER-reported
    # raw hit that was otherwise absent. This does not change select_domains().
    existing={r['protein_id'] for r in rows}
    raw_by_protein={}
    for d in ledger: raw_by_protein.setdefault(d['protein_id'],[]).append(d)
    blast_by_protein={}
    for h in raw_blast: blast_by_protein.setdefault(h['qseqid'],[]).append(h)
    target_by_protein={}
    for h in raw_targets: target_by_protein.setdefault(h['protein_id'],[]).append(h)
    for h in auxiliary_targets: target_by_protein.setdefault(h['protein_id'],[]).append(h)
    for row in rows:
        pid=row['protein_id']
        if pid in raw_by_protein or pid in blast_by_protein or pid in target_by_protein:
            row['raw_signal_evidence']={'family_domain':raw_by_protein.get(pid,[]),
                                        'reference_blast':blast_by_protein.get(pid,[]),
                                        'whole_target':target_by_protein.get(pid,[])}
    for pid in sorted(set(raw_by_protein)|set(blast_by_protein)|set(target_by_protein)):
        hits=raw_by_protein.get(pid,[]); blast_hits=blast_by_protein.get(pid,[])
        if pid in existing: continue
        if pid not in genes: raise ValueError('Raw-domain candidate lacks source gene ownership: '+pid)
        gene=genes[pid]
        coordinate={k:gene.get(k,'') for k in ['contig','start','end','strand','bin_id','coordinate_evidence']}
        if run.get('mode')=='transcript' or coordinate['coordinate_evidence']=='unavailable_translated_transcript':
            coordinate.update(start='',end='',strand='')
        rows.append(dict(protein_id=pid,source_sets=['raw_discovery_signal'],
            discovery_status='weak_or_partial_signal',
            family_candidates=sorted({d['family'] for d in hits}|{refs.get(h['sseqid'],{}).get('center','') for h in blast_hits}|{h['family'] for h in target_by_protein.get(pid,[])}-{''}), GA_profiles=[],
            coordinates=coordinate, core_identity='not_assessed', core_function='unassigned',
            rejected_reason='', potential_only=True, domain_evidence=hits,
            raw_signal_evidence={'family_domain':hits,'reference_blast':blast_hits,
                                 'whole_target':target_by_protein.get(pid,[])}))
    if recovery_manifest and (annotation/'family.tbl').exists():
        recovery_manifest_path=out/'family_target_recovery/manifest.json'
        failed=json.loads(recovery_manifest_path.read_text())
        failed.update(status='failed',error='Source family.tbl appeared before screen completion')
        recovery_manifest_path.write_text(json.dumps(failed,ensure_ascii=False,indent=2)+'\n')
        raise ValueError('Source family.tbl appeared before screen completion')
    if before!={p:sha(p) for p in before}: raise ValueError('Candidate-screen inputs changed')
    write(out/'candidate_screen.tsv',rows,FIELDS)
    write(out/'family_domain_ledger.tsv',ledger,['protein_id','family','profile','independent_evalue','domain_bitscore','query_length','query_start','query_end','model_length','model_start','model_end','profile_coverage','retained_by_core_discovery','filter_reason'])
    report=dict(status='complete',annotation=str(annotation),domains=str(domains),input_sha256=before,candidate_count=len(rows),GA_candidate_count=len(ga),
        raw_signal_counts=dict(family_domain=len(ledger),reference_blast=len(raw_blast),whole_target=len(raw_targets)),
        whole_target_scope='evaluated_and_bound',whole_target_source=target_source,domain_source=domain_source,
        family_target_recovery_manifest_sha256=sha(out/'family_target_recovery/manifest.json') if recovery_manifest else None,
        added_outside_core_anchor_count=sum(r['protein_id'] not in anchors for r in rows),raw_domain_rows=len(ledger),
        output_sha256={str(p.relative_to(out)):sha(p) for p in [out/'candidate_screen.tsv',out/'family_domain_ledger.tsv'] +
                       ([out/'family_target_recovery'/name for name in ('family.domtbl','family.tbl','hmmsearch.log','manifest.json')] if recovery_manifest else [])},
        scope='High-recall queue over structurally valid HMMER-reported domains, selected whole-target core profiles, and reported catalytic/unresolved/decoy BLAST homologs. Strict discovery thresholds remain unchanged; weak raw signals never promote identity or function, and unreported sequences remain outside this queue.')
    (out/'candidate_screen_manifest.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotation',type=Path,required=True);parser.add_argument('--domains',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();print(json.dumps(screen(args.annotation,args.domains,args.out),ensure_ascii=False))
