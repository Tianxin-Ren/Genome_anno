#!/usr/bin/env python3
"""Actual source-protein controls for identity/function separation, negatives and incomplete sites."""
import argparse
import copy
import json
from pathlib import Path
from annotate import ROOT,annotate,extract_genbank,write_tsv
from catalytic_confirmation import sha,table


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args();out=args.out.resolve()
    if (out/'summary.json').exists():parser.error('Use a fresh validation directory')
    out.mkdir(parents=True,exist_ok=True);checks=[];policy_path=ROOT/'config/catalytic_confirmation.json';before=sha(policy_path)
    refs=json.loads((ROOT/'references/processed/references.json').read_text());units=json.loads((ROOT/'references/validation/round09_identity_corpus/units.json').read_text())
    def check(name,passed,detail):checks.append(dict(test=name,passed=bool(passed),observed=json.dumps(detail)))
    def run(name,genes,mode='genome'):
        d=out/name;d.mkdir();source=d/'input_genes.json';source.write_text(json.dumps(genes,indent=2)+'\n')
        result=annotate(copy.deepcopy(genes),d,mode=mode,threads=1,input_files=[source])
        return result,table(d/'catalytic_confirmation.tsv')
    ref=next(r for r in refs if r['reference_id']=='unresolved__cautoethanogenum_DSM10061_CAETHG_1576')
    gene=next(g for g in extract_genbank(ROOT/'references/raw'/(ref['sample']+'.gb'))[0] if g['id']==ref['locus_tag'])
    result,confirmation=run('known_identity_incomplete_context',[gene]);a=next(r for r in result if r['protein_id']==gene['id'])
    check('known_source_identity_does_not_bypass_insufficient_context',a['catalytic_identity_confidence']=='not_high_confidence' and a['active_center_status']=='insufficient_alignment' and not a['potential_function_ids'],a)
    ref=next(r for r in refs if r['reference_id']=='unresolved__cautoethanogenum_DSM10061_CAETHG_3569')
    gene=next(g for g in extract_genbank(ROOT/'references/raw'/(ref['sample']+'.gb'))[0] if g['id']==ref['locus_tag'])
    result,confirmation=run('known_identity_unknown_function',[gene]);a=next(r for r in result if r['protein_id']==gene['id'])
    check('exact_exposed_identity_is_supported_without_function',a['catalytic_identity_confidence']=='calibrated_high' and a['catalytic_status']=='supported' and not a['potential_function_ids'],a)
    check('original_reference_kind_is_still_unresolved',next(r for r in json.loads((ROOT/'references/processed/references.json').read_text()) if r['reference_id']==ref['reference_id'])['kind']=='unresolved',ref['reference_id'])
    for accession in ['Q6DHP6','WP_028130118.1']:
        unit=next(u for u in units if any(s['accession']==accession for s in u['sources']))
        g=dict(id='reviewed_'+accession,sequence=unit['sequence'],partial=False,gene='',aliases='',protein_accession='',contig='unlinked',start='',end='',strand='',coordinate_evidence='translated_source_protein')
        result,rows=run(accession,[g],'transcript');r=rows[0]
        check(accession+':conserved_ligands_do_not_prove_identity',r['identity_confidence']=='not_high_confidence' and all(not a['potential_function_ids'] for a in result),r)
        check(accession+':negative_source_preserved',unit['label']=='negative',unit['sources'])
    ref=next(r for r in refs if r['reference_id']=='NiFe_Ecoli_Hyd1__HyaB')
    gene=next(g for g in extract_genbank(ROOT/'references/raw'/(ref['sample']+'.gb'))[0] if g['id']==ref['locus_tag'])
    truncated=dict(gene,sequence=gene['sequence'][:240],partial=True,coordinate_evidence='synthetic_source_truncation')
    result,rows=run('partial_catalytic_protein',[truncated])
    check('partial_CDS_cannot_be_high_identity',rows[0]['identity_confidence']=='not_high_confidence' and all(not a['potential_function_ids'] for a in result),rows)
    changed=list(gene['sequence']);changed[75]='X';ambiguous=dict(gene,sequence=''.join(changed),coordinate_evidence='synthetic_unknown_ligand')
    result,rows=run('unknown_ligand',[ambiguous])
    check('unknown_ligand_cannot_be_high_identity',rows[0]['identity_confidence']=='not_high_confidence' and all(not a['potential_function_ids'] for a in result),rows)
    check('policy_unchanged',before==sha(policy_path),before)
    write_tsv(out/'checks.tsv',checks)
    summary=dict(status='development_checks_passed' if all(c['passed'] for c in checks) else 'failed',checks=len(checks),passed=sum(c['passed'] for c in checks),
        failed=[c for c in checks if not c['passed']],policy_sha256=before,unchanged=before==sha(policy_path),
        scope='Six actual source-protein/input-state runs; covered and context-limited exposed identities, two development negatives and two synthetic boundaries. Not independent test accuracy or per-protein activity truth.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='failed'},indent=2))
    if summary['failed']:raise SystemExit(1)


if __name__=='__main__':main()
