#!/usr/bin/env python3
"""Test real split domains and deliberately damaged sequences with unchanged source models."""
import json
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from active_centers import annotate_active_centers
from audit_maturation_domains import scan_profile_domains
from build_evidence_six import P, ROOT, WORK
from evidence_domain_fragments import summarize_domains, structure_correspondences
from query_domain_organization import load_models
from localization import sha


def main():
    data=json.loads((WORK/'catalog.json').read_text()); existing=json.loads((WORK/'mapped_features.json').read_text())
    seqs=data['sheets']['sequences.tsv'];out=ROOT/'hydrogenase_annotation/experiments/domain_fragments_20261112'
    out.mkdir(parents=True,exist_ok=True)
    policy,models,combined,resources=load_models(P)
    grouped=summarize_domains(existing['domains'],policy)
    old={d['序列ID'] for d in existing['domains'] if d['passes']}
    new={qid for qid,rs in grouped.items() if any(r['supported'] for r in rs)}
    def protein(sid,role):return next(r for r in seqs if sid==r['体系ID'] and r['亚基或基因名']==role)
    huc=protein('NiFe_HucLMS','HucL');ech=protein('NiFe_EchABCDEF','EchE')
    shi=protein('NiFe_HyhNADPH','SHI_alpha');frh=protein('NiFe_FrhABG','FrhA')
    examples=[]
    for r,acc,want in [(huc,'PF00374',428/508),(ech,'PF00346',238/271),(shi,'PF00374',279/508),(frh,'PF00374',169/508)]:
        matches=[d for d in grouped[r['序列ID']] if d['accession']==acc]
        assert len(matches)==1,(r['亚基或基因名'],acc)
        d=matches[0];assert abs(d['model_coverage']-want)<1e-8,(r['亚基或基因名'],d['model_coverage'],want)
        examples.append({'protein':r['亚基或基因名'],'sequence_id':r['序列ID'],**d})
    structures=structure_correspondences([shi,frh],P)
    assert all(r['序列ID'] in structures for r in [shi,frh])
    s=huc['完整CDS蛋白序列']
    cases={
        'HucL_native':s,
        'HucL_Nterminal_deletion':s[130:],
        'HucL_Cterminal_deletion':s[:-100],
        'HucL_internal_deletion':s[:200]+s[370:],
        'HucL_repeated_Nterminal_fragment':s[:140]+s[:140],
        'HucL_excessive_spacer':s[:140]+'G'*800+s[140:],
    }
    references=json.loads((P/'references/processed/references.json').read_text())
    source_sequences={r.id:str(r.seq) for r in SeqIO.parse(P/'references/processed/reference_source_proteins.faa','fasta')}
    decoys=[r for r in references if r['kind']=='decoy']
    for r in decoys:cases['noncatalytic_'+r['reference_id']]=source_sequences[r['reference_id']]
    faa=out/'challenges.faa'
    SeqIO.write([SeqRecord(Seq(s),id=qid,description='natural_negative_reference' if qid.startswith('noncatalytic_') else 'natural_positive_reference' if qid=='HucL_native' else 'synthetic_control_not_natural_genome') for qid,s in cases.items()],faa,'fasta')
    ds=scan_profile_domains(combined,faa,models,policy['min_model_coverage'],out/'domains',2)
    rows=[{'序列ID':qid,'模型集':'reference',**d} for qid,hits in ds.items() for d in hits]
    support=summarize_domains(rows,policy)
    center,_=annotate_active_centers(faa,out,threads=2,anchors=set(cases))
    bycenter={r['protein_id']:r for r in center}
    results=[]
    for qid in cases:
        pf=[d for d in support.get(qid,[]) if d['accession']=='PF00374']
        covered=any(d['supported'] for d in pf)
        any_fold=any(d['supported'] for d in support.get(qid,[]))
        c=bycenter.get(qid,{})
        intact=c.get('active_center_status')=='template_ligand_set_supported' and c.get('template_span_alignment_status')=='supported'
        results.append({'case':qid,'NiFeSe_Hases_profile_covered':covered,'any_panel_profile_covered':any_fold,'NiFeSe_Hases_model_coverage':max([d['model_coverage'] for d in pf],default=0),
                        'center_status':c.get('active_center_status','not_mapped'),
                        'template_span_status':c.get('template_span_alignment_status','not_mapped'),
                        'profile_and_center_and_span':any_fold and intact})
    assert results[0]['NiFeSe_Hases_profile_covered'],results
    assert results[0]['profile_and_center_and_span'],results
    assert not any(r['NiFeSe_Hases_profile_covered'] for r in results if r['case'] in {'HucL_repeated_Nterminal_fragment','HucL_excessive_spacer'}),results
    assert not any(r['profile_and_center_and_span'] for r in results[1:]),results
    assert any(r['any_panel_profile_covered'] for r in results if r['case'].startswith('noncatalytic_')),results
    assert all(sha(Path(path))==digest for path,digest in resources.items())
    report={'status':'passed','existing_min_model_coverage':policy['min_model_coverage'],'existing_chain_policy':policy['chain_policy'],
        'proteins_with_single_segment_support':len(old),'proteins_with_single_or_collinear_support':len(new),
        'rescued_from_single_segment_false_missing':len(new-old),'real_examples':examples,'challenge_results':results,
        'natural_negative_references':len(decoys),'synthetic_damage_controls':5,
        'source_resources_sha256':resources,'source_models_unchanged':True,
        'scope':'Domain coverage and site/span checks only; constructed sequence challenge is not an independent natural-genome accuracy benchmark.'}
    (out/'results.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in {'source_resources_sha256','real_examples'}},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
