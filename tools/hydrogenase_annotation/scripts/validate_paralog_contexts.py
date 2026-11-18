#!/usr/bin/env python3
"""Primary-structure correspondences and full native-context paralog deployment checks."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import SimpleLocation
from annotate import ROOT,annotate,extract_genbank,write_tsv
from active_centers import structure_atoms
from catalytic_confirmation import sha,table
from run_acceptance_checks import state_hashes


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    raw=ROOT/'references/raw/round18_paralog_sources';before=state_hashes();checks=[];sources=[];runs=[]
    for record in json.loads((raw/'download_manifest.json').read_text()):
        if sha(raw/record['file'])!=record['sha256']:raise ValueError('Paralog source changed: '+record['file'])
        sources.append(record)
    def check(name,passed,observed):
        checks.append(dict(test=name,passed=bool(passed),observed=json.dumps(observed)))
    def run(name,genes,source,mode='genome'):
        anonymous=copy.deepcopy(genes);mapping={}
        for i,g in enumerate(anonymous):
            mapping[g['id']]='query_'+str(i);g.update(id=mapping[g['id']],gene='',aliases='',protein_accession='')
        directory=out/name;directory.mkdir();(directory/'source_correspondence.json').write_text(json.dumps(mapping,indent=2)+'\n')
        result=annotate(anonymous,directory,threads=2,mode=mode,input_files=[source])
        manifest=json.loads((directory/'run_manifest.json').read_text())
        runs.append(dict(name=name,genes=len(genes),mode=mode,input_sha256=sha(source),run_manifest_sha256=sha(directory/'run_manifest.json'),
            source_correspondence_sha256=sha(directory/'source_correspondence.json'),candidates=len(result),
            high_identity=manifest['catalytic_confirmation']['high_identity_outputs']))
        return {r['protein_id']:r for r in result},mapping

    pfur_path=ROOT/'references/raw/pfuriosus_DSM3638.gb';pfur=extract_genbank(pfur_path)[0]
    seqs=structure_atoms(raw/'6U8Y.cif')[0];mbs=[]
    for entity,sequence in seqs.items():
        matches=[g for g in pfur if g['sequence']==sequence]
        check('MBS_6U8Y_entity_'+entity+'_exact_native_CDS',len(matches)==1,[(g['id'],g['protein_accession']) for g in matches])
        if len(matches)!=1:raise ValueError('MBS native correspondence is not unique')
        g=matches[0];mbs.append(dict(entity=entity,locus=g['id'],old_locus=g['aliases'],accession=g['protein_accession'],start=g['start'],end=g['end'],strand=g['strand'],sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest()))
    write_tsv(out/'MBS_6U8Y_native_members.tsv',mbs)
    check('MBS_13_distinct_structural_members',len(mbs)==13 and len({r['locus'] for r in mbs})==13,len(mbs))
    found,mapping=run('Pfur_complete_genome',pfur,pfur_path)
    mbs_l=found.get(mapping['PF_RS07225'],{})
    check('native_MBS_L_has_no_hydrogenase_function',not mbs_l.get('potential_function_ids') and mbs_l.get('catalytic_identity_confidence')!='calibrated_high',mbs_l or 'not retained as catalytic anchor')
    mbh_l=found.get(mapping['PF_RS07185'],{})
    check('same_genome_MBH_positive_control_is_retained',mbh_l.get('catalytic_family')=='NiFe' and mbh_l.get('active_center_status')=='template_ligand_set_supported',mbh_l)

    dth_path=raw/'dthermolithotrophum_DSM11699.gb';dth=extract_genbank(dth_path)[0]
    dth_gene=next(g for g in dth if g['id']=='Dester_1504');native=dth_gene['sequence'];construct=structure_atoms(raw/'6HUZ.cif')[0]['1']
    check('DTH_active_HmdII_6HUZ_exact_native_sequence',len(native)==357 and construct[:357]==native and len(construct)==370,dict(native_length=len(native),polymer_length=len(construct),native_range=[1,357],terminal_construct_extension=construct[357:]))
    found,mapping=run('DTH_complete_genome',dth,dth_path);a=found.get(mapping['Dester_1504'],{})
    check('active_DTH_HmdII_retained_as_Fe_candidate',a.get('catalytic_family')=='Fe' and a.get('catalytic_status')=='candidate',a)
    check('THF_activity_does_not_transfer_H4MPT_or_high_identity',bool(a) and not a.get('potential_function_ids') and a.get('catalytic_identity_confidence')=='not_high_confidence',a)

    mj_path=raw/'mjannaschii_DSM2661.gb';mj=extract_genbank(mj_path)[0];gene=next(g for g in mj if 'MJ1338' in g['aliases'])
    construct=structure_atoms(raw/'4YT4.cif')[0]['1'];native=construct[:353]
    record=next(r for r in SeqIO.parse(mj_path,'genbank') if r.id==gene['contig'])
    extended=SimpleLocation(gene['start']-1-33,gene['end'],strand=1).extract(record.seq).translate(table=11,to_stop=True)
    check('MJ1338_alternative_start_explains_sequence_difference',gene['strand']==1 and len(gene['sequence'])==342 and native[11:]==gene['sequence'] and str(extended)==native,
        dict(current_accession=gene['protein_accession'],current_length=len(gene['sequence']),experimental_native_length=353,experimental_native_range=[1,353],polymer_length=len(construct),native_extension=11,extended_start=gene['start']-33,source_end=gene['end'],scope='Observed source correspondence; current annotation retained, no universal start-site correction'))
    found,mapping=run('MJ_complete_genome',mj,mj_path);a=found.get(mapping[gene['id']],{})
    check('MJ1338_genomic_annotation_has_no_transferred_Hmd_function',not a.get('potential_function_ids') and a.get('catalytic_identity_confidence')!='calibrated_high',a or 'not retained')
    native_gene=dict(gene,id='experimental_MJ1338',sequence=native,partial=False,start='',end='',strand='',contig='unlinked',coordinate_evidence='experimental_native_protein')
    found,mapping=run('MJ1338_experimental_native',[native_gene],raw/'4YT4.cif',mode='transcript');a=found.get(mapping['experimental_MJ1338'],{})
    check('FeGP_bound_MJ1338_does_not_prove_Hmd_reaction',not a.get('potential_function_ids') and a.get('catalytic_identity_confidence')!='calibrated_high',a or 'not retained')

    after=state_hashes();check('model_and_source_unchanged_during_runs',before==after,dict(changed=[k for k in set(before)|set(after) if before.get(k)!=after.get(k)]))
    write_tsv(out/'checks.tsv',checks)
    report=dict(status='development_context_checks_passed' if all(c['passed'] for c in checks) else 'failed',checks=len(checks),passed=sum(c['passed'] for c in checks),failed=[c for c in checks if not c['passed']],
        source_files=sources,additional_source_sha256={str(pfur_path):sha(pfur_path)},runs=runs,before_hashes=before,after_hashes=after,
        evidence_urls=['https://pubmed.ncbi.nlm.nih.gov/30181217/','https://doi.org/10.1038/s41467-020-19697-7','https://pubmed.ncbi.nlm.nih.gov/30600878/','https://pubmed.ncbi.nlm.nih.gov/26094576/'],
        scope='Three complete native genomes plus one experimental native protein; MBS is a known reference/training exposure, Hmd cases have been inspected for development. No independent accuracy, universal HmdII inactivity or measured query activity claim.')
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(status=report['status'],checks=len(checks),passed=report['passed'],failed=report['failed']),indent=2))
    if report['failed']:raise SystemExit(1)


if __name__=='__main__':main()
