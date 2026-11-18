#!/usr/bin/env python3
"""Build source-reviewed post-identity structure templates; never rebuild L1/identity models."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, extract_genbank
from catalytic_confirmation import sha
from evaluate_natural_complexes import verify_corpus


def source_case(spec,truth,sources,root=ROOT):
    if 'native_case' in spec:
        case=spec['native_case']
        for name,digest in case['source_files'].items():
            if sha(Path(root)/name)!=digest:raise ValueError('Reviewed native source changed: '+name)
        return case,Path(root)/case['input_file']
    case=truth[spec['case_id']]
    return case,Path(root)/sources[case['source']]['input_file']


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    review_path=ROOT/'config/complex_extension_sources.json';review=json.loads(review_path.read_text())
    corpus=ROOT/review['source_corpus'];freeze=json.loads((corpus/'freeze.json').read_text());verify_corpus(corpus,freeze)
    truth={c['id']:c for c in json.loads((corpus/'truth.json').read_text())};sources=json.loads((corpus/'sources.json').read_text())
    refs=[];systems=[];graphs=[];checks=[];source_files={}
    for spec in review['systems']:
        case,source=source_case(spec,truth,sources);genes={g['id']:g for g in extract_genbank(source)[0]}
        source_files[str(source.relative_to(ROOT))]=sha(source);source_files[case['paper']]=sha(ROOT/case['paper'])
        source_files.update(case.get('source_files',{}))
        sid=spec['structure_id'];members={m['role']:m for m in case['members']};nodes=[];loci={}
        if spec['center']!=('NiFe' if sid.startswith('NiFe_') else 'FeFe'):raise ValueError('Unreviewed catalytic center')
        if set(spec['role_categories'])!=set(members):raise ValueError('Unreviewed or missing physical member')
        for role,member in members.items():
            g=genes[member['source_locus']];sequence=g['sequence'];digest=hashlib.sha256(sequence.encode()).hexdigest()
            if digest!=member['sequence_sha256'] or g['protein_accession']!=member['protein_accession']:raise ValueError('Native source changed: '+g['id'])
            rid=sid+'__'+role;loci[role]=g['id'];catalytic=role==spec['catalytic_role']
            refs.append(dict(reference_id=rid,structure_id=sid,role=role,kind='catalytic' if catalytic else 'module',center=spec['center'],
                sample=case['source'],locus_tag=g['id'],aliases='',protein_accession=g['protein_accession'],contig=g['contig'],start=g['start'],end=g['end'],strand=g['strand'],
                sequence=sequence,sequence_sha256=digest,source_sequence_sha256=digest,source_protein_length=len(sequence),protein_region_start=1,protein_region_end=len(sequence),
                cysteine_positions=[i+1 for i,a in enumerate(sequence) if a=='C'],selenocysteine_positions=[i+1 for i,a in enumerate(sequence) if a=='U'],
                diagnostic_positions=[],diagnostic_residues={},diagnostic_scope='No new L2 seed model; independent pre-existing active-center/span evidence is required by the post-identity selector.',
                source_accession=g['contig'],source_url='https://www.ebi.ac.uk/ena/browser/view/'+g['contig'],protein_source_url='https://www.ncbi.nlm.nih.gov/protein/'+g['protein_accession'],
                evidence_url=case['url'],evidence_scope='reviewed_purified_reference_member_exact_native_CDS_correspondence',coverage_scope='full_reference_protein',
                domain_validation='experimental_member_correspondence_not_new_domain_classifier',domain_specificity_scope=spec['review']))
            nodes.append(dict(node_id=role,role=role,reference_id=rid,reference_sequence_sha256=digest,reference_start=1,reference_end=len(sequence),
                segment_scope='full_experimentally_corresponding_subunit',module_category=spec['role_categories'][role],required_for_structure=True,
                required_for_function=[],linkage='local',source=case['url'],source_scope=spec['review']))
            checks.append(dict(structure_id=sid,role=role,locus=g['id'],protein=g['protein_accession'],uniprot=member.get('uniprot',''),length=len(sequence),sequence_sha256=digest,source_verified=True))
        systems.append(dict(structure_id=sid,sample=case['source'],center=spec['center'],reference_clade=spec['reference_clade'],catalytic_role=spec['catalytic_role'],local=loci,remote={},context={},
            function_ids=[],reference_function_records=[],reference_contract_ids=spec['contract_ids'],direct_partner=spec['direct_partner'],mechanism=spec['mechanism'],
            evidence_url=case['url'],reference_small_subunit_topology='not_assigned_by_structure_extension; independent_FeS_mapping_required',reference_energy_coupling='unresolved'))
        graphs.append(dict(structure_id=sid,nodes=nodes))
    if (len(refs)!=sum(len(s['role_categories']) for s in review['systems']) or len(systems)!=len(review['systems']) or
        sum(r['kind']=='catalytic' for r in refs)!=len(systems)):raise ValueError('Unexpected reviewed bank composition')
    for name,data in [('references.json',refs),('systems.json',systems),('graphs.json',graphs),('source_checks.json',checks)]:
        (out/name).write_text(json.dumps(data,indent=2)+'\n')
    SeqIO.write([SeqRecord(Seq(r['sequence']),id=r['reference_id'],description='') for r in refs],out/'references.faa','fasta')
    with (out/'makeblastdb.log').open('w') as log:
        subprocess.run(['makeblastdb','-in',str(out/'references.faa'),'-dbtype','prot','-out',str(out/'seed')],stdout=log,stderr=log,check=True)
    manifest=dict(status='reviewed_development_structure_extension_bank',scope=review['scope'],policy_sha256=sha(review_path),builder_sha256=sha(Path(__file__)),
        source_corpus_freeze_sha256=sha(corpus/'freeze.json'),source_files=source_files,reference_replay_cases=[s['case_id'] for s in review['systems']],
        structures=len(systems),proteins=len(refs),catalytic_templates=len(systems),legacy_identity_model_modified=False,
        artifacts={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({k:manifest[k] for k in ['status','structures','proteins','catalytic_templates','legacy_identity_model_modified']},indent=2))


if __name__=='__main__':main()
