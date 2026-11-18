#!/usr/bin/env python3
"""Extract genome-local maturation candidates with existing independent profile rules."""
import csv
import hashlib
import json
import re
import subprocess
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from audit_maturation_domains import scan_profile_domains, architecture_supported
from audit_reference_domains import model_evidence
from build_evidence_six import P, CURATION, OUT, WORK
import reorganize_evidence as catalog


def main():
    data=json.loads((WORK/'catalog.json').read_text())
    sources=defaultdict(set)
    for s in data['sheets']['systems_genomes.tsv']:
        if s['蛋白数']:
            rel=str(__import__('pathlib').Path(s['本地GenBank文件']).relative_to(P))
            sources[rel].add(s['活性中心'])
    for s in json.loads((CURATION/'added_systems.json').read_text()):sources['references/raw/'+s['source']].add(s['center'])
    directory=OUT/'genomic_helper_search';directory.mkdir(parents=True,exist_ok=True)
    queries={};owners=defaultdict(list);hcgs=[]
    for source,centers in sorted(sources.items()):
        for rec in catalog.genbank(source):
            for g in rec['genes']:
                digest=hashlib.sha256(g['sequence'].encode()).hexdigest()
                queries[digest]=g['sequence'];owners[digest].append((source,g,centers))
                if 'Fe' in centers and re.search(r'FeGP|cofactor biosynthesis protein Hmd',g['product'],re.I):
                    match=re.search(r'Hcg[A-G]|Hmd[BC]',g['gene']+' '+g['product'],re.I)
                    hcgs.append((source,g,match.group() if match else 'FeGP_biosynthesis_candidate'))
    query=directory/'genomic_proteins_U_to_X.faa'
    SeqIO.write([SeqRecord(Seq(s.replace('U','X')),id=qid,description='genomewide_helper_search') for qid,s in queries.items()],query,'fasta')
    reference=json.loads((P/'references/processed/maturation_domain_audit/summary.json').read_text())
    policy=json.loads((P/'config/maturation_domain_policy.json').read_text())
    assert hashlib.sha256((P/'references/processed/maturation_domain_audit/maturation_domains.hmm').read_bytes()).hexdigest()==reference['combined_model_sha256']
    domains=scan_profile_domains(P/'references/processed/maturation_domain_audit/maturation_domains.hmm',query,reference['profile_models'],policy['min_model_coverage'],directory,2)
    refs=json.loads((P/'references/processed/references.json').read_text())
    names={r['reference_id']:r['role'].split('_')[-1].lower() for r in refs if r['kind']=='maturation'}
    source_records={r.id:r for r in SeqIO.parse(P/'references/processed/reference_source_proteins.faa','fasta')}
    ref_faa=directory/'curated_maturation_references.faa'
    SeqIO.write([source_records[rid] for rid in names],ref_faa,'fasta')
    candidate_faa=directory/'domain_candidates.faa'
    SeqIO.write([SeqRecord(Seq(queries[qid]),id=qid,description='profile_candidate') for qid in domains],candidate_faa,'fasta')
    align=directory/'curated_maturation_matches.tsv'
    subprocess.run(['blastp','-query',str(candidate_faa),'-subject',str(ref_faa),'-evalue','1e-8','-max_hsps','1',
        '-outfmt','6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue','-out',str(align)],check=True)
    homologs=defaultdict(set)
    for line in align.read_text().splitlines():
        qid,rid,pident,length,qlen,slen,qs,qe,ss,se,evalue=line.split()
        if float(pident)>=25 and (int(qe)-int(qs)+1)/int(qlen)>=.65 and (int(se)-int(ss)+1)/int(slen)>=.65:
            homologs[qid].add(names[rid])
    chain_policy=json.loads((P/'config/reference_domain_policy.json').read_text())['chain_policy']
    rows=[];mapping=defaultdict(list)
    for qid,hits in domains.items():
        folded=model_evidence(hits,{'min_model_coverage':policy['min_model_coverage'],'chain_policy':chain_policy})
        accepted={d['accession'].split('.')[0] for d in folded if d['supported']}
        for role,rule in policy['rules'].items():
            if not architecture_supported(accepted,rule):continue
            for source,g,centers in owners[qid]:
                if (role.lower() in {'hyde','hydf','hydg'} and 'FeFe' not in centers) or (role.lower() not in {'hyde','hydf','hydg'} and 'NiFe' not in centers):continue
                named=any(a.lower()==role.lower() for a in g['aliases']) or re.search(r'\b'+re.escape(role)+r'\b',g['product'],re.I)
                if role.lower() not in homologs[qid] and not named:continue
                mapping[source].append([role,g['locus'],'maturation'])
                rows.append({'source':source,'role':role,'locus_tag':g['locus'],'protein_accession':g['protein_accession'],
                    'sequence_sha256':qid,'domains':';'.join(sorted(accepted)),'scope':'genomewide_profile_and_curated_homology_or_named_candidate; shared maturation; operon assignment unproven'})
    for source,g,role in hcgs:
        mapping[source].append([role,g['locus'],'maturation'])
        rows.append({'source':source,'role':role,'locus_tag':g['locus'],'protein_accession':g['protein_accession'],
            'sequence_sha256':hashlib.sha256(g['sequence'].encode()).hexdigest(),'domains':'GenBank FeGP biosynthesis annotation',
            'scope':'genomewide_annotation_candidate; not an experimentally proven Hmd-specific interaction'})
    for source in mapping:mapping[source]=[list(x) for x in sorted({tuple(r) for r in mapping[source]})]
    (CURATION/'genomic_helpers.json').write_text(json.dumps(mapping,ensure_ascii=False,indent=2)+'\n')
    with (OUT/'genomic_helper_candidates.tsv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    print(json.dumps({'unique_genomic_proteins':len(queries),'source_files':len(sources),'helper_candidates':len(rows),'sources_with_helpers':len(mapping)},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
