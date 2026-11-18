"""Exercise raw Hnd preparation and claim-specific context/relay boundaries."""
import csv,hashlib,json,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'research_workflow')]
from pyridine_function import prepare_extension,load_model,annotation_evidence,predict
from annotate import BLAST_FIELDS
SID='FeFe_Sfructo_HndABCD'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
class HndTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=load_model();cls.refs=[r for r in json.loads((ROOT/'references/processed/complex_extensions/references.json').read_text()) if r['structure_id']==SID]
    def fixture(self,path,omit=(),partial=(),duplicate_b=False):
        ext=path/'structure_extensions';ext.mkdir(); genes=[];fa=[];blast=[];domains=[]
        ridmap={r['role']:r for r in self.refs}
        for i,role in enumerate(['HndA','HndB','HndC','HndD'],1):
            if role in omit:continue
            ref=ridmap[role];seq=ref['sequence'];pid=f'q{i:04d}';fa.append(f'>{pid}\n{seq}\n')
            genes.append(dict(id=pid,contig='ctg1',start=ref["start"],end=ref["end"],strand=1,bin_id='b1',partial=role in partial))
            h=dict(qseqid=pid,sseqid=ref['reference_id'],pident=100,length=len(seq),qstart=1,qend=len(seq),sstart=1,send=len(seq),evalue=0,bitscore=len(seq)*2,qlen=len(seq),slen=len(seq),qseq=seq,sseq=seq);blast.append(h)
            specs={'HndA':[('PF01257',30,171)],'HndC':[('PF01512',53,223),('PF10589',335,418)],'HndD':[('PF13510',3,83),('PF10588',90,128),('PF02906',229,518)]}.get(role,[])
            domains.extend(dict(protein_id=pid,accession=a,query_start=lo,query_end=hi,passes='True') for a,lo,hi in specs)
        if duplicate_b:
            ref=ridmap['HndB'];seq=ref['sequence'];fa.append(f'>duplicate\n{seq}\n');genes.append(dict(id='duplicate',contig='ctg1',start=950,end=1327,strand=1,bin_id='b1',partial=False));blast.append({**next(h for h in blast if h['sseqid'].endswith('__HndB')),'qseqid':'duplicate'})
        (path/'proteins.faa').write_text(''.join(fa));(ext/'queries.faa').write_text(''.join(fa));(path/'blast.tsv').write_text('')
        (ext/'blast.tsv').write_text(''.join('\t'.join(str(h[k]) for k in BLAST_FIELDS)+'\n' for h in blast))
        with (path/'genes.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(genes[0]),delimiter='\t');w.writeheader();w.writerows(genes)
        th=json.loads((ROOT/'config/thresholds.json').read_text());(path/'run_manifest.json').write_text(json.dumps(dict(thresholds=th,input_mode='genome')))
        m=dict(status='post_identity_structure_extension_search_complete',model_manifest_sha256=sha(ROOT/'references/processed/complex_extensions/manifest.json'),reaction_registry_sha256=sha(ROOT/'config/reaction_contracts.json'),script_sha256=sha(ROOT/'scripts/complex_extensions.py'),mode='genome',eligible_query_proteins=len(genes),source_proteins=len(genes),hits=len(blast),query_fasta_sha256=sha(ext/'queries.faa'),raw_alignment_sha256=sha(ext/'blast.tsv'))
        (ext/'search_manifest.json').write_text(json.dumps(m))
        return dict(protein_id='q0004',structure_id=SID,input_mode='genome',catalytic_identity_confidence='calibrated_high'),domains
    def evaluate(self,**kwargs):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);a,d=self.fixture(p,**kwargs);return prepare_extension(p,a,d,self.model,ROOT,SID)[0]
    def test_native_full_roles_support_coupling(self):
        r=self.evaluate();self.assertEqual(r['coupling_potential_status'],'local_coupling_potential_uncalibrated')
    def test_missing_B_retains_core_and_module(self):
        r=self.evaluate(omit=['HndB']);self.assertEqual(r['hydrogen_reaction_core_status'],'sequence_supported_candidate');self.assertEqual(r['pyridine_nucleotide_module_status'],'sequence_supported_candidate');self.assertEqual(r['coupling_potential_status'],'unresolved');self.assertIn('HndB_reference_context_unresolved',r['coupling_blockers'])
    def test_partial_B_only_blocks_coupling(self):
        r=self.evaluate(partial=['HndB']);self.assertEqual(r['hydrogen_reaction_core_status'],'sequence_supported_candidate');self.assertEqual(r['pyridine_nucleotide_module_status'],'sequence_supported_candidate');self.assertEqual(r['coupling_potential_status'],'unresolved')
    def test_ambiguous_B_only_blocks_coupling(self):
        r=self.evaluate(duplicate_b=True);self.assertEqual(r['hydrogen_reaction_core_status'],'sequence_supported_candidate');self.assertEqual(r['pyridine_nucleotide_module_status'],'sequence_supported_candidate');self.assertEqual(r['coupling_potential_status'],'unresolved')
    def test_partial_C_blocks_module_not_core(self):
        r=self.evaluate(partial=['HndC']);self.assertEqual(r['hydrogen_reaction_core_status'],'sequence_supported_candidate');self.assertEqual(r['pyridine_nucleotide_module_status'],'unresolved')
    def test_relay_loss_retains_H2_core(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);a,d=self.fixture(p);d=[r for r in d if r['accession']!='PF10588'];r=prepare_extension(p,a,d,self.model,ROOT,SID)[0];self.assertEqual(r['hydrogen_reaction_core_status'],'sequence_supported_candidate');self.assertEqual(r['coupling_potential_status'],'unresolved')
    def test_raw_hash_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);a,d=self.fixture(p);raw=p/'structure_extensions/blast.tsv';raw.write_text(raw.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'provenance'):prepare_extension(p,a,d,self.model,ROOT,SID)
    def test_bound_query_residue_mismatch_rejected_even_if_hash_rebound(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);a,d=self.fixture(p);raw=p/'structure_extensions/blast.tsv';raw.write_text(raw.read_text().replace('MQNSTCQ','AQNSTCQ'));m=p/'structure_extensions/search_manifest.json';s=json.loads(m.read_text());s['raw_alignment_sha256']=sha(raw);m.write_text(json.dumps(s))
            with self.assertRaisesRegex(ValueError,'sequence'):prepare_extension(p,a,d,self.model,ROOT,SID)
if __name__=='__main__':unittest.main()
