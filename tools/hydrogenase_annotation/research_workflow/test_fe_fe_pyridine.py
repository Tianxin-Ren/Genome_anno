"""Revision candidate: exercise FeFe normalization from raw members/genes/domains/hits."""
import copy, json, unittest
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'research_workflow')]
from pyridine_function import annotation_evidence, load_model, predict

class FeFeNormalizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.model=load_model()

    def fixture(self, structure='FeFe_TmHydABC', members=None, gene_updates=None, domain_updates=None, raw_extra=()):
        roles={'FeFe_TmHydABC':('HydA','HydB','HydC'),'FeFe_SwHyd1ABC':('Hyd1A','Hyd1B','Hyd1C')}[structure]
        anchor={'protein_id':'hyda','structure_id':structure,'input_mode':'genome','bin_id':'bin1'}
        genes=[{'id':'hyda','contig':'ctg','start':1000,'end':1300,'strand':1,'bin_id':'bin1','partial':False,'gene':'anchor','product':'HydA'}]
        mem=[]; domains=[]; raw=[]
        for i,role in enumerate(roles):
            pid='hyda' if i==0 else role.lower(); genes.append({'id':pid,'contig':'ctg','start':1500+i*600,'end':1800+i*600,'strand':1,'bin_id':'bin1','partial':False,'gene':'label_'+role,'product':'product_'+role}) if pid!='hyda' else None
            mem.append({'anchor':'hyda','structure_id':structure,'role':role,'protein_id':pid,'identity':55,'query_coverage':.9,'reference_coverage':.9,'query_domain_start':1,'query_domain_end':700,'cross_anchor_pairing_status':'no_overlapping_fixed_member_claim_detected'})
            raw.append({'qseqid':pid,'sseqid':structure+'__'+role,'pident':55.,'qcov':.9,'scov':.9,'evalue':1e-50,'bitscore':100.,'qlen':700,'qstart':1,'qend':700})
            if role in ('HydB','Hyd1B'):
                domains += [{'protein_id':pid,'accession':a,'query_start':lo,'query_end':hi,'passes':'True'} for a,lo,hi in (('PF01512',10,180),('PF10589',220,400))]
            elif role in ('HydC','Hyd1C'): domains.append({'protein_id':pid,'accession':'PF01257','query_start':10,'query_end':180,'passes':'True'})
            else: domains += [{'protein_id':pid,'accession':a,'query_start':lo,'query_end':hi,'passes':'True'} for a,lo,hi in (('PF13510',10,70),('PF10588',80,140),('PF22117',150,210),('PF02906',250,500))]
        if members is not None: mem=members
        for row in genes:
            if row['id'] in (gene_updates or {}): row.update(gene_updates[row['id']])
        for row in domains:
            row.update((domain_updates or {}).get(row['protein_id'],{}))
        raw.extend(raw_extra)
        return anchor,mem,genes,domains,raw

    def evaluate(self, identity='calibrated_high', **kw):
        a,m,g,d,r=self.fixture(**kw); e=annotation_evidence(a,m,g,identity,domains=d,model=self.model,raw_hits=r); return e,predict(e,self.model)

    def test_both_full_architectures(self):
        for sid in ('FeFe_TmHydABC','FeFe_SwHyd1ABC'):
            _,o=self.evaluate(structure=sid, identity='sequence_supported_uncalibrated_noncanonical' if sid.startswith('FeFe_') else 'calibrated_high'); self.assertEqual(o['coupling_potential_status'],'local_coupling_potential_uncalibrated')
        _,o=self.evaluate(identity='evidence_insufficient'); self.assertEqual(o['hydrogen_reaction_core_status'],'unresolved')

    def test_missing_roles_recomputed_from_members(self):
        a,m,g,d,r=self.fixture(); m=[x for x in m if x['role']!='HydB']; e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r); o=predict(e,self.model); self.assertEqual(o['hydrogen_reaction_core_status'],'sequence_supported_candidate'); self.assertEqual(o['coupling_potential_status'],'unresolved')
        a,m,g,d,r=self.fixture(); m=[x for x in m if x['role']!='HydA']; e=annotation_evidence(a,m,g,'evidence_insufficient',domains=d,model=self.model,raw_hits=r); o=predict(e,self.model); self.assertEqual(o['hydrogen_reaction_core_status'],'unresolved'); self.assertEqual(o['pyridine_nucleotide_module_status'],'sequence_supported_candidate')

    def test_relay_topology_is_raw_domain_evidence(self):
        a,m,g,d,r=self.fixture(); d=[x for x in d if x['accession']!='PF10588']; e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r); o=predict(e,self.model); self.assertEqual(o['hydrogen_reaction_core_status'],'sequence_supported_candidate'); self.assertEqual(o['coupling_potential_status'],'unresolved')
        a,m,g,d,r=self.fixture(); d=[dict(x,query_start=600,query_end=650) if x['accession']=='PF13510' else x for x in d]; e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r); self.assertEqual(predict(e,self.model)['coupling_potential_status'],'unresolved')

    def test_raw_rivals_and_resource_scope(self):
        rival={'qseqid':'hydb','sseqid':'NiFe_MvhHdr__MvhA','pident':50.,'qcov':.9,'scov':.9,'evalue':1e-40,'bitscore':95.,'qlen':700,'qstart':1,'qend':700}
        a,m,g,d,r=self.fixture(raw_extra=[rival]); e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r); self.assertEqual(predict(e,self.model)['coupling_potential_status'],'unresolved')
        a,m,g,d,r=self.fixture(gene_updates={'hydb':{'bin_id':'bin2'}}); a['input_mode']='bin'; e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r); self.assertEqual(predict(e,self.model)['coupling_potential_status'],'unresolved')
        a,m,g,d,r=self.fixture(domain_updates={'hydb':{'query_start':600,'query_end':650}}); e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r); self.assertEqual(predict(e,self.model)['pyridine_nucleotide_module_status'],'unresolved')

    def test_compatible_sw_rival_and_unknown_labels_do_not_change_validity(self):
        rival={'qseqid':'hydb','sseqid':'FeFe_SwHyd1ABC__Hyd1B','pident':50.,'qcov':.9,'scov':.9,'evalue':1e-40,'bitscore':95.,'qlen':700,'qstart':1,'qend':700}
        a,m,g,d,r=self.fixture(raw_extra=[rival]);
        for gene in g: gene['gene']='changed_label'; gene['product']='changed_product'
        e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r); self.assertEqual(predict(e,self.model)['coupling_potential_status'],'local_coupling_potential_uncalibrated')

    def test_hyt_aw_short_hyd_a_lacks_internal_relay(self):
        a,m,g,d,r=self.fixture(); d=[x for x in d if x['protein_id']!='hyda' or x['accession']=='PF02906']
        e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r); o=predict(e,self.model); self.assertEqual(o['hydrogen_reaction_core_status'],'sequence_supported_candidate'); self.assertEqual(o['pyridine_nucleotide_module_status'],'sequence_supported_candidate'); self.assertEqual(o['coupling_potential_status'],'unresolved')

    def test_raw_resource_reuse_blocks_coupling(self):
        a,m,g,d,r=self.fixture(); c=next(x for x in m if x['role']=='HydC'); c['protein_id']='hydb'; d.append({'protein_id':'hydb','accession':'PF01257','query_start':450,'query_end':600,'passes':'True'}); r.append({'qseqid':'hydb','sseqid':'FeFe_TmHydABC__HydC','pident':55.,'qcov':.9,'scov':.9,'evalue':1e-50,'bitscore':90.,'qlen':700,'qstart':1,'qend':700}); e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r); self.assertEqual(predict(e,self.model)['coupling_potential_status'],'unresolved')

    def test_candidate_count_ignores_noncompetitive_unfeasible_paralog(self):
        a,m,g,d,r=self.fixture(); nodes=[{'anchor':'hyda','structure_id':'FeFe_TmHydABC','role':'HydB','candidate_count':'2','graph_search_status':'complete','candidate_realizations':json.dumps([{'protein_ids':['hydb'],'competitive':True,'jointly_feasible':True},{'protein_ids':['decoy'],'competitive':False,'jointly_feasible':False}])}]
        e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r,nodes=nodes); self.assertNotIn('multiple_role_partner_candidates_HydB',e['module_competition_blockers'])

    def test_candidate_count_blocks_two_valid_competitive_copies(self):
        a,m,g,d,r=self.fixture(); nodes=[{'anchor':'hyda','structure_id':'FeFe_TmHydABC','role':'HydB','candidate_count':'2','graph_search_status':'complete','candidate_realizations':json.dumps([{'protein_ids':['hydb'],'competitive':True,'jointly_feasible':True},{'protein_ids':['decoy'],'competitive':True,'jointly_feasible':True}])}]
        e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r,nodes=nodes); self.assertIn('multiple_role_partner_candidates_HydB',e['module_competition_blockers'])

    def test_truncated_or_missing_graph_status_remains_conservative(self):
        for status in ('truncated', None):
            a,m,g,d,r=self.fixture(); node={'anchor':'hyda','structure_id':'FeFe_TmHydABC','role':'HydB','candidate_count':'2','candidate_realizations':json.dumps([{'protein_ids':['hydb'],'competitive':True,'jointly_feasible':True},{'protein_ids':['decoy'],'competitive':False,'jointly_feasible':False}])}
            if status is not None: node['graph_search_status']=status
            e=annotation_evidence(a,m,g,'calibrated_high',domains=d,model=self.model,raw_hits=r,nodes=[node]); self.assertIn('multiple_role_partner_candidates_HydB',e['module_competition_blockers'])

if __name__=='__main__': unittest.main()
