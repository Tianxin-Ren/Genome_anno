import copy
import unittest

from function_rules import alignment_evidence,decide,donor_comparisons,evaluate_architecture,load_rules,overlaps


class FunctionRuleTests(unittest.TestCase):
    def setUp(self):
        self.cls=dict(id='Fe08',label='NADP reaction',unresolved_qualifiers=['Fd_nonrequirement'])
        self.arch=dict(id='test',system='test',members=[
            dict(role='HydA',reference='A',category='催化亚基',requirement='required'),
            dict(role='HydE',reference='E',category='电子传递或其他相关亚基',requirement='required'),
            dict(role='hydG',reference='G',category='成熟辅助因子',requirement='maturation')],
            association='reference_layout',support_mode='reference_architecture_prediction',citations=[],limitations='')
        self.model=dict(search=dict(max_local_span_bp=30000,catalytic_identity_margin=5,min_identity=35),references={
            x:dict(accession='chromosome',start=i*1000,end=i*1000+300) for i,x in enumerate('AEG')})
        self.genes={x:dict(id=x,contig='chromosome',start=i*1000+1,end=i*1000+300,strand=1,contig_length=10000,bin_id='bin1') for i,x in enumerate('AEG')}
        self.hits={x:[dict(protein_id=x,reference_id=x,identity=100,reference_coverage=1,bitscore=100,
                         query_interval=[1,100],reference_interval=[1,100],admissible=True,failures=[])] for x in 'AEG'}

    def evaluate(self,identity='sequence_supported_uncalibrated_core'):
        return evaluate_architecture(self.cls,self.arch,'A',self.hits,self.genes,identity,self.model,100,'bin')

    def test_compiled_34_classes_and_sequence_bank(self):
        model,_,_=load_rules()
        self.assertEqual(len(model['classes']),34)
        a=next(a for c in model['classes'] for a in c['architectures'] if a['system']=='FeFe_Mthermo_HydABCDEF')
        self.assertEqual({m['role'] for m in a['members'] if m['requirement']=='required'},set(['HydA','HydB','HydC','HydD','HydE','HydF']))
        fragment=model['references']['FRAG_BG04_163_BMR_473_1049']
        self.assertEqual(len(fragment['sequence']),577)
        self.assertEqual(fragment['construct_range'],[473,1049])

    def test_complete_members_predict_function(self):
        self.assertEqual(self.evaluate()['status'],'reference_architecture_supported')

    def test_missing_maturase_does_not_revoke_function(self):
        self.hits['G']=[]
        r=self.evaluate();self.assertEqual(r['status'],'reference_architecture_supported')
        self.assertEqual(r['maturation_candidates'],[])

    def test_deleting_complex_hydE_revokes_complete_architecture(self):
        self.hits['E']=[]
        r=self.evaluate();self.assertEqual(r['status'],'partial_architecture')
        self.assertEqual(r['missing_required'],['HydE'])
        self.assertIn('HydA',[m['role'] for m in r['members']])

    def test_cross_bin_member_cannot_complete_function(self):
        self.genes['E']['bin_id']='bin2'
        self.assertEqual(self.evaluate()['status'],'partial_architecture')

    def test_unconfirmed_identity_cannot_be_promoted_by_members(self):
        self.assertEqual(self.evaluate('evidence_insufficient')['status'],'excluded_identity_or_center')

    def test_ligand_conflict_blocks_function(self):
        self.hits['A'][0].update(admissible=False,failures=['incompatible_catalytic_ligand'])
        self.assertIn('incompatible_catalytic_ligand',self.evaluate()['blockers'])

    def test_overlapping_roles_not_counted_twice(self):
        self.hits['E'][0]['protein_id']='A'
        self.assertEqual(self.evaluate()['status'],'partial_architecture')

    def test_nonoverlapping_fusion_can_supply_two_roles(self):
        self.hits['E'][0].update(protein_id='A',query_interval=[101,200])
        self.assertEqual(self.evaluate()['status'],'reference_architecture_supported')

    def test_distant_members_in_local_template_reduce_claim(self):
        self.genes['E']['start']=70000
        self.assertIn('member_linkage_unresolved',self.evaluate()['blockers'])

    def test_weak_reference_does_not_prove_full_mechanism(self):
        self.arch['support_mode']='candidate_only'
        r=self.evaluate();self.assertEqual(r['status'],'candidate_only')
        self.assertEqual(r['function_confidence'],'low')

    def test_external_reconstitution_not_native_complex(self):
        self.arch['association']='external'
        r=self.evaluate();self.assertEqual(r['status'],'external_reconstitution_potential')

    def test_weak_reference_cannot_override_missing_members(self):
        self.arch['support_mode']='candidate_only';self.hits['E']=[]
        self.assertEqual(self.evaluate()['status'],'partial_architecture')

    def test_external_reference_cannot_override_missing_members(self):
        self.arch['association']='external';self.hits['E']=[]
        self.assertEqual(self.evaluate()['status'],'partial_architecture')

    def test_external_reference_cannot_promote_unconfirmed_core(self):
        self.arch['association']='external'
        self.assertEqual(self.evaluate('evidence_insufficient')['status'],'excluded_identity_or_center')

    def test_internal_deletion_not_hidden_by_alignment_span(self):
        ref=dict(sequence='ACDEFGHIKL',metal_sites=[[2,'C']],domains=[dict(accession='domain',segments=[[3,6]])])
        raw=dict(qseqid='q',sseqid='r',qseq='AC----HIKL',sseq='ACDEFGHIKL',qstart=1,qend=6,sstart=1,send=10,pident=100,bitscore=20)
        result=alignment_evidence(raw,ref,'ACHIKL',dict(min_reference_coverage=.8))
        self.assertEqual(result['reference_coverage'],.6)
        self.assertIn('reference_domain_span_incomplete',result['failures'])

    def test_ligand_mapping_uses_actual_sec_residue(self):
        ref=dict(sequence='ACDE',metal_sites=[[2,'U']],domains=[])
        raw=dict(qseqid='q',sseqid='r',qseq='AXDE',sseq='ACDE',qstart=1,qend=4,sstart=1,send=4,pident=75,bitscore=20)
        self.assertEqual(alignment_evidence(raw,ref,'AUDE',dict(min_reference_coverage=.8))['ligand_changes'],[])

    def test_competing_mechanisms_are_not_both_reported_as_supported(self):
        one=self.evaluate();one['function_id']='Fe03';one['score']=.90
        two=copy.deepcopy(one);two['function_id']='Fe05';two['score']=.89;two['system']='other_reference'
        summary=decide([one,two])
        self.assertEqual(summary['function_rule_status'],'ambiguous_architecture')
        self.assertEqual(len(summary['function_class_ids']),1)
        self.assertEqual(summary['function_rule_confidence'],'low')

    def test_reverse_directions_can_coexist(self):
        one=self.evaluate();one['function_id']='Fe03'
        two=copy.deepcopy(one);two['function_id']='Fe04'
        self.assertEqual(set(decide([one,two])['function_class_ids']),{'Fe03','Fe04'})

    def test_one_source_can_have_multiple_cofactor_capabilities(self):
        one=self.evaluate();one['function_id']='Ni12R'
        two=copy.deepcopy(one);two['function_id']='Ni14'
        self.assertEqual(set(decide([one,two])['function_class_ids']),{'Ni12R','Ni14'})

    def test_stronger_nonformate_partner_is_excluded(self):
        self.arch['members'][1]['role']='FdhF'
        self.model['role_counters']={'Fdh_catalytic':['negative']}
        self.hits['negative']=[dict(self.hits['E'][0],reference_id='negative',bitscore=300)]
        r=self.evaluate();self.assertEqual(r['status'],'partial_architecture')
        self.assertIn('stronger_nonformate_partner_reference:negative',r['member_rejections']['FdhF'])

    def test_hmdii_cannot_be_overwritten_by_h4mpt_rule(self):
        self.cls['id']='H01'
        r=evaluate_architecture(self.cls,self.arch,'A',self.hits,self.genes,'sequence_supported_uncalibrated_core',self.model,100,'bin',{'preferred_family':'HmdII'})
        self.assertEqual(r['status'],'excluded_incompatible_function')


class MvhModuleAssociationTests(unittest.TestCase):
    def setUp(self):
        self.roles={'MvhA':'core','MvhG':'small','HdrA':'flavin','MvhD':'relay','HdrB':'hdrB','HdrC':'hdrC'}
        locations={'core':(300001,301500,1),'small':(299000,300000,1),'flavin':(50000,52000,-1),
                   'relay':(49580,49996,-1),'hdrB':(55000,55800,-1),'hdrC':(55800,56350,-1)}
        self.genes={p:dict(id=p,contig='chromosome',start=a,end=b,strand=s,bin_id='one',contig_length=500000,topology='linear') for p,(a,b,s) in locations.items()}
        members=[dict(role=r,reference=r,category='催化亚基' if r=='MvhA' else '电子传递或其他相关亚基',requirement='required') for r in self.roles]
        self.arch=dict(id='test_mvh',system='test_mvh',members=members,association='reference_layout',support_mode='reference_architecture_prediction',citations=[],limitations='')
        self.cls=dict(id='Ni13',label='H2/Hdr reaction',unresolved_qualifiers=['obligatory_bifurcation'])
        self.model=dict(search=dict(max_local_span_bp=30000,catalytic_identity_margin=5,min_identity=35),
            references={r:dict(accession='source',start=i*1000,end=i*1000+500) for i,r in enumerate(self.roles)},
            module_prediction_labels={'Ni13':'Mvh-Hdr coupling prediction'},
            _donor_comparisons={'flavin':dict(donor_state='donor_unresolved',identity_margin=2.465,alternative_bitscore_ratio=579/617)})
        self.hits={r:[dict(protein_id=p,reference_id=r,identity=100,reference_coverage=1,bitscore=100,
                          query_interval=[1,650] if r=='HdrA' else [1,100],reference_interval=[1,100],admissible=True,failures=[])] for r,p in self.roles.items()}

    def evaluate(self,identity='calibrated_high'):
        return evaluate_architecture(self.cls,self.arch,'core',self.hits,self.genes,identity,self.model,100,'bin')

    def test_adjacent_separate_HdrA_MvhD_support_complete_function(self):
        row=self.evaluate()
        self.assertEqual(row['status'],'modular_function_prediction')
        self.assertEqual(row['function_confidence'],'medium')
        self.assertEqual(row['prediction_scope'],'genome_level_coupling_prediction')
        self.assertIn('HdrA_donor_specificity_not_resolved_by_two_reference_sequences',row['unmeasured_qualifiers'])
        self.assertTrue(row['modular_association']['HdrA_MvhD_local_pair'])

    def test_donor_comparison_is_owned_by_HdrA_not_separate_MvhD(self):
        self.model['_donor_comparisons']['relay']=dict(donor_state='alternative_favored',alternative_bitscore_ratio=2)
        self.assertEqual(self.evaluate()['status'],'modular_function_prediction')

    def test_fused_and_separate_modules_have_same_function(self):
        before=decide([self.evaluate()])['function_class_ids']
        self.hits['MvhD'][0].update(protein_id='flavin',query_interval=[700,840])
        self.assertEqual(decide([self.evaluate()])['function_class_ids'],before)

    def test_five_point_hint_does_not_veto_complete_modules(self):
        self.model['HdrA_donor_comparison']={'H2':'positive','F420_Fd_alternative':'alternative'}
        for margin in [-4.999,0,2.465,4.999,5.0]:
            with self.subTest(margin=margin):
                p=dict(protein_id='flavin',identity=50+margin,reference_coverage=1,bitscore=100,query_interval=[1,650])
                a=dict(p,identity=50)
                self.model['_donor_comparisons']=donor_comparisons({'positive':[p],'alternative':[a]},self.model)
                self.assertEqual(self.evaluate()['status'],'modular_function_prediction')

    def test_clear_alternative_donor_remains_excluded(self):
        self.model['_donor_comparisons']['flavin'].update(donor_state='alternative_favored',identity_margin=-15,alternative_bitscore_ratio=1.3)
        row=self.evaluate()
        self.assertIn('strong_competing_F420_Fd_HdrA_reference',row['blockers'])
        self.assertFalse(any('Near-tied' in risk for risk in row['risks']))

    def test_small_bitscore_difference_does_not_establish_alternative_donor(self):
        self.model['_donor_comparisons']['flavin'].update(donor_state='alternative_favored',identity_margin=-5.1,alternative_bitscore_ratio=1.01)
        self.assertEqual(self.evaluate()['status'],'modular_function_prediction')

    def test_unlinked_remote_relay_still_blocks_coupling(self):
        self.genes['relay'].update(start=40000,end=40400)
        self.assertIn('remote_MvhD_requires_linked_or_fused_HdrA',self.evaluate()['blockers'])

    def test_opposite_strand_does_not_supply_adjacent_module(self):
        self.genes['relay']['strand']=1
        self.assertIn('remote_MvhD_requires_linked_or_fused_HdrA',self.evaluate()['blockers'])

    def test_cross_bin_relay_cannot_complete_function(self):
        self.genes['relay']['bin_id']='different'
        self.assertIn('MvhD',self.evaluate()['missing_required'])

    def test_missing_HdrC_cannot_complete_function(self):
        self.hits['HdrC']=[]
        self.assertEqual(self.evaluate()['status'],'partial_architecture')

    def test_destroyed_HdrC_ligands_cannot_complete_function(self):
        self.hits['HdrC'][0]['coordination_evidence']={'HdrC':[dict(observed='S') for _ in range(8)]}
        self.assertIn('HdrC',self.evaluate()['missing_required'])

    def test_competing_hydrogenase_cores_keep_pairing_ambiguous(self):
        self.genes['second_core']=dict(self.genes['core'],id='second_core',start=330000,end=331500)
        self.hits['MvhA'].append(dict(self.hits['MvhA'][0],protein_id='second_core'))
        self.assertIn('remote_donor_pairing_ambiguous_between_cores',self.evaluate()['blockers'])

    def test_unconfirmed_core_is_not_promoted_by_complete_modules(self):
        self.assertEqual(self.evaluate('evidence_insufficient')['status'],'excluded_identity_or_center')


if __name__=='__main__':unittest.main()
