import copy
import csv
import json
import hashlib
import tempfile
import unittest
from Bio import SeqIO
from pathlib import Path
from identity_scope import local_domain_evidence,sequence_core_support,load_lineage_review,complete_fefe_core_support,reference_backbones,nife_backbone_support,nife_core_regions,nife_core_region_support,NIFE_CORE_REGION_STATUS
from fefe_profile_ligands import prepare as prepare_fefe_profile
from coding_completeness import prepare as coding_prepare,load_evidence as load_coding


class IdentityScopeTests(unittest.TestCase):
    def setUp(self):
        self.f=dict(family='NiFe',eligibility_reasons=[],positive_prototype='p',positive_score=70,competitor_margin=20,
            core_domain_coverage=.33,source_active_center=dict(direct_ligand_status='present',active_center_status='template_ligand_set_supported',template_span_alignment_status='supported'),
            domain_evidence=[dict(accession='PF00374',profile='NiFeSe_Hases',model_coverage=.33,supported=False,query_core_overlap=.9,
                components=[dict(sequence_bitscore=100,sequence_ga=25,domain_bitscore=95,domain_ga=25)])])
        self.p=dict(min_domain_coverage=.8,min_core_similarity=30,min_competitor_margin=0)

    def test_local_family_support_does_not_require_complete_global_profile(self):
        before=copy.deepcopy(self.f)
        self.assertTrue(sequence_core_support(self.f,{},self.p));self.assertEqual(self.f,before)
        self.assertFalse(local_domain_evidence(self.f)[0]['full_domain_model_supported'])

    def test_GA_and_core_colocation_cannot_be_bypassed(self):
        self.f['domain_evidence'][0]['components'][0]['domain_bitscore']=24
        self.assertFalse(sequence_core_support(self.f,{},self.p))
        self.f['domain_evidence'][0]['components'][0]['domain_bitscore']=95
        self.f['domain_evidence'][0]['query_core_overlap']=.69
        self.assertFalse(sequence_core_support(self.f,{},self.p))

    def test_full_ligation_context_is_not_a_family_identity_requirement(self):
        self.f['source_active_center'].update(active_center_status='insufficient_alignment',template_span_alignment_status='insufficient')
        self.f['eligibility_reasons']=['reviewed_single_chain_center_or_span_not_supported']
        self.assertTrue(sequence_core_support(self.f,{},self.p))
        self.assertFalse(sequence_core_support(self.f,{},self.p,keep_center_context=True))

    def test_partial_competitor_and_unmapped_ligands_still_block(self):
        self.f['eligibility_reasons']=['partial_CDS'];self.assertFalse(sequence_core_support(self.f,{},self.p))
        self.f['eligibility_reasons']=[];self.f['competitor_margin']=0;self.assertFalse(sequence_core_support(self.f,{},self.p))
        self.f['competitor_margin']=20;self.f['source_active_center']['direct_ligand_status']='not_covered'
        self.assertFalse(sequence_core_support(self.f,{},self.p))

    def test_partial_feature_can_use_bound_source_proof_with_gene_aa_hash(self):
        d=Path(tempfile.mkdtemp());ann=d/'annotation';ann.mkdir();source=d/'input.fna';source.write_text('>ctg\nATGGTTTAA\n')
        (ann/'proteins.faa').write_text('>p1\nMV\n');(ann/'genes.tsv').write_text('id\tcontig\tstart\tend\tstrand\tpartial\n'+'p1\tctg\t1\t9\t1\tTrue\n')
        (ann/'prodigal.gff').write_text('##gff-version 3\n# Sequence Data: seqhdr\n# Model Data: x;transl_table=11;\nctg\tProdigal\tCDS\t1\t9\t.\t+\t0\tID=p1\n');(ann/'run_manifest.json').write_text(json.dumps({'status':'complete','input_files':[{'path':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()}]}))
        rows,_=coding_prepare(ann,d/'proof');_,loaded=load_coding(d/'proof',annotation=ann);row=loaded['p1']
        feature=copy.deepcopy(self.f);feature.update(protein_id='p1',eligibility_reasons=['partial_CDS']);gene={'contig':'ctg','start':1,'end':9,'strand':1,'protein_sha256':hashlib.sha256(b'MV').hexdigest()}
        self.assertTrue(sequence_core_support(feature,{},self.p,coding_completeness=row,gene=gene))
        bad=dict(gene,protein_sha256='0'*64);self.assertFalse(sequence_core_support(feature,{},self.p,coding_completeness=row,gene=bad))
        self.assertFalse(sequence_core_support(feature,{},self.p,coding_completeness=row))

    def test_family_profiles_are_not_interchangeable(self):
        self.f['family']='FeFe';self.assertFalse(sequence_core_support(self.f,{},self.p))
        self.f['domain_evidence'][0]['accession']='PF02906';self.assertTrue(sequence_core_support(self.f,{},self.p))

    def test_reviewed_noncanonical_pattern_is_separate_from_bare_substitution(self):
        self.f['source_active_center']['direct_ligand_status']='substituted'
        self.assertFalse(sequence_core_support(self.f,{},self.p))
        alt=dict(candidate_model_id='Aw_H_cubane_3C1Ser_literature',pattern_status='literature_ligand_pattern_candidate',query_partial='False')
        self.assertTrue(sequence_core_support(self.f,alt,self.p))
        alt['query_partial']='True';self.assertFalse(sequence_core_support(self.f,alt,self.p))

    def test_unconfirmed_reference_lineage_cannot_enable_stronger_transfer(self):
        review={'uncertain_positive_prototypes':[{'prototype_id':'p'}]}
        self.assertFalse(sequence_core_support(self.f,{},self.p,source_review=review))
        self.assertTrue(sequence_core_support(self.f,{},self.p))

    def test_source_scope_review_rejects_changed_primary_or_prototype(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'research_workflow').mkdir();source=root/'source.html';prototype=root/'prototypes.json'
            source.write_text('primary');prototype.write_text('[]')
            review=dict(source_path='source.html',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                prototype_path='prototypes.json',prototype_sha256=hashlib.sha256(prototype.read_bytes()).hexdigest())
            (root/'research_workflow/identity_scope_source_review.json').write_text(json.dumps(review))
            self.assertEqual(load_lineage_review(root),review)
            for path in [source,prototype]:
                original=path.read_bytes();path.write_bytes(original+b'changed')
                with self.assertRaises(ValueError):load_lineage_review(root)
                path.write_bytes(original)

    def complete_fefe(self):
        f=copy.deepcopy(self.f)
        f.update(family='FeFe',positive_score=31.8,competitor_margin=3.5,
            eligibility_reasons=['reviewed_single_chain_center_or_span_not_supported','required_domain_architecture_not_colocated_with_center'],
            positive_evidence=dict(core_coverage=.95,evalue=1e-25))
        f['source_active_center'].update(active_center_status='insufficient_alignment',template_span_alignment_status='insufficient')
        f['domain_evidence'][0].update(accession='PF02906',profile='Fe_hyd_lg_C',supported=True,model_coverage=.979,evidence_kind='single_domain')
        return f,dict(query_partial=False),dict(min_domain_coverage=.8,min_core_similarity=40,min_competitor_margin=5)

    def test_complete_FeFe_profile_can_support_family_without_upgrading_context(self):
        f,a,p=self.complete_fefe();before=copy.deepcopy(f)
        self.assertTrue(complete_fefe_core_support(f,a,p))
        self.assertTrue(sequence_core_support(f,a,p))
        self.assertFalse(sequence_core_support(f,a,p,keep_center_context=True))
        self.assertEqual(f,before)

    def test_complete_profile_path_is_bounded_to_FeFe_and_known_ligands(self):
        for change in ['NiFe','ligand','multichain','partial','unknown_partial','untrusted','competitor','weak_match','low_match_coverage']:
            f,a,p=self.complete_fefe();review=None
            if change=='NiFe':f['family']='NiFe';f['domain_evidence'][0]['accession']='PF00374'
            if change=='ligand':f['source_active_center']['direct_ligand_status']='substituted'
            if change=='multichain':f['source_active_center']['center_organization']='multiple_polypeptides'
            if change=='partial':a['query_partial']=True
            if change=='unknown_partial':a={}
            if change=='untrusted':review={'uncertain_positive_prototypes':[{'prototype_id':'p'}]}
            if change=='competitor':f['competitor_margin']=0
            if change=='weak_match':f['positive_evidence']['evalue']=1.1e-10
            if change=='low_match_coverage':f['positive_evidence']['core_coverage']=.69
            self.assertFalse(sequence_core_support(f,a,p,source_review=review),change)

    def test_complete_profile_path_does_not_join_fragments_or_waive_GA(self):
        for change in ['coverage','core_overlap','GA','fragments','not_supported','partial_reason']:
            f,a,p=self.complete_fefe();d=f['domain_evidence'][0]
            if change=='coverage':d['model_coverage']=.799
            if change=='core_overlap':d['query_core_overlap']=.699
            if change=='GA':d['components'][0]['domain_bitscore']=24
            if change=='fragments':d['components']*=2;d['evidence_kind']='collinear_domain_segments'
            if change=='not_supported':d['supported']=False
            if change=='partial_reason':f['eligibility_reasons'].append('partial_CDS')
            self.assertFalse(sequence_core_support(f,a,p),change)

    def _actual_fefe_profile_evidence(self):
        """Run the profile helper against the real R80 candidate, not a mock row."""
        root = Path(__file__).resolve().parents[1]
        annotation = root/'results/round80_anonymous955_v1/annotation'
        feature = next(f for f in json.loads((annotation/'confirmation/features.json').read_text())
                       if f['protein_id'] == 'anon955_0023')
        seq = next(r for r in annotation.joinpath('proteins.faa').read_text().split('>') if r.startswith('anon955_0023\n'))
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory); ann = tmp/'annotation'; ann.mkdir()
            (ann/'proteins.faa').write_text('>'+seq)
            evidence, _ = prepare_fefe_profile(ann, {'anon955_0023': feature}, tmp/'profile', {}, root=root)
            return feature, evidence['anon955_0023']

    def test_real_anon955_profile_fallback_is_bounded_and_keeps_center_uncovered(self):
        feature, evidence = self._actual_fefe_profile_evidence()
        params = dict(min_domain_coverage=.8, min_core_similarity=40, min_competitor_margin=5)
        alternative = dict(query_partial=False, FeFe_profile_ligands=evidence)
        self.assertEqual(feature['source_active_center']['direct_ligand_status'], 'not_covered')
        self.assertTrue(sequence_core_support(feature, alternative, params))
        self.assertFalse(sequence_core_support(feature, alternative, params, keep_center_context=True))

    def test_real_profile_fallback_rejects_old_conflicts_partial_lineage_and_bad_evidence(self):
        feature, evidence = self._actual_fefe_profile_evidence()
        params = dict(min_domain_coverage=.8, min_core_similarity=40, min_competitor_margin=5)
        base = dict(query_partial=False, FeFe_profile_ligands=evidence)
        cases = []
        f = copy.deepcopy(feature); f['source_active_center']['direct_ligand_status'] = 'substituted'; cases.append((f, base, None))
        f = copy.deepcopy(feature); f['source_active_center']['direct_ligand_status'] = 'ambiguous'; cases.append((f, base, None))
        f = copy.deepcopy(feature); a = copy.deepcopy(base); a['query_partial'] = True; cases.append((f, a, None))
        f = copy.deepcopy(feature); f['competitor_margin'] = 0; cases.append((f, base, None))
        f = copy.deepcopy(feature); cases.append((f, base, {'uncertain_positive_prototypes':[{'prototype_id':feature['positive_prototype']}]}))
        for label, mutate in [('posterior', lambda e: e['sites'][0].update(pp='.', posterior_status='low_or_missing')),
                              ('hash', lambda e: e.update(hmm_sha256='wrong')),
                              ('domain', lambda e: e['domain_ownership'].update(query_start=0))]:
            a = copy.deepcopy(base); mutate(a['FeFe_profile_ligands']); cases.append((copy.deepcopy(feature), a, None))
        for f, a, review in cases:
            self.assertFalse(sequence_core_support(f, a, params, source_review=review))


class NiFeIndeterminateCompetitionTests(unittest.TestCase):
    def setUp(self):
        self.root=Path(__file__).resolve().parents[1]
        path=self.root/'results/round140_TH1_counter_genome_v1/annotation/confirmation/features.json'
        self.feature=next(f for f in json.loads(path.read_text()) if f['protein_id']=='AP018558.1_1409')
        self.parameters=json.loads((self.root/'config/catalytic_confirmation_thresholds.json').read_text())['families']['NiFe']
        self.alternative={'query_partial':False}

    def test_actual_TH1_intact_core_survives_unknown_reference_competition(self):
        from explain_candidates import sequence_status,identity_reasons
        before=copy.deepcopy(self.feature)
        self.assertLess(self.feature['competitor_margin'],0)
        self.assertTrue(sequence_core_support(self.feature,self.alternative,self.parameters))
        self.assertEqual(sequence_status({},self.feature,self.alternative,self.parameters),'sequence_supported_uncalibrated_core')
        reason=next(r for r in identity_reasons(self.feature,self.parameters) if r['code']=='nonpositive_competitor_margin')
        self.assertEqual(reason['type'],'truth_unavailable')
        self.assertEqual(reason['claim'],'calibrated_high_identity')
        self.assertEqual(self.feature,before)

    def test_confirmed_or_untyped_competitors_still_block(self):
        for label,family in [('negative','NiFe'),('positive','NiFe'),('indeterminate','FeFe'),('','NiFe')]:
            f=copy.deepcopy(self.feature);f['competitor_evidence'].update(label=label,family=family)
            self.assertFalse(sequence_core_support(f,self.alternative,self.parameters),(label,family))
        f=copy.deepcopy(self.feature);f['competitor_evidence']=None
        self.assertFalse(sequence_core_support(f,self.alternative,self.parameters))

    def test_unknown_competition_never_substitutes_for_catalytic_evidence(self):
        cases=[]
        for status in ['substituted','not_covered','ambiguous']:
            f=copy.deepcopy(self.feature);f['source_active_center']['direct_ligand_status']=status;cases.append(f)
        for field,value in [('template_span_alignment_status','insufficient'),('context_alignment_status','insufficient'),
                            ('center_organization','multiple_polypeptides'),('site_contexts','[]')]:
            f=copy.deepcopy(self.feature);f['source_active_center'][field]=value;cases.append(f)
        f=copy.deepcopy(self.feature);f['domain_evidence'][0]['accession']='PF00346';cases.append(f)
        f=copy.deepcopy(self.feature);f['domain_evidence'][0]['components'][0]['domain_bitscore']=0;cases.append(f)
        f=copy.deepcopy(self.feature);f['domain_evidence'][0]['supported']=False;cases.append(f)
        f=copy.deepcopy(self.feature);f['positive_evidence']['core_coverage']=.7;cases.append(f)
        f=copy.deepcopy(self.feature);f['core_domain_coverage']=.7;cases.append(f)
        f=copy.deepcopy(self.feature);f['positive_score']=0;cases.append(f)
        f=copy.deepcopy(self.feature);f['eligibility_reasons']=['partial_CDS'];cases.append(f)
        f=copy.deepcopy(self.feature);f['positive_evidence']['prototype_id']='other';cases.append(f)
        for f in cases:self.assertFalse(sequence_core_support(f,self.alternative,self.parameters))

    def test_uncertain_positive_lineage_does_not_gain_catalytic_transfer(self):
        from explain_candidates import sequence_status
        review={'uncertain_positive_prototypes':[{'prototype_id':self.feature['positive_prototype']}]}
        self.assertFalse(sequence_core_support(self.feature,self.alternative,self.parameters,source_review=review))
        self.assertEqual(sequence_status({},self.feature,self.alternative,self.parameters,source_review=review),'hydrogenase_like_lineage_unconfirmed')


class FeFeIndeterminateCompetitionTests(unittest.TestCase):
    def setUp(self):
        root=Path(__file__).resolve().parents[1]
        path=root/'results/round117_function_conversion_v4/CpIII_B_whole_genome/annotation/confirmation/features.json'
        self.feature=next(f for f in json.loads(path.read_text()) if f['protein_id']=='CP009268.1_3397')
        self.parameters=json.loads((root/'config/catalytic_confirmation_thresholds.json').read_text())['families']['FeFe']

    def test_complete_FeFe_core_is_not_denied_by_exact_unclassified_reference(self):
        from explain_candidates import sequence_status,identity_reasons
        original=copy.deepcopy(self.feature)
        self.assertEqual(self.feature['competitor_evidence']['similarity'],100)
        self.assertTrue(sequence_core_support(self.feature,{},self.parameters))
        self.assertEqual(sequence_status({},self.feature,{},self.parameters),'sequence_supported_uncalibrated_core')
        self.assertTrue(any(r['type']=='truth_unavailable' for r in identity_reasons(self.feature,self.parameters)))
        self.assertEqual(self.feature,original)

    def test_FeFe_core_still_requires_canonical_catalytic_evidence(self):
        changes=[('competitor','negative'),('family','NiFe'),('ligand','substituted'),('span','insufficient'),('profile','PF02256'),('partial',True)]
        for key,value in changes:
            f=copy.deepcopy(self.feature)
            if key=='competitor':f['competitor_evidence']['label']=value
            elif key=='family':f['competitor_evidence']['family']=value
            elif key=='ligand':f['source_active_center']['direct_ligand_status']=value
            elif key=='span':f['source_active_center']['template_span_alignment_status']=value
            elif key=='profile':f['domain_evidence']=[d for d in f['domain_evidence'] if d['accession']==value]
            else:f['eligibility_reasons']=['partial_CDS']
            self.assertFalse(sequence_core_support(f,{},self.parameters),key)


class RealCpIIIC222SProfileIndeterminateTests(unittest.TestCase):
    """A published active engineered variant exercises the shared profile path."""

    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parents[1]
        run = cls.root/'results/round151_CpIII_C222S_v1'
        cls.feature = next(f for f in json.loads(
            (run/'annotation/confirmation/features.json').read_text())
            if f['protein_id'] == 'CpIII_C222S')
        cls.parameters = json.loads(
            (cls.root/'config/catalytic_confirmation_thresholds.json').read_text()
        )['families']['FeFe']
        # candidate_progression is the recorded pre-fix result; it contains
        # the exact hash-bound profile proof emitted for this run.
        with (run/'functions/candidate_progression.tsv').open(newline='') as handle:
            cells = next(r for r in csv.reader(handle, delimiter='\t')
                         if r and r[0] == 'CpIII_C222S')
        cls.progression = json.loads(cells[5])
        cls.evidence = cls.progression['FeFe_profile_ligands']
        with (run/'annotation/proteins.faa').open() as handle:
            query=SeqIO.read(handle,'fasta')
        cls.gene = {'protein_sha256': hashlib.sha256(str(query.seq).encode()).hexdigest()}
        assert cls.gene['protein_sha256']==cls.evidence['query_sha256']

    def alternative(self):
        return dict(query_partial=False, FeFe_profile_ligands=copy.deepcopy(self.evidence))

    def test_real_CpIII_C222S_profile_proof_supports_uncalibrated_core(self):
        alternative = self.alternative()
        self.assertTrue(sequence_core_support(self.feature, alternative, self.parameters,
                                              gene=self.gene))
        from explain_candidates import sequence_status
        self.assertEqual(sequence_status({}, self.feature, alternative, self.parameters,
                                         gene=self.gene), 'sequence_supported_uncalibrated_core')

    def test_real_profile_path_rejects_negative_or_foreign_competitor(self):
        for label, family in [('negative', 'FeFe'), ('indeterminate', 'NiFe')]:
            feature = copy.deepcopy(self.feature)
            feature['competitor_evidence'].update(label=label, family=family)
            alternative=self.alternative()
            from fefe_profile_ligands import _digest
            alternative['FeFe_profile_ligands']['feature_sha256']=_digest(feature)
            self.assertFalse(sequence_core_support(feature, alternative,
                                                    self.parameters, gene=self.gene),
                             (label, family))

    def test_real_profile_path_rejects_ligand_conflict_or_missing_site(self):
        for mutate in (
            lambda f, e: f['source_active_center'].update(direct_ligand_status='substituted'),
            lambda f, e: e['sites'].pop(),
            lambda f, e: e['sites'][0].update(residue='S', residue_status='substituted'),
        ):
            feature = copy.deepcopy(self.feature)
            alternative = self.alternative()
            mutate(feature, alternative['FeFe_profile_ligands'])
            from fefe_profile_ligands import _digest
            alternative['FeFe_profile_ligands']['feature_sha256']=_digest(feature)
            self.assertFalse(sequence_core_support(feature, alternative, self.parameters,
                                                    gene=self.gene))

    def test_real_profile_path_rejects_stale_proof_and_unproven_partial(self):
        for field, value in [('query_sha256', '0'*64), ('feature_sha256', '0'*64),
                             ('hmm_sha256', '0'*64), ('reference_model_sha256', '0'*64)]:
            alternative = self.alternative()
            alternative['FeFe_profile_ligands'][field] = value
            self.assertFalse(sequence_core_support(self.feature, alternative, self.parameters,
                                                    gene=self.gene), field)
        alternative = self.alternative()
        alternative['query_partial'] = True
        self.assertFalse(sequence_core_support(self.feature, alternative, self.parameters,
                                                gene=self.gene))

    def test_current_bound_domain_still_requires_GA_and_completeness(self):
        from fefe_profile_ligands import _digest
        for kind in ['GA','coverage']:
            feature=copy.deepcopy(self.feature);alternative=self.alternative()
            domain=feature['domain_evidence'][0]
            if kind=='GA':domain['components'][0]['domain_bitscore']=0
            else:domain['model_coverage']=.5
            alternative['FeFe_profile_ligands']['feature_sha256']=_digest(feature)
            alternative['FeFe_profile_ligands']['domain_ownership']['domain_sha256']=_digest(domain)
            self.assertFalse(sequence_core_support(feature,alternative,self.parameters,gene=self.gene),kind)


class NiFeBackboneTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=Path(__file__).resolve().parents[1]
        cls.annotation=cls.root/'results/round74_HsuL2_A_v1/annotation'
        cls.feature=json.loads((cls.annotation/'confirmation/features.json').read_text())[0]
        cls.parameters=json.loads((cls.root/'config/catalytic_confirmation_thresholds.json').read_text())['families']['NiFe']
        cls.backbones,cls.bound=reference_backbones(cls.annotation,{'query001':cls.feature},cls.root)

    def inputs(self):
        return copy.deepcopy(self.feature),dict(query_partial=False,NiFe_reference_backbone=copy.deepcopy(self.backbones['query001']))

    def test_actual_exposed_no_GA_case_has_supported_alternate_but_not_complete_center(self):
        f,a=self.inputs();before=copy.deepcopy(f)
        self.assertFalse(local_domain_evidence(f))
        self.assertAlmostEqual(a['NiFe_reference_backbone']['paired_reference_coverage'],568/596)
        self.assertTrue(sequence_core_support(f,a,self.parameters))
        self.assertFalse(sequence_core_support(f,a,self.parameters,keep_center_context=True))
        self.assertEqual(f,before)

    def test_reference_competition_ligands_context_and_partial_are_independent_requirements(self):
        for change in ['unknown_partial','partial','competitor','score','missing_backbone','backbone_fragment','wrong_reference','ligand','context','split','site_count','lineage']:
            f,a=self.inputs();review=None
            if change=='unknown_partial':a.pop('query_partial')
            if change=='partial':a['query_partial']=True
            if change=='competitor':f['competitor_margin']=0
            if change=='score':f['positive_score']=0
            if change=='missing_backbone':a.pop('NiFe_reference_backbone')
            if change=='backbone_fragment':a['NiFe_reference_backbone']['paired_reference_coverage']=.79
            if change=='wrong_reference':a['NiFe_reference_backbone']['reference_id']='other'
            if change=='ligand':f['source_active_center']['direct_ligand_status']='substituted'
            if change=='context':f['source_active_center']['context_alignment_status']='insufficient'
            if change=='split':f['source_active_center']['center_organization']='multiple_polypeptides'
            if change=='site_count':f['source_active_center']['site_contexts']='[]'
            if change=='lineage':review={'uncertain_positive_prototypes':[{'prototype_id':f['positive_prototype']}]}
            self.assertFalse(sequence_core_support(f,a,self.parameters,source_review=review),change)

    def test_raw_backbone_must_match_real_query_sequence(self):
        with tempfile.TemporaryDirectory() as directory:
            ann=Path(directory);(ann/'confirmation').mkdir()
            (ann/'proteins.faa').write_bytes((self.annotation/'proteins.faa').read_bytes())
            rows=(self.annotation/'confirmation/prototype_blast.tsv').read_text().splitlines()
            for i,row in enumerate(rows):
                cells=row.split('\t')
                if cells[1]==self.feature['positive_prototype']:
                    cells[12]='X'+cells[12][1:];rows[i]='\t'.join(cells)
            (ann/'confirmation/prototype_blast.tsv').write_text('\n'.join(rows)+'\n')
            with self.assertRaises(ValueError):reference_backbones(ann,{'query001':self.feature},self.root)


class NiFeCoreRegionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=Path(__file__).resolve().parents[1]
        cls.annotation=cls.root/'results/round120_core_region_baseline_v2/Tlitoralis_core_stop_missing/annotation'
        cls.features={f['protein_id']:f for f in json.loads((cls.annotation/'confirmation/features.json').read_text())}
        cls.parameters=dict(min_domain_coverage=.8,min_core_similarity=30,min_competitor_margin=5)
        cls.proofs,cls.inputs=nife_core_regions(cls.annotation,cls.features,cls.root)

    def test_actual_partial_stoploss_has_bound_region_proof(self):
        pid=next(iter(self.proofs)); proof=self.proofs[pid]; feature=self.features[pid]
        self.assertEqual(proof['status'],NIFE_CORE_REGION_STATUS)
        self.assertEqual(proof['recomputed_metrics']['window_coverage'],1.0)
        self.assertTrue(nife_core_region_support(feature,{'NiFe_core_region':proof},self.parameters,gene={'protein_sha256':proof['protein_sha256']}))
        self.assertFalse(nife_core_region_support(feature,{'NiFe_core_region_proof':proof},self.parameters,gene={'protein_sha256':proof['protein_sha256']}))

    def test_region_proof_rejects_hash_metrics_ligands_and_lineage(self):
        pid=next(iter(self.proofs)); feature=self.features[pid]; proof=self.proofs[pid]
        for key,value in [('protein_sha256','0'*64),('feature_sha256','0'*64)]:
            bad=copy.deepcopy(proof);bad[key]=value
            self.assertFalse(nife_core_region_support(feature,{'NiFe_core_region':bad},self.parameters,gene={'protein_sha256':proof['protein_sha256']}))
        for key,value in [('window_coverage',.9),('template_span_alignment_status','insufficient'),('context_alignment_status','insufficient'),('direct_ligand_status','not_covered')]:
            bad=copy.deepcopy(proof);bad['recomputed_metrics'][key]=value
            self.assertFalse(nife_core_region_support(feature,{'NiFe_core_region':bad},self.parameters,gene={'protein_sha256':proof['protein_sha256']}))
        for key,value in [('rule_version','wrong'),('center_organization','multiple_polypeptides')]:
            bad=copy.deepcopy(proof)
            if key=='rule_version': bad[key]=value
            else: bad['recomputed_metrics'][key]=value
            self.assertFalse(nife_core_region_support(feature,{'NiFe_core_region':bad},self.parameters,gene={'protein_sha256':proof['protein_sha256']}))
        bad=copy.deepcopy(feature);bad['competitor_margin']=4.9
        self.assertFalse(nife_core_region_support(bad,{'NiFe_core_region':proof},self.parameters,gene={'protein_sha256':proof['protein_sha256']}))
        bad=copy.deepcopy(proof);next(s for s in bad['sites'] if s.get('residue_role')=='direct_metal_ligand')['status']='substituted'
        self.assertFalse(nife_core_region_support(feature,{'NiFe_core_region':bad},self.parameters,gene={'protein_sha256':proof['protein_sha256']}))
        review={'uncertain_positive_prototypes':[{'prototype_id':feature['positive_prototype']}]}
        self.assertFalse(nife_core_region_support(feature,{'NiFe_core_region':proof},self.parameters,source_review=review,gene={'protein_sha256':proof['protein_sha256']}))

    def test_region_proof_keeps_partial_and_rejects_foreign_or_weak_feature(self):
        pid=next(iter(self.proofs)); feature=self.features[pid]; proof=self.proofs[pid]
        self.assertIn('partial_CDS',feature['eligibility_reasons'])
        for mutate in [lambda f:f.update(family='FeFe'),lambda f:f.update(positive_score=0),lambda f:f.update(competitor_margin=0),lambda f:f.update(eligibility_reasons=['partial_CDS','family_not_calibrated'])]:
            bad=copy.deepcopy(feature);mutate(bad)
            self.assertFalse(nife_core_region_support(bad,{'NiFe_core_region':proof},self.parameters,gene={'protein_sha256':proof['protein_sha256']}))

    def test_actual_hsp_revalidation_rejects_deleted_protein_region(self):
        import shutil
        with tempfile.TemporaryDirectory() as directory:
            ann=Path(directory)/'annotation';shutil.copytree(self.annotation,ann)
            fasta=ann/'proteins.faa'
            with fasta.open() as handle: records=list(SeqIO.parse(handle,'fasta'))
            target=records[0];target.seq=target.seq[:100]+target.seq[160:]
            with fasta.open('w') as handle: SeqIO.write(records,handle,'fasta')
            features={f['protein_id']:f for f in json.loads((ann/'confirmation/features.json').read_text())}
            proofs,_=nife_core_regions(ann,features,self.root)
            self.assertEqual(proofs,{})


if __name__=='__main__':unittest.main()
