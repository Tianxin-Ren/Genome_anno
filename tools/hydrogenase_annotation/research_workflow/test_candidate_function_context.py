import csv, json, hashlib, tempfile, unittest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from explain_candidates import _candidate_context, sha
from candidate_screen import catalog

def write_tsv(path, rows):
    fields=sorted({k for r in rows for k in r})
    with path.open('w') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)

class CandidateFunctionContextTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup); self.r=Path(self.tmp.name); self.a=self.r/'a';self.s=self.r/'s';self.c=self.r/'c';self.a.mkdir();self.s.mkdir();self.c.mkdir()
        (self.a/'run_manifest.json').write_text(json.dumps({'input_mode':'genome'})); write_tsv(self.a/'catalytic_anchors.tsv',[{'protein_id':'core'}]); write_tsv(self.a/'genes.tsv',[{'id':'extra','contig':'ctg','start':'1','end':'9','strand':'1','partial':'False'}]);
        SeqIO.write([SeqRecord(Seq('MKK'),id='extra')],self.a/'proteins.faa','fasta')
        screen=self.s/'candidate_screen.tsv';write_tsv(screen,[{'protein_id':'extra','source_sets':'["verified_GA_domain"]','discovery_status':'discovery_only','family_candidates':'["FeFe_related"]','rejected_reason':''}]);
        sm_inputs={str((self.a/'run_manifest.json').resolve()):sha(self.a/'run_manifest.json')}
        sm={'status':'complete','annotation':str(self.a.resolve()),'input_sha256':sm_inputs,'output_sha256':{'candidate_screen.tsv':sha(screen)} };(self.s/'candidate_screen_manifest.json').write_text(json.dumps(sm))
        self.feature={'protein_id':'extra','family':'FeFe','eligibility_reasons':[],'core_domain_coverage':1.0,'positive_score':100.0,'competitor_margin':10.0,'positive_prototype':'p','competitor_prototype':'n','source_active_center':{}}
        subset=self.c/'candidate_annotation';subset.mkdir();write_tsv(subset/'genes.tsv',[{'id':'extra','contig':'ctg','start':'1','end':'9','strand':'1','partial':'False'}]);SeqIO.write([SeqRecord(Seq('MKK'),id='extra')],subset/'proteins.faa','fasta');(subset/'candidate_subset_manifest.json').write_text(json.dumps({'source_annotation':str(self.a.resolve()),'candidate_proteins':['extra'],'genes_sha256':sha(subset/'genes.tsv'),'proteins_sha256':sha(subset/'proteins.faa')}))
        fp=subset/'features.json';fp.write_text(json.dumps([self.feature]));
        code=self.r/'code';model=self.r/'model';code.write_text('code');model.write_text('model')
        report={'status':'complete','annotation':str(self.a.resolve()),'source_screen_manifest_sha256':sha(self.s/'candidate_screen_manifest.json'),'source_screen_inputs':sm_inputs,'code_sha256':{str(code):sha(code)},'model_sha256':{str(model):sha(model)},'decisions':[{'candidate_id':'extra','family':'FeFe','normalized_families':['FeFe'],'sequence_status':'sequence_supported_uncalibrated_core','blocking_reasons':[],'features':self.feature}], 'evidence':{'extra':{'annotation':str(subset),'features':str(fp)}}}
        rp=self.c/'candidate_confirmation.json';rp.write_text(json.dumps(report));(self.c/'candidate_confirmation_manifest.json').write_text(json.dumps({'status':'complete','report_sha256':sha(rp)}))
    def tearDown(self): self.tmp.cleanup()
    def test_unconfirmed_without_evidence_is_ignored(self):
        d=json.loads((self.c/'candidate_confirmation.json').read_text());d['decisions']=[{'candidate_id':'extra','family':'FeFe','normalized_families':['FeFe'],'sequence_status':'unconfirmed','blocking_reasons':[]}];d['evidence']={};rp=self.c/'candidate_confirmation.json';rp.write_text(json.dumps(d));m=json.loads((self.c/'candidate_confirmation_manifest.json').read_text());m['report_sha256']=sha(rp);(self.c/'candidate_confirmation_manifest.json').write_text(json.dumps(m));self.assertEqual(_candidate_context(self.a.resolve(),self.s,self.c)[0],[])
    def test_accepted_missing_feature_fails_closed(self):
        fp=self.c/'candidate_annotation/features.json';fp.unlink();
        with self.assertRaises(ValueError): _candidate_context(self.a.resolve(),self.s,self.c)
    def _report(self): return json.loads((self.c/'candidate_confirmation.json').read_text())
    def _save_report(self, report):
        p=self.c/'candidate_confirmation.json';p.write_text(json.dumps(report));(self.c/'candidate_confirmation_manifest.json').write_text(json.dumps({'status':'complete','report_sha256':sha(p)}))
    def test_accepted_helper_loads_and_binds_subset(self):
        accepted,features,alts,bound=_candidate_context(self.a.resolve(),self.s,self.c)
        self.assertEqual([d['candidate_id'] for d in accepted],['extra']);self.assertIn(self.c/'candidate_annotation/genes.tsv',bound)
    def test_wrong_annotation_and_source_change_fail(self):
        d=self._report();d['annotation']=str((self.r/'other').resolve());self._save_report(d)
        with self.assertRaises(ValueError): _candidate_context(self.a.resolve(),self.s,self.c)
        d['annotation']=str(self.a.resolve());self._save_report(d);(self.a/'run_manifest.json').write_text('{}')
        with self.assertRaises(ValueError): _candidate_context(self.a.resolve(),self.s,self.c)
    def test_feature_report_and_declared_alternative_fail(self):
        d=self._report();d['decisions'][0]['features']['positive_score']=99;self._save_report(d)
        with self.assertRaises(ValueError): _candidate_context(self.a.resolve(),self.s,self.c)
        self.setUp();d=self._report();d['evidence']['extra']['noncanonical']=str(self.r/'missing.tsv');self._save_report(d)
        with self.assertRaises(ValueError): _candidate_context(self.a.resolve(),self.s,self.c)
    def test_duplicate_features_and_partial_mismatch_fail(self):
        fp=self.c/'candidate_annotation/features.json';fp.write_text(json.dumps([self.feature,self.feature]));
        with self.assertRaises(ValueError): _candidate_context(self.a.resolve(),self.s,self.c)
    def test_existing_ni_fe_rule_and_resource_guard_are_real(self):
        from pyridine_function import load_model, predict, guard_cross_anchor_coupling
        def hit(role, start=1, end=100): return {'protein_id':role,'identity':55,'query_coverage':.95,'reference_coverage':.95,'partial':False,'evalue_verified':True,'domain_organization_supported':True,'query_domain_start':start,'query_domain_end':end}
        roles={r:[hit(r)] for r in ('SHI_alpha','SHI_delta','SHI_beta','SHI_gamma')}
        ev={'structure_id':'NiFe_HyhNADPH','identity_status':'sequence_supported_uncalibrated_core','roles':roles,'association_supported':True,'coupling_linkage_supported':True,'partial_members':False,'module_competition_audited':True}
        result=predict(ev,load_model());self.assertEqual(result['coupling_potential_status'],'local_coupling_potential_uncalibrated')
        rival={'evidence':{'roles':roles},'coupling_potential_status':'local_coupling_potential_uncalibrated','coupling_blockers':[]}
        guard_cross_anchor_coupling({'extra':result,'rival':rival})
        self.assertEqual(result['coupling_potential_status'],'unresolved')
        self.setUp();write_tsv(self.c/'candidate_annotation/genes.tsv',[{'id':'extra','contig':'ctg','start':'1','end':'9','strand':'1','partial':'True'}]);
        m=json.loads((self.c/'candidate_annotation/candidate_subset_manifest.json').read_text());m['genes_sha256']=sha(self.c/'candidate_annotation/genes.tsv');(self.c/'candidate_annotation/candidate_subset_manifest.json').write_text(json.dumps(m))
        with self.assertRaises(ValueError): _candidate_context(self.a.resolve(),self.s,self.c)

    def _function_fixture(self, identity='sequence_supported_uncalibrated_core'):
        report=self._report()
        report['decisions'][0].update(identity_status=identity,sequence_status=identity,function_status='unassigned')
        self._save_report(report)
        folder=self.r/'functions';folder.mkdir()
        write_tsv(folder/'functional_annotations.tsv',[dict(candidate_id='extra',identity_status=identity,
            function_prediction='candidate_function_rule_missing',confidence='unresolved',supported_attributes='[]')])
        sm=json.loads((self.s/'candidate_screen_manifest.json').read_text())
        metadata=dict(input_sha256=sm['input_sha256'],functional_annotations_sha256=sha(folder/'functional_annotations.tsv'),
            candidate_scoped_research_predictions=['extra'],candidate_context_metadata=dict(
                confirmation=str(self.c),confirmation_manifest_sha256=sha(self.c/'candidate_confirmation_manifest.json')))
        (folder/'manifest.json').write_text(json.dumps(metadata))
        return folder

    def test_catalog_joins_nonempty_extra_rows_and_requires_context(self):
        folder=self._function_fixture()
        rows=catalog(self.s,folder,self.r/'catalog.tsv',self.c)
        self.assertEqual(rows[0]['candidate_id'],'extra')
        self.assertEqual(rows[0]['function_confidence'],'unresolved')
        self.assertNotIn('candidate_scoped_research_prediction;',rows[0]['evidence'])
        metadata=json.loads((folder/'manifest.json').read_text());metadata.pop('candidate_context_metadata')
        (folder/'manifest.json').write_text(json.dumps(metadata))
        with self.assertRaisesRegex(ValueError,'context missing'):
            catalog(self.s,folder,self.r/'invalid.tsv',self.c)

    def test_catalog_rejects_matching_but_unconfirmed_extra_identity(self):
        folder=self._function_fixture(identity='unconfirmed')
        with self.assertRaisesRegex(ValueError,'Unconfirmed candidate'):
            catalog(self.s,folder,self.r/'invalid.tsv',self.c)

if __name__=='__main__': unittest.main()
