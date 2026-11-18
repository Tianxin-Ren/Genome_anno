import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from candidate_screen import screen, catalog, ROOT, sha, _raw_blast_hits, _raw_target_hits
from identity_scope import NIFE_CORE_REGION_RULE, NIFE_CORE_REGION_STATUS
from discovery_evidence import load_domains
from Bio import SeqIO


def save(path,rows):
    with path.open('w') as handle:
        fields=sorted({k for r in rows for k in r});writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)


class CandidateScreenTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();root=Path(self.tmp.name).resolve()
        self.a=root/'a';self.d=root/'d';self.out=root/'out';self.a.mkdir();self.d.mkdir()
        self.dman={'source_runs':[str(self.a)]}
        run=dict(status='complete',mode='genome',thresholds=dict(hmm_domain_evalue=1e-5,hmm_coverage=.55),
            discovery_profile_manifest_sha256=sha(ROOT/'references/processed/profile_manifest.json'),discovery_evidence_policy_sha256=sha(ROOT/'config/discovery_evidence.json'))
        (self.a/'run_manifest.json').write_text(json.dumps(run));(self.d/'query_domain_manifest.json').write_text(json.dumps(self.dman))
        save(self.a/'catalytic_anchors.tsv',[dict(protein_id='anchor',catalytic_identity_confidence='calibrated_high',function_status='unassigned')])
        save(self.a/'family_candidates.tsv',[dict(protein_id='anchor',family='NiFe')])
        save(self.a/'rejected.tsv',[dict(protein_id='rejected',reason='non_catalytic_homolog_wins')])
        save(self.a/'genes.tsv',[dict(id=p,contig='c',start=1,end=1000,strand=1,coordinate_evidence='source_CDS') for p in ['anchor','ga','nuo','shared','weak','rejected']])
        (self.a/'proteins.faa').write_text(''.join(f'>{p}\n{"M" * n}\n' for p,n in [('anchor',400),('ga',417),('nuo',200),('shared',200),('weak',100),('rejected',200)]))
        # An empty saved target table models the normal, already-evaluated
        # cache. Recovery-specific tests remove it explicitly.
        (self.a/'family.tbl').write_text('# empty\n')
        chain=dict(accession='PF02906',supported=True,model_coverage=.8321,components=[dict(query_start=46,query_end=197,domain_bitscore=120.8),dict(query_start=253,query_end=342,domain_bitscore=59.5)])
        save(self.d/'query_domain_organizations.tsv',[
            dict(source_run=str(self.a),protein_id='ga',profile_evidence=json.dumps([chain])),
            dict(source_run=str(self.a),protein_id='nuo',profile_evidence=json.dumps([dict(accession='PF00346',supported=True,components=[dict(query_start=1,query_end=200)])])),
            dict(source_run=str(self.a),protein_id='shared',profile_evidence=json.dumps([dict(accession='PF00175',supported=True,components=[])]))])
        def raw(pid,profile,qlen,mlen,ms,me,qs,qe,ev):
            return f'{pid} - {qlen} {profile} - {mlen} 1e-30 100 0 1 1 1e-30 {ev} 100 0 {ms} {me} {qs} {qe} {qs} {qe} 0.99 -\n'
        (self.a/'family.domtbl').write_text(raw('ga','FeFe',417,570,1,185,10,199,'5e-35')+raw('ga','FeFe',417,570,300,414,294,408,'7e-20')+
            raw('weak','FeFe',100,570,1,20,1,20,'1e-2')+raw('anchor','NiFe',400,466,1,220,1,190,'1e-30')+raw('anchor','NiFe',400,466,260,459,201,390,'1e-30'))
    def tearDown(self): self.tmp.cleanup()
    def run_screen(self):
        with patch('candidate_screen.verify_completed',return_value=self.dman): return screen(self.a,self.d,self.out)
    def test_combined_GA_retains_candidate_below_whole_HMM_coverage(self):
        self.run_screen();rows={r['protein_id']:r for r in self.read('candidate_screen.tsv')}
        self.assertEqual(rows['ga']['core_identity'],'not_assessed')
        self.assertEqual(rows['ga']['potential_only'],'True')
        self.assertEqual(len(json.loads(rows['ga']['domain_evidence'])[0]['components']),2)
        self.assertNotIn('shared',rows);self.assertEqual(rows['anchor']['core_identity'],'calibrated_high')
    def read(self,name):
        with (self.out/name).open() as handle:return list(csv.DictReader(handle,delimiter='\t'))
    def test_Nuo_related_and_rejected_are_not_confirmed_hydrogenases(self):
        self.run_screen();rows={r['protein_id']:r for r in self.read('candidate_screen.tsv')}
        self.assertEqual(json.loads(rows['nuo']['family_candidates']),['NiFe_or_Nuo_related'])
        self.assertEqual(rows['nuo']['core_function'],'unassigned')
        self.assertEqual(rows['rejected']['discovery_status'],'noncatalytic_homolog_favored_retained')
    def test_raw_parser_and_saved_chain_rules_keep_real_coordinates(self):
        self.run_screen();rows=self.read('family_domain_ledger.tsv')
        self.assertEqual(rows[0]['model_length'],'570');self.assertEqual(rows[0]['query_length'],'417')
        self.assertAlmostEqual(float(rows[0]['profile_coverage']),185/570)
        self.assertTrue(all(r['retained_by_core_discovery']=='True' for r in rows if r['protein_id']=='anchor'))
        self.assertTrue(all(r['retained_by_core_discovery']=='False' for r in rows if r['protein_id']=='ga'))
        weak={r['protein_id']:r for r in self.read('candidate_screen.tsv')}['weak']
        self.assertEqual(weak['discovery_status'],'weak_or_partial_signal')
        self.assertEqual(json.loads(weak['family_candidates']),['FeFe'])
        self.assertEqual(weak['core_identity'],'not_assessed')
        self.assertEqual(weak['core_function'],'unassigned')
        self.assertEqual(weak['potential_only'],'True')
        self.assertIn('raw_discovery_signal',json.loads(weak['source_sets']))
        self.assertEqual(json.loads(weak['domain_evidence'])[0]['filter_reason'],
                         'above_saved_E_threshold;individual_domain_below_saved_coverage')
    def test_source_mismatch_fails_before_output(self):
        with patch('candidate_screen.verify_completed',return_value={'source_runs':['elsewhere']}):
            with self.assertRaises(ValueError): screen(self.a,self.d,self.out)
        self.assertFalse(self.out.exists())
    def test_no_overwrite_and_output_hashes(self):
        report=self.run_screen()
        self.assertTrue(all(sha(self.out/p)==d for p,d in report['output_sha256'].items()))
        with self.assertRaises(FileExistsError): self.run_screen()

    def test_catalog_accepts_native_or_extra_region_decision_without_function_promotion(self):
        self.run_screen()
        confirmation = self.out / 'confirmation'; confirmation.mkdir()
        features = confirmation / 'features.json'; features.write_text(json.dumps([{'protein_id': 'anchor'}]))
        proof = {'status': NIFE_CORE_REGION_STATUS, 'rule_version': NIFE_CORE_REGION_RULE,
                 'protein_id': 'anchor', 'protein_sha256': 'x'}
        decision = {'candidate_id': 'anchor', 'identity_status': NIFE_CORE_REGION_STATUS,
                    'sequence_status': NIFE_CORE_REGION_STATUS, 'function_status': 'unassigned',
                    'blocking_reasons': [], 'NiFe_core_region': proof}
        report = {'source_screen_manifest_sha256': sha(self.out/'candidate_screen_manifest.json'),
                  'decisions': [decision], 'evidence': {'anchor': {'features': str(features)}},
                  'NiFe_core_region_input_sha256': {str(self.a/'genes.tsv'): sha(self.a/'genes.tsv')}}
        report_path = confirmation/'candidate_confirmation.json'; report_path.write_text(json.dumps(report))
        (confirmation/'candidate_confirmation_manifest.json').write_text(json.dumps({'report_sha256': sha(report_path)}))
        function = self.out / 'functions'; function.mkdir()
        function_table = function/'functional_annotations.tsv'; function_table.write_text('candidate_id\n')
        source_manifest = str(self.a/'run_manifest.json')
        fmanifest = {'input_sha256': {source_manifest: sha(self.a/'run_manifest.json')},
                     'functional_annotations_sha256': sha(function_table)}
        (function/'manifest.json').write_text(json.dumps(fmanifest))
        with patch('candidate_screen.nife_core_region_support', return_value=True), patch('candidate_screen.nife_core_regions', return_value=({'anchor': proof}, {})):
            rows = catalog(self.out, function, self.out/'catalog', confirmation_dir=confirmation)
        anchor = next(r for r in rows if r['candidate_id'] == 'anchor')
        self.assertEqual(anchor['function_prediction'], 'NiFe_catalytic_core_region_supported_function_unresolved')
        self.assertEqual(anchor['function_confidence'], 'uncalibrated_core_region')
        self.assertNotEqual(anchor['function_prediction'], 'reference_sequence_matched_reaction_potential')
        decision['candidate_id'] = 'ga'; proof['protein_id'] = 'ga'
        features.write_text(json.dumps([{'protein_id': 'ga'}]))
        report['decisions'] = [decision]; report['evidence'] = {'ga': {'features': str(features)}}
        report_path.write_text(json.dumps(report)); (confirmation/'candidate_confirmation_manifest.json').write_text(json.dumps({'report_sha256': sha(report_path)}))
        with patch('candidate_screen.nife_core_region_support', return_value=True), patch('candidate_screen.nife_core_regions', return_value=({'ga': proof}, {})):
            rows = catalog(self.out, function, self.out/'catalog_extra', confirmation_dir=confirmation)
        extra = next(r for r in rows if r['candidate_id'] == 'ga')
        self.assertEqual(extra['function_prediction'], 'NiFe_catalytic_core_region_supported_function_unresolved')
        proof['protein_sha256'] = '0' * 64
        report_path.write_text(json.dumps(report)); (confirmation/'candidate_confirmation_manifest.json').write_text(json.dumps({'report_sha256': sha(report_path)}))
        with self.assertRaisesRegex(ValueError, 'does not validate'):
            catalog(self.out, function, self.out/'catalog_bad_proof', confirmation_dir=confirmation)

    def test_raw_domain_rejects_nonfinite_scores_and_invalid_coordinates(self):
        def row(evalue='nan', model_start='1'):
            fields=['p','-', '100','PF02906','-','570','1e-30','100','0','1','1','1e-30',
                    evalue,'20','0',model_start,'20','1','20','1','20','-']
            return ' '.join(fields)+'\n'
        path=self.a/'malformed.domtbl'
        path.write_text(row())
        with self.assertRaisesRegex(ValueError,'Invalid HMM coordinates or scores'):
            load_domains(path,{'PF02906':'FeFe'})
        path.write_text(row(evalue='1e-20',model_start='0'))
        with self.assertRaisesRegex(ValueError,'Invalid HMM coordinates or scores'):
            load_domains(path,{'PF02906':'FeFe'})

    def test_reported_blast_and_whole_target_signals_are_typed_and_bound(self):
        refs=json.loads((ROOT/'references/processed/references.json').read_text())
        ref=next(r for r in refs if r.get('kind') in {'catalytic','unresolved','decoy'})
        refseq=next(r for r in SeqIO.parse(ROOT/'references/processed/references.faa','fasta') if r.id==ref['reference_id'])
        run=dict(thresholds=dict(blast_evalue=1e-5,min_identity=0,min_query_coverage=0,min_reference_coverage=0),
                 references_json_sha256=sha(ROOT/'references/processed/references.json'),reference_fasta_sha256=sha(ROOT/'references/processed/references.faa'))
        values=['weak',ref['reference_id'],'50','10','1','10','1','10','1e-2','20','100',str(len(refseq.seq)),'M'*10,str(refseq.seq[:10])]
        (self.a/'blast.tsv').write_text('\t'.join(values)+'\n')
        hits,_,_= _raw_blast_hits(self.a,run,{'weak':{}})
        self.assertEqual(hits[0]['qseqid'],'weak')
        values[10]='10'
        (self.a/'blast.tsv').write_text('\t'.join(values)+'\n')
        with self.assertRaisesRegex(ValueError,'Invalid BLAST'):
            _raw_blast_hits(self.a,run,{'weak':{}})
        (self.a/'family.tbl').write_text('weak - FeFe - 1e-2 20 0 1 1 1 1 1 1 1 1 1 1 1 1 1\n')
        targets,_=_raw_target_hits(self.a,{'weak':{}})
        self.assertEqual(targets[0]['family'],'FeFe')

    def test_reported_blast_rejects_nonfinite_score_and_unknown_owner(self):
        refs=json.loads((ROOT/'references/processed/references.json').read_text())
        ref=next(r for r in refs if r.get('kind') in {'catalytic','unresolved','decoy'})
        refseq=next(r for r in SeqIO.parse(ROOT/'references/processed/references.faa','fasta') if r.id==ref['reference_id'])
        run=dict(thresholds=dict(blast_evalue=1e-5,min_identity=0,min_query_coverage=0,min_reference_coverage=0),
                 references_json_sha256=sha(ROOT/'references/processed/references.json'),reference_fasta_sha256=sha(ROOT/'references/processed/references.faa'))
        values=['unknown',ref['reference_id'],'50','10','1','10','1','10','nan','20','100',str(len(refseq.seq)),'M'*10,str(refseq.seq[:10])]
        (self.a/'blast.tsv').write_text('\t'.join(values)+'\n')
        with self.assertRaisesRegex(ValueError,'ownership|score'):
            _raw_blast_hits(self.a,run,{'weak':{}})

    def test_auxiliary_whole_target_signal_enters_weak_queue(self):
        save(self.d/'query_target_signals.tsv',[dict(source_run=str(self.a),protein_id='shared',query_id='q',accession='PF00374.1',profile='NiFeSe_Hases',sequence_evalue='1e-4',sequence_bitscore='12',sequence_bias='0',model_length='466',query_length='200',model_sha256='x',source='aux')])
        self.run_screen(); rows={r['protein_id']:r for r in self.read('candidate_screen.tsv')}
        self.assertIn('shared',rows)
        self.assertEqual(rows['shared']['discovery_status'],'weak_or_partial_signal')
        self.assertIn('whole_target',json.loads(rows['shared']['raw_signal_evidence']))
        self.assertEqual(rows['shared']['core_function'],'unassigned')

    def _bind_recovery_run(self, model=None):
        run=json.loads((self.a/'run_manifest.json').read_text())
        run.update(model_sha256=model or sha(ROOT/'references/processed/catalytic_families.hmm'),
                   discovery_evidence_script_sha256=sha(ROOT/'scripts/discovery_evidence.py'))
        (self.a/'run_manifest.json').write_text(json.dumps(run))

    def _fake_recovery(self, fasta, out, root, threads, thresholds):
        out.mkdir(parents=True,exist_ok=True)
        profile=json.loads((ROOT/'references/processed/profile_manifest.json').read_text())[0]['profile']
        (out/'family.tbl').write_text(f'weak - {profile} - 1e-2 20 0 1 1 1 1 1 1 1 1 1 1 1 1 1\n')
        (out/'family.domtbl').write_text('# no domain rows\n')
        (out/'hmmsearch.log').write_text('fake\n')
        return [],ROOT/'references/processed/profile_manifest.json',ROOT/'config/discovery_evidence.json'

    def test_missing_family_table_recovers_only_target_signal_without_source_mutation(self):
        self._bind_recovery_run();(self.a/'family.tbl').unlink()
        source_hash=sha(self.a/'proteins.faa')
        with patch('candidate_screen.verify_completed',return_value=self.dman), patch('candidate_screen.search_discovery',side_effect=self._fake_recovery):
            report=screen(self.a,self.d,self.out)
        rows={r['protein_id']:r for r in self.read('candidate_screen.tsv')}
        self.assertIn('weak',rows)
        self.assertEqual(report['whole_target_source'],'recovered_table')
        self.assertEqual(json.loads(rows['weak']['raw_signal_evidence'])['whole_target'][0]['family'],json.loads((ROOT/'references/processed/profile_manifest.json').read_text())[0]['family'])
        recovery=json.loads((self.out/'family_target_recovery/manifest.json').read_text())
        self.assertEqual(recovery['status'],'complete');self.assertEqual(recovery['source_kind'],'recovered_table')
        self.assertEqual(source_hash,sha(self.a/'proteins.faa'))

    def test_recovered_domain_only_query_enters_weak_ledger(self):
        self._bind_recovery_run();(self.a/'family.tbl').unlink()
        def recovery_with_domain(fasta,out,root,threads,thresholds):
            self._fake_recovery(fasta,out,root,threads,thresholds)
            profile=json.loads((ROOT/'references/processed/profile_manifest.json').read_text())[0]['profile']
            (out/'family.domtbl').write_text(f'shared - 200 {profile} - 570 1e-30 100 0 1 1 1e-30 1e-2 100 0 1 20 1 20 1 20 0.99 -\n')
        with patch('candidate_screen.verify_completed',return_value=self.dman), patch('candidate_screen.search_discovery',side_effect=recovery_with_domain):
            report=screen(self.a,self.d,self.out)
        rows={r['protein_id']:r for r in self.read('candidate_screen.tsv')}
        self.assertIn('shared',rows)
        self.assertEqual(report['domain_source'],'recovered_table')
        self.assertEqual(rows['shared']['discovery_status'],'weak_or_partial_signal')

    def test_existing_family_table_does_not_rescan(self):
        with patch('candidate_screen.verify_completed',return_value=self.dman), patch('candidate_screen.search_discovery',side_effect=AssertionError('unexpected recovery')):
            report=screen(self.a,self.d,self.out)
        self.assertEqual(report['whole_target_source'],'saved_table')
        self.assertFalse((self.out/'family_target_recovery').exists())

    def test_recovery_model_mismatch_rejects_without_complete_manifest(self):
        self._bind_recovery_run('wrong');(self.a/'family.tbl').unlink()
        with patch('candidate_screen.verify_completed',return_value=self.dman), self.assertRaisesRegex(ValueError,'model_sha256'):
            screen(self.a,self.d,self.out)
        self.assertFalse((self.out/'candidate_screen_manifest.json').exists())

    def test_failed_recovery_is_recorded_failed_not_complete(self):
        self._bind_recovery_run();(self.a/'family.tbl').unlink()
        with patch('candidate_screen.verify_completed',return_value=self.dman), patch('candidate_screen.search_discovery',side_effect=RuntimeError('scan failed')), self.assertRaisesRegex(RuntimeError,'scan failed'):
            screen(self.a,self.d,self.out)
        recovery=json.loads((self.out/'family_target_recovery/manifest.json').read_text())
        self.assertEqual(recovery['status'],'failed')
        self.assertFalse((self.out/'candidate_screen_manifest.json').exists())

    def test_recovery_rejects_source_mutation_during_scan(self):
        self._bind_recovery_run();(self.a/'family.tbl').unlink()
        def mutating_scan(fasta, out, root, threads, thresholds):
            self._fake_recovery(fasta,out,root,threads,thresholds)
            with fasta.open('a') as handle: handle.write('>tampered\nM\n')
        with patch('candidate_screen.verify_completed',return_value=self.dman), patch('candidate_screen.search_discovery',side_effect=mutating_scan), self.assertRaisesRegex(ValueError,'changed during scan'):
            screen(self.a,self.d,self.out)
        recovery=json.loads((self.out/'family_target_recovery/manifest.json').read_text())
        self.assertEqual(recovery['status'],'failed')
        self.assertFalse((self.out/'candidate_screen_manifest.json').exists())

    def test_recovery_rejects_source_family_table_created_during_scan(self):
        self._bind_recovery_run();(self.a/'family.tbl').unlink()
        def creates_source_table(fasta,out,root,threads,thresholds):
            self._fake_recovery(fasta,out,root,threads,thresholds)
            (self.a/'family.tbl').write_text('# unexpected source mutation\n')
        with patch('candidate_screen.verify_completed',return_value=self.dman), patch('candidate_screen.search_discovery',side_effect=creates_source_table), self.assertRaisesRegex(ValueError,'family.tbl appeared'):
            screen(self.a,self.d,self.out)
        recovery=json.loads((self.out/'family_target_recovery/manifest.json').read_text())
        self.assertEqual(recovery['status'],'failed')
        self.assertFalse((self.out/'candidate_screen_manifest.json').exists())
if __name__=='__main__': unittest.main()
