import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from candidate_confirmation import confirm, file_sha


def write_table(path, rows):
    fields=sorted({k for r in rows for k in r})
    with path.open('w', newline='') as h:
        w=csv.DictWriter(h, fieldnames=fields, delimiter='\t');w.writeheader();w.writerows(rows)

def fake_extract_empty(*args, **kwargs):
    out = Path(args[1]); out.mkdir(parents=True, exist_ok=True)
    (out / 'features.json').write_text('[]\n')
    return []


class CandidateConfirmationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.root=Path(self.tmp.name)
        self.annotation=self.root/'annotation'; self.screen=self.root/'screen'; self.annotation.mkdir(); self.screen.mkdir()
        write_table(self.annotation/'catalytic_anchors.tsv',[dict(protein_id='anchor')])
        write_table(self.annotation/'genes.tsv',[dict(id='anchor'),dict(id='ga'),dict(id='hp')])
        (self.annotation/'proteins.faa').write_text('>anchor\nMKK\n>ga\nMKK\n>hp\nMKK\n')
        (self.annotation/'run_manifest.json').write_text(json.dumps(dict(status='complete',input_files=[])))
        rows=[dict(protein_id='anchor',family_candidates='["NiFe_related"]',source_sets='[]',discovery_status='existing_anchor',rejected_reason=''),
              dict(protein_id='ga',family_candidates='["FeFe_related"]',source_sets='["verified_GA_domain"]',discovery_status='discovery_only',rejected_reason=''),
              dict(protein_id='hp',family_candidates='["NiFe_or_Nuo_related"]',source_sets='["core_rejected_retained"]',discovery_status='noncatalytic_homolog_favored_retained',rejected_reason='non_catalytic_homolog_wins')]
        write_table(self.screen/'candidate_screen.tsv',rows)
        (self.screen/'candidate_screen_manifest.json').write_text(json.dumps({'status':'complete','annotation':str(self.annotation),'scope':'High-recall homolog queue only. No identity or function promotion','input_sha256':{str(self.annotation/'catalytic_anchors.tsv'):file_sha(self.annotation/'catalytic_anchors.tsv')},'output_sha256':{'candidate_screen.tsv':file_sha(self.screen/'candidate_screen.tsv')}}))
    def tearDown(self): self.tmp.cleanup()
    def test_background_core_rejected_is_not_upgraded(self):
        with patch('active_centers.annotate_active_centers', self.fake_active), patch('noncanonical_center_evidence.annotate_evidence', lambda *a,**k: None), patch('catalytic_confirmation.extract_features', fake_extract_empty):
            report=confirm(self.annotation,self.screen,self.root/'out')
        hp=next(r for r in report['decisions'] if r['candidate_id']=='hp')
        self.assertEqual(hp['sequence_status'],'unconfirmed')
        self.assertIn('core_rejected_non_catalytic_homolog',hp['blocking_reasons'])
    def test_multi_family_is_blocked_without_selection(self):
        with (self.screen/'candidate_screen.tsv').open() as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        rows.append(dict(protein_id='shared',family_candidates='["NiFe_related", "FeFe_related"]',source_sets='[]',discovery_status='discovery_only',rejected_reason=''))
        write_table(self.screen/'candidate_screen.tsv',rows)
        p=json.loads((self.screen/'candidate_screen_manifest.json').read_text());p['output_sha256']['candidate_screen.tsv']=file_sha(self.screen/'candidate_screen.tsv');(self.screen/'candidate_screen_manifest.json').write_text(json.dumps(p))
        with patch('active_centers.annotate_active_centers', self.fake_active), patch('noncanonical_center_evidence.annotate_evidence', lambda *a,**k: None), patch('catalytic_confirmation.extract_features', fake_extract_empty):
            report=confirm(self.annotation,self.screen,self.root/'out')
        shared=next(r for r in report['decisions'] if r['candidate_id']=='shared')
        self.assertEqual(shared['sequence_status'],'unconfirmed');self.assertIn('multi_family_ambiguity',shared['blocking_reasons'])
    def test_source_mismatch_rejected(self):
        other=self.root/'other';other.mkdir();p=json.loads((self.screen/'candidate_screen_manifest.json').read_text());p['annotation']=str(other);(self.screen/'candidate_screen_manifest.json').write_text(json.dumps(p))
        with self.assertRaises(ValueError): confirm(self.annotation,self.screen,self.root/'out')
    def test_missing_bound_input_rejected_before_confirmation(self):
        (self.annotation/'catalytic_anchors.tsv').unlink()
        with self.assertRaisesRegex(ValueError,'Candidate screen input changed'):
            confirm(self.annotation,self.screen,self.root/'out')
        self.assertFalse((self.root/'out').exists())
    def test_confirmed_candidate_cannot_emit_high_or_function(self):
        # Exercise the isolated candidate path with all expensive evidence calls mocked.
        def fake_active(fasta,out,**kwargs):
            out.mkdir(exist_ok=True);(out/'active_center_blast.tsv').write_text('');write_table(out/'active_center_summary.tsv',[dict(protein_id='ga',active_center_status='template_ligand_set_supported',template_span_alignment_status='supported')])
        def fake_extract(*args,**kwargs):
            out=args[1];out.mkdir(parents=True,exist_ok=True);f=dict(protein_id='ga',family='FeFe',eligible=True,eligibility_reasons=[],core_domain_coverage=1,positive_score=100,competitor_score=0,competitor_margin=100,positive_prototype='p',competitor_prototype='',source_active_center={});(out/'features.json').write_text(json.dumps([f]));return [f]
        with patch('active_centers.annotate_active_centers',fake_active), patch('noncanonical_center_evidence.annotate_evidence',lambda *a,**k: None), patch('catalytic_confirmation.extract_features',fake_extract), patch('explain_candidates.load_lineage_review',return_value={}), patch('explain_candidates.sequence_status',return_value='sequence_supported_uncalibrated_core'):
            report=confirm(self.annotation,self.screen,self.root/'out')
        ga=next(r for r in report['decisions'] if r['candidate_id']=='ga')
        self.assertTrue(ga['sequence_status'].startswith('sequence_supported_uncalibrated'))
        self.assertNotEqual(ga['sequence_status'], 'sequence_supported_uncalibrated_core_region')
        self.assertNotEqual(ga['sequence_status'],'calibrated_high');self.assertEqual(ga['function_status'],'unassigned')
        self.assertNotIn('noncanonical',report['evidence']['ga'])

    def test_raw_discovery_signal_stays_unconfirmed_and_unassigned(self):
        write_table(self.annotation/'genes.tsv',[dict(id=p) for p in ['anchor','ga','hp','raw']])
        with (self.annotation/'proteins.faa').open('a') as handle: handle.write('>raw\nMKK\n')
        with (self.screen/'candidate_screen.tsv').open() as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        rows.append(dict(protein_id='raw',family_candidates='["FeFe"]',
                         source_sets='["raw_discovery_signal"]',
                         discovery_status='weak_or_partial_signal',rejected_reason=''))
        write_table(self.screen/'candidate_screen.tsv',rows)
        manifest=json.loads((self.screen/'candidate_screen_manifest.json').read_text())
        manifest['output_sha256']['candidate_screen.tsv']=file_sha(self.screen/'candidate_screen.tsv')
        (self.screen/'candidate_screen_manifest.json').write_text(json.dumps(manifest))
        with patch('active_centers.annotate_active_centers',self.fake_active), \
             patch('noncanonical_center_evidence.annotate_evidence',lambda *a,**k: None), \
             patch('catalytic_confirmation.extract_features',fake_extract_empty):
            report=confirm(self.annotation,self.screen,self.root/'out')
        raw=next(r for r in report['decisions'] if r['candidate_id']=='raw')
        self.assertEqual(raw['sequence_status'],'unconfirmed')
        self.assertEqual(raw['function_status'],'unassigned')
        self.assertTrue(raw['blocking_reasons'])

    def fake_active(self, fasta, out, **kwargs):
        out.mkdir(exist_ok=True); (out/'active_center_blast.tsv').write_text('')
        write_table(out/'active_center_summary.tsv',[dict(protein_id='ga',active_center_status='template_ligand_set_supported',template_span_alignment_status='supported')])


if __name__=='__main__': unittest.main()
