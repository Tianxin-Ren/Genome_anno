"""Contract tests for adding unresolved DNA loci to the catalog."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

from candidate_screen import catalog, dna_only_catalog, sha


def write(path, rows):
    fields=list(rows[0])
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)


class CatalogDnaRescueTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();root=Path(self.tmp.name)
        self.s=root/'screen';self.f=root/'functions';self.a=root/'annotation';self.d=root/'dna';self.s.mkdir();self.f.mkdir();self.a.mkdir();self.d.mkdir()
        self.out=root/'catalog.tsv';self.run=self.a/'run_manifest.json';self.run.write_text('{}')
        write(self.s/'candidate_screen.tsv',[dict(protein_id='protein_1',source_sets='["existing_anchor"]',discovery_status='existing_anchor',family_candidates='["NiFe"]',rejected_reason='')])
        (self.s/'candidate_screen_manifest.json').write_text(json.dumps(dict(annotation=str(self.a),output_sha256={'candidate_screen.tsv':sha(self.s/'candidate_screen.tsv')},input_sha256={str(self.run):sha(self.run)})))
        write(self.f/'functional_annotations.tsv',[dict(candidate_id='protein_1',identity_status='not_confirmed',function_prediction='potential_homolog_function_unresolved',supported_attributes='[]',confidence='unresolved')])
        (self.f/'manifest.json').write_text(json.dumps(dict(functional_annotations_sha256=sha(self.f/'functional_annotations.tsv'),input_sha256={str(self.run):sha(self.run)})))
        self.source=self.a/'input.fna';self.source.write_text('>ctg\nATGAAA\n')
        self.run.write_text(json.dumps({'status':'complete','input_files':[{'path':str(self.source.resolve()),'sha256':sha(self.source)}]}))
        (self.a/'genes.tsv').write_text('id\n');(self.a/'proteins.faa').write_text('');(self.a/'prodigal.gff').write_text('##gff-version 3\n')
        self.candidate=self.d/'dna_rescue_candidates.json';self.candidate.write_text(json.dumps([dict(candidate_id='DNA_locus_1',contig='ctg',start=1,end=6,strand=1,frame=1,translation_tables=[11],source_hsp_ids=['h1'],families=['NiFe'],DNA_sha256='d'*64,translation_sha256='a'*64,translation_has_stop=False,identity_status='unconfirmed',function_prediction='potential_homolog_function_unresolved')]))
        self._write_manifest()
    def tearDown(self): self.tmp.cleanup()
    def _write_manifest(self, status='complete'):
        root=Path(__file__).resolve().parents[1]
        artifacts={'dna_rescue_candidates.json':self.candidate,'dna_hsp_ledger.tsv':self.d/'dna_hsp_ledger.tsv','dna_rescue_translations.faa':self.d/'dna_rescue_translations.faa','reference_queries.faa':self.d/'reference_queries.faa','tblastn_version.txt':self.d/'tblastn_version.txt','tblastn_table11.tsv':self.d/'tblastn_table11.tsv','tblastn_table11.log':self.d/'tblastn_table11.log','tblastn_table4.tsv':self.d/'tblastn_table4.tsv','tblastn_table4.log':self.d/'tblastn_table4.log'}
        for p in artifacts.values(): p.write_text(p.read_text() if p.exists() else '')
        artifacts['dna_hsp_ledger.tsv'].write_text('hsp_id\ttranslation\n'+'h1\tM\n')
        bound_paths={'annotation_run_manifest':self.run,'annotation_genes':self.a/'genes.tsv','annotation_proteins':self.a/'proteins.faa','annotation_prodigal_gff':self.a/'prodigal.gff','screen_manifest':self.s/'candidate_screen_manifest.json','screen_table':self.s/'candidate_screen.tsv','source_nucleotide':self.source,'dna_rescue_script':root/'research_workflow/dna_rescue.py','reference_metadata':root/'references/processed/references.json','reference_fasta':root/'references/processed/references.faa','tblastn_binary':root/'.env/bin/tblastn'}
        manifest=dict(status=status,source_mode='annotation_bound',profile_confirmation='not_performed_empty',annotation=str(self.a.resolve()),screen=str(self.s.resolve()),annotation_run_manifest_sha256=sha(self.run),screen_manifest_sha256=sha(self.s/'candidate_screen_manifest.json'),input=str(self.source),bound_input_paths={k:str(v) for k,v in bound_paths.items()},bound_input_sha256={k:sha(v) for k,v in bound_paths.items()},input_sha256={'source_nucleotide':sha(self.source)},output_sha256={k:sha(v) for k,v in artifacts.items()},translation_tables=[11,4],commands=[['tblastn','-db_gencode','11'],['tblastn','-db_gencode','4']],candidate_count=0 if status!='complete' else 1,candidates=0 if status!='complete' else 1)
        (self.d/'dna_rescue_manifest.json').write_text(json.dumps(manifest))
    def test_bound_dna_candidate_is_appended_as_unresolved_locus(self):
        rows=catalog(self.s,self.f,self.out,dna_rescue_dir=self.d)
        dna=next(r for r in rows if r['candidate_id'].startswith('DNA_'))
        self.assertEqual(dna['candidate_kind'],'dna_locus');self.assertEqual(dna['identity_status'],'unconfirmed');self.assertEqual(dna['function_prediction'],'potential_homolog_function_unresolved')
        with self.out.open() as handle:
            fields=next(csv.DictReader(handle,delimiter='\t'));self.assertEqual(fields['candidate_kind'],'protein')
    def test_changed_dna_output_hash_is_rejected(self):
        self.candidate.write_text(self.candidate.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'output changed'):catalog(self.s,self.f,self.out,dna_rescue_dir=self.d)
    def test_other_annotation_is_rejected(self):
        manifest=json.loads((self.d/'dna_rescue_manifest.json').read_text());manifest['annotation']=str(self.a/'other');(self.d/'dna_rescue_manifest.json').write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError,'ownership'):catalog(self.s,self.f,self.out,dna_rescue_dir=self.d)

    def test_noncomplete_manifest_is_not_treated_as_zero_hit_without_binding(self):
        manifest=json.loads((self.d/'dna_rescue_manifest.json').read_text());manifest.update(status='not_applicable',reason='input_mode_is_transcript',annotation=str(self.a/'other'),output_sha256={},candidate_count=0);(self.d/'dna_rescue_manifest.json').write_text(json.dumps(manifest))
        for p in self.d.glob('*.tsv'): p.unlink()
        for p in self.d.glob('*.faa'): p.unlink()
        self.candidate.unlink()
        with self.assertRaisesRegex(ValueError,'another annotation|ownership'):catalog(self.s,self.f,self.out,dna_rescue_dir=self.d)

    def test_bad_coordinate_or_translation_table_is_rejected(self):
        data=json.loads(self.candidate.read_text());data[0]['strand']=-1;self.candidate.write_text(json.dumps(data));self._write_manifest()
        with self.assertRaisesRegex(ValueError,'locus evidence'):catalog(self.s,self.f,self.out,dna_rescue_dir=self.d)
        data[0]['strand']=1;data[0]['translation_tables']=[False];self.candidate.write_text(json.dumps(data));self._write_manifest()
        with self.assertRaisesRegex(ValueError,'translation tables'):catalog(self.s,self.f,self.out,dna_rescue_dir=self.d)

    def test_raw_dna_catalog_cannot_be_consumed_as_annotation_bound(self):
        root=Path(__file__).resolve().parents[1]
        bound={'raw_nucleotide':self.source,'dna_rescue_script':root/'research_workflow/dna_rescue.py','reference_metadata':root/'references/processed/references.json','reference_fasta':root/'references/processed/references.faa','tblastn_binary':root/'.env/bin/tblastn'}
        manifest=json.loads((self.d/'dna_rescue_manifest.json').read_text());manifest.update(source_mode='raw_nucleotide_only',annotation=None,screen=None,input=str(self.source.resolve()),bound_input_paths={k:str(v) for k,v in bound.items()},bound_input_sha256={k:sha(v) for k,v in bound.items()});(self.d/'dna_rescue_manifest.json').write_text(json.dumps(manifest))
        rows=dna_only_catalog(self.d,self.out,self.source)
        self.assertEqual(rows[0]['candidate_kind'],'dna_locus')
        with self.assertRaisesRegex(ValueError,'Raw DNA rescue cannot be consumed|ownership mismatch'):
            catalog(self.s,self.f,self.tmp.name+'/ordinary.tsv',dna_rescue_dir=self.d)
    def test_upgraded_or_duplicate_dna_candidate_is_rejected(self):
        data=json.loads(self.candidate.read_text());data[0]['function_prediction']='hydrogenase_family_supported_function_unresolved';self.candidate.write_text(json.dumps(data));self._write_manifest()
        with self.assertRaisesRegex(ValueError,'promoted|upgraded'):catalog(self.s,self.f,self.out,dna_rescue_dir=self.d)
        data[0]['function_prediction']='potential_homolog_function_unresolved';data[0]['candidate_id']='protein_1';self.candidate.write_text(json.dumps(data));self._write_manifest()
        with self.assertRaisesRegex(ValueError,'collision'):catalog(self.s,self.f,self.out,dna_rescue_dir=self.d)


if __name__=='__main__': unittest.main()
