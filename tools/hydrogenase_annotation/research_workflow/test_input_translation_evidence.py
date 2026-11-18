import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from Bio.Seq import Seq
sys.path[:0] = [str(Path(__file__).resolve().parents[1]/'scripts'), str(Path(__file__).resolve().parent)]
from input_translation_evidence import prodigal_tables, audit, sha


class TranslationModelEvidenceTests(unittest.TestCase):
    def test_model_is_associated_with_actual_contig_not_global_default(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/'x.gff'
            p.write_text('# Sequence Data: a\n# Model Data: version=x;transl_table=4;uses_sd=1\na\tProdigal\tCDS\t1\t12\t.\t+\t0\tID=a\n# Sequence Data: b\n# Model Data: version=x;transl_table=11;uses_sd=1\nb\tProdigal\tCDS\t1\t12\t.\t+\t0\tID=b\n')
            self.assertEqual(prodigal_tables(p), ({'a': 4, 'b': 11}, {'a': 1, 'b': 1}))

    def test_missing_models_remain_unknown_and_conflicting_models_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/'x.gff'
            p.write_text('a\tProdigal\tCDS\t1\t12\t.\t+\t0\tID=a\n')
            self.assertEqual(prodigal_tables(p), ({'a': None}, {'a': 1}))
            p.write_text('# Model Data: version=x;transl_table=4;uses_sd=1\na\tProdigal\tCDS\t1\t12\t.\t+\t0\tID=a\n# Model Data: version=x;transl_table=11;uses_sd=1\na\tProdigal\tCDS\t20\t40\t.\t+\t0\tID=b\n')
            with self.assertRaises(ValueError): prodigal_tables(p)

    def test_missing_headers_do_not_block_other_annotation_layers(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory); gff = p/'provided.gff'
            gff.write_text('a\tProdigal\tCDS\t1\t12\t.\t+\t0\tID=a\n')
            (p/'run_manifest.json').write_text(json.dumps(dict(translation_table=11,input_files=[dict(path=str(gff),sha256=sha(gff))])))
            r, _ = audit(p, [], [])
            self.assertEqual(r['status'], 'Prodigal_model_headers_incomplete')
            self.assertEqual(r['effective_tables_by_contig'], {'a': None})
            self.assertIsNone(r['requested_effective_disagreement'])

    def test_tga_called_w_is_model_sensitive_without_being_relabelled_sec(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory); dna = p/'input.fna'; dna.write_text('>a\nATGGCTTGAGGTTAA\n')
            (p/'prodigal.gff').write_text('# Model Data: version=x;transl_table=4;uses_sd=1\na\tProdigal\tCDS\t1\t15\t.\t+\t0\tID=q\n')
            (p/'run_manifest.json').write_text(json.dumps(dict(translation_table=11, input_files=[dict(path=str(dna), sha256=sha(dna))])))
            gene=dict(id='q',contig='a',start='1',end='15',strand='1')
            site=dict(protein_id='q',status='substituted',observed_residue='W',expected_residue='CU',query_position='3',reference_position='20',residue_role='direct_metal_ligand')
            r, bindings = audit(p, [gene], [site]); self.assertTrue(r['requested_effective_disagreement'])
            self.assertEqual(r['sensitive_sites'][0]['genomic_codon'], 'TGA')
            self.assertEqual(r['sensitive_sites'][0]['observed_called_residue'], 'W')
            self.assertIn('no U is assigned', r['sensitive_sites'][0]['interpretation'])
            self.assertIn(str(dna), bindings)
            context, _ = audit(p, [gene], [dict(site, residue_role='L1_L2_sequence_context')])
            self.assertEqual(context['sensitive_sites'], [])
            compressed = p/'input.fna.gz'
            compressed.write_bytes(gzip.compress(dna.read_bytes()))
            manifest = p/'run_manifest.json'
            original = manifest.read_text()
            manifest.write_text(json.dumps(dict(translation_table=11, input_files=[
                dict(path=str(compressed), sha256=sha(compressed))])))
            compressed_result, compressed_bindings = audit(p, [gene], [site])
            self.assertEqual(compressed_result, r)
            self.assertEqual(compressed_bindings[str(compressed)], sha(compressed))
            manifest.write_text(original)
            dna.write_text('>a\nATGGCTTGGGGTTAA\n')
            with self.assertRaises(ValueError): audit(p, [gene], [site])

    def test_annotated_gff_without_prodigal_header_has_no_invented_code(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory); gff = p/'genes.gff'; gff.write_text('a\texplicit_source_CDS\tCDS\t1\t15\t.\t+\t0\tID=q\n')
            (p/'run_manifest.json').write_text(json.dumps(dict(translation_table=None,input_files=[dict(path=str(gff),sha256=sha(gff))])))
            r, _ = audit(p, [], []); self.assertEqual(r['status'], 'GFF_has_no_Prodigal_model_header')

    def test_circular_codon_observations_are_rotation_invariant_on_both_strands(self):
        coding = 'ATGGCTTGAGGTTAA'
        site = dict(protein_id='q', status='substituted', observed_residue='W', expected_residue='CU',
                    query_position='3', reference_position='20', residue_role='direct_metal_ligand')
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory); source = p/'input.fna'
            for strand in (1, -1):
                native = (coding if strand == 1 else str(Seq(coding).reverse_complement())) + 'C'*15
                baseline = None
                for rotation in (0, 7, 10):
                    with self.subTest(strand=strand, rotation=rotation):
                        dna = native[rotation:] + native[:rotation]
                        source.write_text('>a\n'+dna+'\n')
                        start = (-rotation) % len(dna) + 1; end = start + len(coding) - 1
                        (p/'prodigal.gff').write_text('# Model Data: version=x;transl_table=4;uses_sd=1\n'
                            f'a\tProdigal\tCDS\t{start}\t{end}\t.\t{"+" if strand == 1 else "-"}\t0\tID=q\n')
                        (p/'run_manifest.json').write_text(json.dumps(dict(translation_table=11,
                            input_files=[dict(path=str(source), sha256=sha(source))])))
                        gene = dict(id='q', contig='a', start=start, end=end, strand=strand,
                                    topology='circular', contig_length=len(dna))
                        result, _ = audit(p, [gene], [site])
                        self.assertEqual(len(result['sensitive_sites']), 1)
                        self.assertEqual(result['sensitive_sites'][0]['genomic_codon'], 'TGA')
                        self.assertEqual(result['sensitive_sites'][0]['observed_called_residue'], 'W')
                        self.assertIn('no U is assigned', result['sensitive_sites'][0]['interpretation'])
                        if baseline is None: baseline = result['sensitive_sites']
                        self.assertEqual(result['sensitive_sites'], baseline)

    def test_invalid_coordinate_or_topology_cannot_create_a_codon_observation(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory); source = p/'input.fna'; source.write_text('>a\nATGGCTTGAGGTTAA\n')
            (p/'prodigal.gff').write_text('# Model Data: version=x;transl_table=4;uses_sd=1\n'
                                       'a\tProdigal\tCDS\t1\t15\t.\t+\t0\tID=q\n')
            (p/'run_manifest.json').write_text(json.dumps(dict(translation_table=11,
                input_files=[dict(path=str(source), sha256=sha(source))])))
            gene = dict(id='q', contig='a', start=1, end=15, strand=1, topology='circular', contig_length=15)
            site = dict(protein_id='q', status='substituted', observed_residue='W', expected_residue='CU',
                        query_position='3', reference_position='20', residue_role='direct_metal_ligand')
            for changes in (dict(contig_length=16), dict(strand=0), dict(end=30),
                            dict(end=30, topology='linear')):
                with self.subTest(changes=changes):
                    result, _ = audit(p, [dict(gene, **changes)], [site])
                    self.assertEqual(result['sensitive_sites'], [])


if __name__ == '__main__': unittest.main()
