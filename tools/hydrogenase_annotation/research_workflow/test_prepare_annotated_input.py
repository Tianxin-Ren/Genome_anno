import sys
import unittest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_annotated_input import annotated_translation


class AnnotatedSecTests(unittest.TestCase):
    def fixture(self, reverse=False):
        dna = Seq('GTGGCTTGAGGTTAA'); record = SeqRecord(dna.reverse_complement() if reverse else dna, id='c')
        position = 'complement(7..9)' if reverse else '7..9'
        feature = SeqFeature(SimpleLocation(0, 15, strand=-1 if reverse else 1), type='CDS', qualifiers=dict(translation=['MAXG'],
            transl_except=['(pos:'+position+',aa:OTHER)'], note=['location 7..9 selenocysteine, tga no stop']))
        return record, feature

    def test_explicit_source_sec_and_alternative_start(self):
        r, f = self.fixture(); sequence, changes = annotated_translation(r, f)
        self.assertEqual(sequence, 'MAUG'); self.assertEqual(changes[0]['protein_position'], 3)
        self.assertEqual(changes[0]['genomic_positions_1based'], [7, 8, 9])

    def test_reverse_strand_uses_coding_orientation(self):
        r, f = self.fixture(True); sequence, changes = annotated_translation(r, f)
        self.assertEqual(sequence, 'MAUG'); self.assertEqual(changes[0]['genomic_positions_1based'], [9, 8, 7])

    def test_other_without_position_bound_sec_note_remains_unknown(self):
        r, f = self.fixture(); f.qualifiers['note'] = ['unknown translation exception']
        self.assertEqual(annotated_translation(r, f), ('MAXG', []))
        f.qualifiers['note'] = ['location 1..3 selenocysteine']
        self.assertEqual(annotated_translation(r, f), ('MAXG', []))

    def test_sec_amino_acid_exception_is_sufficient_and_preserves_existing_u(self):
        r, f = self.fixture(); f.qualifiers['transl_except'] = ['(pos:7..9,aa:Sec)']; f.qualifiers['note'] = []
        self.assertEqual(annotated_translation(r, f)[0], 'MAUG')
        f.qualifiers['translation'] = ['MAUG']; self.assertEqual(annotated_translation(r, f), ('MAUG', []))
        del f.qualifiers['transl_except']; self.assertEqual(annotated_translation(r, f), ('MAUG', []))

    def test_conflicting_codon_frame_or_amino_acid_cannot_be_recoded(self):
        for value in ['(pos:4..6,aa:Sec)', '(pos:6..8,aa:Sec)', '(pos:20..22,aa:Sec)']:
            r, f = self.fixture(); f.qualifiers['transl_except'] = [value]
            with self.assertRaises(ValueError): annotated_translation(r, f)
        r, f = self.fixture(); f.qualifiers['translation'] = ['MACG']
        with self.assertRaises(ValueError): annotated_translation(r, f)


if __name__ == '__main__': unittest.main()
