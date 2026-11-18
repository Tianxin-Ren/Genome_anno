"""Actual native HSPs and counterfactual ownership/coordinate guards."""
import copy
import csv
import json
from pathlib import Path
import unittest
import shutil
import tempfile
from Bio import SeqIO
from annotate import load_blast
from complementary_core_context import attributes, evaluate, prepare

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'experiments/round219_SoHydAB_challenge/runs/native_encoded/annotation'


class ComplementaryContextTests(unittest.TestCase):
    def setUp(self):
        self.models={m['reference_id']:m for m in json.loads((ROOT/'references/processed/active_centers/models.json').read_text()) if m['family']=='FeFe'}
        self.params=json.loads((ROOT/'references/processed/active_centers/manifest.json').read_text())['parameters']
        self.thresholds=json.loads((RUN/'run_manifest.json').read_text())['thresholds']
        with (RUN/'genes.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.genes={r['id']:{**r,**{k:int(r[k]) for k in ('start','end','strand','contig_length')}} for r in rows}
        with (RUN/'proteins.faa').open() as f:self.sequences={r.id:str(r.seq) for r in SeqIO.parse(f,'fasta')}
        self.hits=load_blast(RUN/'active_center_blast.tsv',dict(blast_evalue=self.params['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0))
        with (RUN/'active_center_summary.tsv').open() as f:self.summaries={r['protein_id']:r for r in csv.DictReader(f,delimiter='\t')}
        self.anchors=[dict(protein_id='q001',catalytic_family='FeFe')]

    def run_case(self):
        return evaluate(self.anchors,self.genes,self.sequences,self.models,self.hits,self.summaries,self.params,self.thresholds)

    def clone(self,pid,new,offset):
        self.genes[new]={**self.genes[pid], 'id':new,'start':self.genes[pid]['start']+offset,'end':self.genes[pid]['end']+offset}
        self.sequences[new]=self.sequences[pid]
        self.hits.extend(dict(h,qseqid=new) for h in list(self.hits) if h['qseqid']==pid)

    def test_actual_partial_context_does_not_pass_the_last_site_guard(self):
        result=self.run_case()['q001'];candidate,=result['candidates']
        self.assertEqual(result['status'],'conditional')
        self.assertEqual(candidate['new_reference_positions'],list(range(524,573)))
        self.assertAlmostEqual(candidate['metrics']['min_required_context_coverage'],23/31)
        self.assertFalse(candidate['metrics']['context_alignment_supported'])
        self.assertIn('required_site_context_incomplete',candidate['blockers'])
        self.assertEqual({r['protein_id'] for p,r in candidate['reference_position_owners'].items() if int(p) in {300,355,499,503}},{'q001'})
        row,=attributes(result)
        self.assertEqual(row['status'],'conditional_partial_context')
        self.assertEqual(result['physical_assembly'],'unknown')
        self.assertEqual(result['query_activity'],'not_tested')

    def test_deleted_tail_withdraws_only_the_context_candidate(self):
        del self.genes['q002'];del self.sequences['q002']
        self.hits=[h for h in self.hits if h['qseqid']!='q002']
        self.assertEqual(self.run_case(),{})
        self.assertEqual(attributes({}),[])

    def test_rename_and_record_order_do_not_change_the_explanation(self):
        expected=self.run_case()
        names={'q001':'opaque_z','q002':'opaque_a'}
        self.genes={names[k]:dict(v,id=names[k]) for k,v in reversed(list(self.genes.items()))}
        self.sequences={names[k]:v for k,v in reversed(list(self.sequences.items()))}
        self.hits=[dict(h,qseqid=names[h['qseqid']]) for h in reversed(self.hits)]
        self.anchors=[dict(protein_id='opaque_z',catalytic_family='FeFe')]
        self.summaries={names[k]:dict(v,protein_id=names[k]) for k,v in self.summaries.items()}
        def restore(x):
            if isinstance(x,dict):return {restore(k):restore(v) for k,v in x.items()}
            if isinstance(x,list):return [restore(v) for v in x]
            return {v:k for k,v in names.items()}.get(x,x) if isinstance(x,str) else x
        self.assertEqual(restore(self.run_case()),expected)

    def test_circular_origin_preserves_context_metrics(self):
        expected=self.run_case()['q001']
        for gene in self.genes.values():
            for k in ('start','end'):gene[k]=(gene[k]-4069055-1)%gene['contig_length']+1
        observed=self.run_case()['q001']
        self.assertEqual(observed,expected)

    def test_far_tail_and_opposite_strand_have_no_local_assignment(self):
        for mutate in ('far','opposite','missing_coordinates','partial','wrong_bin'):
            with self.subTest(mutate=mutate):
                saved=copy.deepcopy(self.genes)
                if mutate=='far':
                    self.genes['q002']['start']+=200000;self.genes['q002']['end']+=200000
                if mutate=='opposite':self.genes['q002']['strand']=-1
                if mutate=='missing_coordinates':self.genes['q002']['coordinate_evidence']='unavailable_translated_transcript'
                if mutate=='partial':self.genes['q002']['partial']='1'
                if mutate=='wrong_bin':self.genes['q002']['bin_id']='other'
                result=self.run_case()['q001']
                self.assertEqual(result['status'],'unresolved')
                self.assertEqual(result['selected_members'],[])
                self.assertEqual(result['candidates'][0]['association_status'],'rejected')
                self.genes=saved

    def test_duplicate_tail_loci_remain_ambiguous(self):
        self.clone('q002','copy_tail',1000)
        result=self.run_case()['q001']
        self.assertEqual(result['status'],'ambiguous')
        self.assertEqual(result['selected_members'],[])
        self.assertEqual(len(result['candidates']),2)

    def test_remote_paralog_does_not_erase_unique_local_linkage(self):
        self.clone('q002','remote_tail',200000)
        result=self.run_case()['q001']
        self.assertEqual(result['status'],'conditional')
        self.assertEqual(result['selected_members'][1]['protein_id'],'q002')

    def test_alternative_core_owners_are_recorded_without_copy_capacity_rule(self):
        self.clone('q001','second_core',-1500)
        self.anchors.append(dict(protein_id='second_core',catalytic_family='FeFe'))
        self.summaries['second_core']=dict(self.summaries['q001'],protein_id='second_core')
        results=self.run_case()
        self.assertEqual(set(results),{'q001','second_core'})
        for result in results.values():
            self.assertEqual(result['status'],'ambiguous')
            self.assertEqual(result['selected_members'],[])
            self.assertEqual(len(result['candidates'][0]['alternative_anchors']),1)

    def test_duplicate_identical_search_row_does_not_invent_another_locus(self):
        self.hits.extend(copy.deepcopy(self.hits))
        self.assertEqual(self.run_case()['q001']['status'],'conditional')

    def test_unrelated_protein_and_metadata_do_not_change_target(self):
        expected=self.run_case()
        self.genes['unrelated']=dict(self.genes['q002'],id='unrelated')
        self.sequences['unrelated']='M'*100
        self.assertEqual(self.run_case(),expected)

    def test_template_name_is_not_an_algorithm_branch(self):
        old=next(iter(self.models));new='reviewed_template_renamed'
        self.models={new:dict(self.models[old],reference_id=new)}
        self.hits=[dict(h,sseqid=new if h['sseqid']==old else h['sseqid']) for h in self.hits]
        self.summaries={k:dict(v,reference_id=new if v['reference_id']==old else v['reference_id']) for k,v in self.summaries.items()}
        self.assertEqual(self.run_case()['q001']['status'],'conditional')

    def test_altered_alignment_and_sequence_length_are_rejected(self):
        target=next(h for h in self.hits if h['qseqid']=='q002')
        original=target['qseq'];target['qseq']='X'+original[1:]
        with self.assertRaises(ValueError):self.run_case()
        target['qseq']=original;target['qlen']+=1
        with self.assertRaises(ValueError):self.run_case()

    def test_saved_producer_tables_and_normalized_faa_are_bound(self):
        with tempfile.TemporaryDirectory() as name:
            folder=Path(name);ann=folder/'annotation';ann.mkdir()
            for n in ('active_center_blast.tsv','active_center_summary.tsv','genes.tsv','proteins.faa','run_manifest.json'):
                shutil.copyfile(RUN/n,ann/n)
            supplied={pid:{k:str(v) for k,v in g.items()} for pid,g in self.genes.items()}
            result,_=prepare(ann,folder/'good',self.anchors,supplied,self.sequences,ROOT)
            self.assertEqual(result['q001']['status'],'conditional')
            for artifact in ('active_center_summary.tsv','genes.tsv','proteins.faa'):
                path=ann/artifact;original=path.read_bytes();path.write_bytes(original+b'\n')
                with self.subTest(artifact=artifact),self.assertRaises(ValueError):
                    prepare(ann,folder/'bad',self.anchors,supplied,self.sequences,ROOT)
                path.write_bytes(original)


if __name__=='__main__':unittest.main()
