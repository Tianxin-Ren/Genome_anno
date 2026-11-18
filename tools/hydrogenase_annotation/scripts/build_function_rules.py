#!/usr/bin/env python3
"""Compile the curated evidence catalog into offline, sequence-bound function rules."""
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent/'outputs/01a18a6b-afa7-7331-abea-680ebf4e2278/evidence_six_20261112'
DEST = ROOT/'references/processed/function_rules'

MEMBER_POLICIES = {
 'HdrB':dict(profiles=['PF02754'],motif='C.{28,45}CC.{28,45}C..C',motif_count=2,min_conserved_ligands=.875),
 'HdrC':dict(profiles=['PF13183'],min_identity=25,min_profile_coverage=.8,min_reference_domain_coverage=.75,
            motif='C..C..C...C',motif_count=2,min_conserved_ligands=.875),
 'EchA':dict(profiles=['PF00361'],min_identity=25,min_profile_coverage=.8,min_reference_domain_coverage=.7),
 'EchB':dict(profiles=['PF00146'],min_identity=25,min_profile_coverage=.8,min_reference_domain_coverage=.7),
 'EchD':dict(profiles=['PF00329'],min_profile_coverage=.8,profile_only_catalytic_fusion=True),
 'EchF':dict(profiles=['PF00037','PF12838','PF13187','PF13237','PF14697'],min_identity=30,min_profile_coverage=.8,
            min_reference_domain_coverage=.7,motif='C..C..C...C',motif_count=2,min_conserved_ligands=.875),
}

# Core module descriptions are functional requirements, not intersections of gene names.
CORES = {
 'Ni01': ['NiFe_H2_core','FeS_small_subunit','cytochrome_c3_partner'],
 'Ni02': ['NiFe_H2_core','FeS_small_subunit','membrane_electron_exit'],
 'Ni03': ['NiFe_H2_core','FeS_small_subunit','quinone_linked_membrane_module'],
 'Ni04': ['H2_core','quinone_linked_electron_transfer_module'],
 'Ni05': ['NiFe_H2_core','oxygen_tolerant_reference_small_subunit','quinone_linked_membrane_module'],
 'Ni06': ['Hhy_like_NiFe_core','Hhy_small_subunit'],
 'Ni07': ['Vht_Vho_like_NiFe_core','FeS_small_subunit','methanophenazine_linked_membrane_module'],
 'Ni08': ['Hup_like_NiFe_core','Hup_small_subunit'],
 'Ni09': ['Huc_like_NiFe_core','HucS','HucM_quinone_delivery_scaffold'],
 'Ni10': ['H2_sensing_NiFe_core','sensor_small_subunit','signal_transduction_context'],
 'Ni11': ['NiFe_H2_core','FeS_small_subunit','FrhB_F420_reductase'],
 'Ni12': ['NiFe_H2_core','FeS_relay','NADP_compatible_pyridine_module'],
 'Ni12R': ['NiFe_H2_core','FeS_relay','NADP_compatible_pyridine_module'],
 'Ni13': ['Mvh_AGD_H2_module','HdrABC_heterodisulfide_module'],
 'Ni14': ['NiFe_H2_core','FeS_relay','NAD_compatible_pyridine_module'],
 'Ni15': ['NiFe_H2_core','formate_dehydrogenase','FeS_relay','FHL_membrane_module'],
 'Ni16': ['NiFe_H2_core','formate_dehydrogenase','membrane_energy_coupling_module'],
 'Ni17': ['NiFe_H2_core','CO_dehydrogenase','FeS_relay','Coo_membrane_module'],
 'Ni18': ['NiFe_H2_core','Fd_compatible_relay','membrane_energy_coupling_module'],
 'Ni19': ['Ech_NiFe_core','Fd_compatible_relay','Ech_membrane_module'],
 'Ni20': ['energy_converting_NiFe_core','Fd_compatible_relay','membrane_module'],
 'Ni23': ['Ech_NiFe_core','Fd_compatible_relay','Ech_membrane_module'],
 'Ni22': ['Hyf_NiFe_core','formate_dehydrogenase','FeS_relay','FHL2_membrane_module'],
 'Ni24': ['NiFe_H2_core','FeS_small_subunit','pyridine_module','additional_Fd_compatible_relay'],
 'Fe01': ['FeFe_H_cluster_core','intrinsic_or_external_Fd_compatible_relay'],
 'Fe02': ['HYDA_like_FeFe_core','reduced_Fd_compatible_electron_input'],
 'Fe03': ['FeFe_H_cluster_core','HydB_pyridine_module','HydC_FeS_relay'],
 'Fe04': ['FeFe_H_cluster_core','HydB_pyridine_module','HydC_FeS_relay'],
 'Fe05': ['FeFe_H_cluster_core','NADH_compatible_pyridine_module','nonbifurcating_reference_architecture'],
 'Fe06': ['FeFe_H_cluster_core','formate_dehydrogenase','FeS_relay'],
 'Fe07': ['FeFe_H_cluster_core','NADP_compatible_pyridine_module','Fd_compatible_relay'],
 'Fe08': ['FeFe_H_cluster_core','NADP_compatible_pyridine_module','Fd_independent_reference_architecture'],
 'H01': ['Hmd_FeGP_core','H4MPT_compatible_substrate_binding'],
 'Ni21': ['H2_core','external_electron_transfer_partner','external_NADP_reductase'],
}
QUALIFIERS = {
 'Ni01':['native_cytochrome_interaction'], 'Ni02':['downstream_terminal_acceptor'],
 'Ni03':['native_MK_DMK_species','in_vivo_net_direction'],
 'Ni04':['quinol_driven_reverse_flux','driving_force'],
 'Ni05':['measured_oxygen_tolerance'], 'Ni06':['atmospheric_H2_affinity','survival_benefit'],
 'Ni07':['native_methanophenazine_turnover'], 'Ni08':['nitrogen_fixation_dependent_H2_recycling'],
 'Ni09':['atmospheric_H2_affinity','native_long_chain_MK_turnover'],
 'Ni10':['H2_responsive_gene_regulation'], 'Ni11':['F420_substrate_selectivity'],
 'Ni12':['NADP_selectivity','in_vivo_net_direction'], 'Ni12R':['NADP_selectivity','in_vivo_net_direction'],
 'Ni13':['obligatory_Fd_heterodisulfide_bifurcation'], 'Ni14':['NAD_selectivity'],
 'Ni15':['in_vivo_net_direction'], 'Ni16':['ion_translocation','ATP_coupling'],
 'Ni17':['CO_dependent_in_vivo_flux'], 'Ni18':['Na_selectivity','ion_translocation'],
 'Ni19':['acetoclastic_methanogenesis_context'], 'Ni20':['autotrophic_methanogenesis_context','in_vivo_net_direction'],
 'Ni23':['ion_identity','in_vivo_net_direction'], 'Ni22':['formate_dependent_in_vivo_flux'],
 'Ni24':['NAD_selectivity','obligatory_Fd_bifurcation'],
 'Fe01':['physiological_Fd_identity','in_vivo_net_direction'], 'Fe02':['photosynthetic_electron_supply'],
 'Fe03':['obligatory_NADH_Fd_confurcation'],
 'Fe04':['obligatory_NAD_Fd_bifurcation'], 'Fe05':['Fd_nonrequirement','NADH_selectivity'],
 'Fe06':['in_vivo_net_direction'], 'Fe07':['NADP_selectivity','obligatory_Fd_bifurcation'],
 'Fe08':['direct_NADP_reduction','Fd_nonrequirement'], 'H01':['FeGP_occupancy'],
 'Ni21':['external_partner_compatibility','native_physiological_relevance'],
}


def table(name):
    with (SOURCE/name).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def compile_rules():
    evidence=table('complete_evidence.tsv'); seqs={r['序列ID']:r for r in table('sequences.tsv')}
    links=defaultdict(dict)
    for r in table('evidence_sequence_links.tsv'):
        links[r['体系ID']][r['亚基或基因名']]=r
    systems={r['体系ID']:r for r in table('systems_genomes.tsv')}
    sites=defaultdict(list)
    for r in table('catalytic_site_details.tsv'):
        if r['residue_role']=='direct_metal_ligand' and r['status']=='present' and r['query_position'].isdigit():
            sites[r['序列ID']].append(dict(position=int(r['query_position']),residue=r['observed_residue']))
    domains=defaultdict(list)
    for r in table('protein_domain_architecture.tsv'):
        if r['supported'].lower()=='true':
            domains[r['序列ID']].append(dict(accession=r['accession'],segments=json.loads(r['query_segments'])))
    fragments=table('construct_fragments.tsv') if (SOURCE/'construct_fragments.tsv').exists() else []
    references={}; classes=[]
    for row in table('completion_summary.tsv'):
        fid=row['功能ID']; assert fid in CORES and fid in QUALIFIERS
        arches=[]
        for sid in sorted({e['体系ID'] for e in evidence if e['功能ID']==fid}):
            members=[]
            for role,link in links[sid].items():
                qid=link['序列ID']; s=seqs[qid]; kind=link['组件类别']
                assert not (sid=='FeFe_Mthermo_HydABCDEF' and role in {'HydE','HydF'} and kind=='成熟辅助因子'), 'Rebuild evidence catalog after role fix'
                reference=dict(sequence=s['完整CDS蛋白序列'],sha256=s['蛋白SHA256'],
                    accession=s['核苷酸登录号'],locus_tag=s['locus_tag'],organism=s['微生物'],
                    start=int(s['起点_1based']),end=int(s['终点_1based']),strand=s['链方向'],
                    domains=domains[qid],metal_sites=sorted({(x['position'],x['residue']) for x in sites[qid]}),
                    reference_domain_status=s.get('结构域检测状态',''),reference_center_status=s.get('催化中心检测状态',''))
                references[qid]=reference
                if role in MEMBER_POLICIES and MEMBER_POLICIES[role].get('motif'):
                    pol=MEMBER_POLICIES[role]
                    matches=list(re.finditer(pol['motif'].replace('C','(C)'),reference['sequence']))
                    if len(matches)==pol['motif_count']:
                        reference.setdefault('role_coordination_sites',{})[role]=[match.start(i)+1 for match in matches for i in range(1,len(match.groups())+1)]
                requirement='maturation' if kind=='成熟辅助因子' else 'optional' if kind in {'调控辅助基因','生理背景基因'} else 'required'
                if fid in {'Fe01','Fe02'} and kind!='催化亚基': requirement='optional' if requirement!='maturation' else requirement
                if role in {'Fd1','Fd2','Fd3','HydD_candidate','HyfR','HydN'}: requirement='optional'
                if fid=='Ni14' and role=='HoxE':requirement='optional'
                if fid=='Fe07' and role=='FdhA':requirement='optional'
                if sid=='NiFe_Tonn_Frh_external' and role in {'Fdh3B','TON_0702','TON_1376'}:requirement='alternative_external_partner'
                if sid=='FeFe_CbA5H_BMR_external' and role=='BMR_CYP102A1':
                    fragment=next(x for x in fragments if x['片段ID']=='FRAG_BG04_163_BMR_473_1049')
                    fragment_id=fragment['片段ID']; start,end=int(fragment['片段起点_1based']),int(fragment['片段终点_1based'])
                    references[fragment_id]={**reference,'sequence':reference['sequence'][start-1:end],
                        'sha256':hashlib.sha256(reference['sequence'][start-1:end].encode()).hexdigest(),
                        'domains':[dict(accession=d['accession'],segments=[[max(a,start)-start+1,min(b,end)-start+1] for a,b in d['segments'] if max(a,start)<=min(b,end)]) for d in reference['domains'] if any(max(a,start)<=min(b,end) for a,b in d['segments'])],
                        'metal_sites':[],'construct_parent':qid,'construct_range':[start,end]}
                    qid=fragment_id
                members.append(dict(role=role,reference=qid,category=kind,requirement=requirement))
            if fid=='Ni10':
                for m in members:
                    if m['category']=='调控辅助基因':m['requirement']='required'
            selected=[e for e in evidence if e['功能ID']==fid and e['体系ID']==sid]
            candidate_only=sid in {'FeFe_Ssphaeroides_KIAC_HydABC','FeFe_Swolfei_Hyd2'}
            archives=dict(id=fid+'__'+sid,system=sid,center=systems[sid]['活性中心'],members=members,
                support_mode='candidate_only' if candidate_only else 'reference_architecture_prediction',
                association='external' if fid=='Ni21' else 'reference_layout',
                citations=[dict(evidence_id=e['证据ID'],citation=e['DOI或文献键'],grade=e['统一证据等级'],measured=e['实际测得'],limitations=e['限制及未测部分']) for e in selected],
                limitations=systems[sid]['范围与限制'])
            assert any(m['category']=='催化亚基' and m['requirement']=='required' for m in members)
            arches.append(archives)
        classes.append(dict(id=fid,label=row['功能类别'],common_core=CORES[fid],architectures=arches,
            maturation=dict(NiFe=['hypA','hypB','hypC','hypD','hypE','hypF'],FeFe=['hydE','hydF','hydG'],Fe=['hcg_maturation_context']),
            exclusions=['unconfirmed_catalytic_identity','incompatible_catalytic_ligand','required_member_incomplete',
                'cross_bin_pairing','overlapping_roles_on_same_query_region','stronger_incompatible_reference_architecture'],
            unresolved_qualifiers=QUALIFIERS[fid],
            interpretation='Sequence-based functional prediction; membership, mechanism and physiology have separate support. Missing maturation candidates do not veto the reaction prediction.'))
    assert len(classes)==len(CORES) and {c['id'] for c in classes}==set(CORES)
    used={m['reference'] for c in classes for a in c['architectures'] for m in a['members']}
    references={k:v for k,v in references.items() if k in used}
    fdh_model=json.loads((ROOT/'research_workflow/fdh_role_model.json').read_text())
    role_counters=[];counter_sources={}
    for spec in fdh_model['source_counters']:
        path=ROOT/spec['fasta'];sequence=''.join(line.strip() for line in path.read_text().splitlines() if not line.startswith('>'))
        assert hashlib.sha256(sequence.encode()).hexdigest()==spec['sequence_sha256']
        rid='COUNTER_FDH__'+spec['reference_id'];role_counters.append(rid);counter_sources[spec['fasta']]=sha(path)
        references[rid]=dict(sequence=sequence,sha256=spec['sequence_sha256'],domains=[],metal_sites=[],counter_class=spec['role_class'])
    donor_bank=ROOT/'research_workflow/mvh_HdrA_references.faa';donors={};key=None;parts=[]
    for line in donor_bank.read_text().splitlines():
        if line.startswith('>'):
            if key:donors[key]=''.join(parts)
            key=line[1:].split()[0];parts=[]
        else:parts.append(line.strip())
    if key:donors[key]=''.join(parts)
    donor_ids={}
    for name,label in [('H2_linked_HdrA_reference','H2'),('F420_Fd_HdrA2_core_reference','F420_Fd_alternative')]:
        sequence=donors[name];rid='DONOR__'+name;donor_ids[label]=rid
        references[rid]=dict(sequence=sequence,sha256=hashlib.sha256(sequence.encode()).hexdigest(),domains=[],metal_sites=[],comparison_class=label)
    counter_sources[str(donor_bank.relative_to(ROOT))]=sha(donor_bank)
    assert all(hashlib.sha256(r['sequence'].encode()).hexdigest()==r['sha256'] for r in references.values())
    DEST.mkdir(parents=True,exist_ok=True)
    fasta=DEST/'members.faa'
    fasta.write_text(''.join('>'+qid+'\n'+r['sequence'].replace('U','X')+'\n' for qid,r in sorted(references.items())))
    payload=dict(version='function-rules-34-v4',classes=classes,references=references,role_counters=dict(Fdh_catalytic=role_counters),counter_source_sha256=counter_sources,
        member_policies=MEMBER_POLICIES,HdrA_donor_comparison=donor_ids,
        modular_associations={'Ni13':{'local_groups':[['MvhA','MvhG'],['HdrB','HdrC']],
            'remote_relay_requires':'nonoverlapping_HdrA_MvhD_fusion_or_same_strand_local_pair; unique_Mvh_core; no_strong_alternative_donor',
            'HdrA_MvhD_max_gap_bp':2500,'identity_margin_hint':5.0,'alternative_bitscore_ratio':1.10,
            'donor_scope':'Near-tied comparable HdrA references limit donor specificity, not the complete-module coupling prediction. A clear alternative identity and bitscore advantage remains a counter.'},
            'Ech':{'local_groups':[['EchA','EchB','EchC','EchD','EchE','EchF']],
                   'fusion_role':'EchD_on_EchE_N_terminal_domain'}},
        module_prediction_labels={'Ni13':'Mvh–Hdr样H₂氧化与Hdr电子传递耦联（电子分叉为预测）',
            'Ni23':'膜结合、Fd相容的能量转换型[NiFe]氢酶（Ech样功能预测）'},
        module_rule_sources=['https://doi.org/10.1111/j.1742-4658.2005.04889.x',
                            'https://doi.org/10.3389/fmicb.2020.635786',
                            'https://doi.org/10.1073/pnas.1016761108',
                            'https://doi.org/10.1126/science.aan0425',
                            'https://www.rcsb.org/structure/5ODC: separate_HdrA_and_MvhD_entities',
                            'research_workflow/mvh_hdr_model.json:region_route;HdrA2_reference'],
        search=dict(min_identity=35.0,min_reference_coverage=.8,max_evalue=1e-6,max_local_span_bp=30000,
                    catalytic_identity_margin=5.0,architecture_margin=.08),
        source_sha256={name:sha(SOURCE/name) for name in ['complete_evidence.tsv','sequences.tsv','evidence_sequence_links.tsv','protein_domain_architecture.tsv','catalytic_site_details.tsv']},
        fasta_sha256=sha(fasta),scope='Offline functional rule layer; no LLM, no changes to frozen catalytic discovery or identity calibration.')
    policy=ROOT/'config/function_rules_34.json';policy.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    fields=['功能ID','功能类别','共同核心','架构ID','必需成员','可选成员','成熟辅助成员','替代外部伙伴','排除条件','未实测限定','支持模式']
    with (DEST/'rules.tsv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for c in classes:
            for a in c['architectures']:
                roles=lambda k:';'.join(m['role'] for m in a['members'] if m['requirement']==k)
                w.writerow(dict(zip(fields,[c['id'],c['label'],';'.join(c['common_core']),a['id'],roles('required'),roles('optional'),roles('maturation'),roles('alternative_external_partner'),';'.join(c['exclusions']),';'.join(c['unresolved_qualifiers']),a['support_mode']])))
    result=dict(status='complete',classes=len(classes),architectures=sum(len(c['architectures']) for c in classes),reference_sequences=len(used),counter_sequences=len(role_counters),donor_comparison_sequences=len(donor_ids),policy_sha256=sha(policy),fasta_sha256=sha(fasta))
    (DEST/'manifest.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    return payload


if __name__=='__main__':compile_rules()
