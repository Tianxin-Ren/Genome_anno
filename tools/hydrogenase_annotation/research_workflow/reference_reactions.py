"""Keep reference-assay capabilities separate from query activity and net flux."""
import csv
import hashlib
import json
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Data import CodonTable
from Bio.Seq import Seq
from annotate import _circular_contig_length, association_allowed, local_cluster
from complex_extensions import load_bank, function_reviews
from coding_completeness import RULE_VERSION as CODING_RULE, _read_genbank, _record_digest, supports as coding_proof_supports
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from circular_coordinates import circular_slice

ROOT=Path(__file__).resolve().parents[1]
VERSION='reference-assay-potential-v9'
# This is deliberately not a second registry of admitted systems.  The
# reviewed model is the sole source of that scope; this is only the legacy
# extension which remains searchable but is not admitted for exact promotion.
LEGACY_EXTENDED_ONLY={'FeFe_Sfructo_HndABCD'}
ACCEPTED={'calibrated_high','sequence_supported_uncalibrated_core','sequence_supported_uncalibrated_noncanonical'}
TUPLE_POTENTIAL_STATUSES={'exact_reference_sequence_potential','exact_reference_coding_locus_potential'}
DNA_INPUT_MODES={'genome','assembly','bin'}
CAP_FIELDS=('contract_id structure_id reaction_key electron_donors electron_acceptors hydrogen_role '
            'ferredoxin_requirement electron_coupling energy_coupling_scope reference_evidence_scope '
            'primary_source molecular_reaction participant_scope direct_acceptor_status limitations').split()

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def table(path):
    if not path.is_file():return []
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))
def _validate_model_scope(m):
    """Validate the model's trust-boundary shape before using its scope."""
    if not isinstance(m,dict) or m.get('version')!=VERSION:
        raise ValueError('Reference reaction review changed')
    admitted=m.get('admitted_systems')
    if not isinstance(admitted,dict) or not admitted:
        raise ValueError('Reference reaction admitted scope is empty or invalid')
    for sid,cids in admitted.items():
        if not isinstance(sid,str) or not sid.strip() or not isinstance(cids,list) or not cids or \
           len(cids)!=len(set(cids)) or any(not isinstance(cid,str) or not cid.strip() for cid in cids):
            raise ValueError('Invalid admitted reference system/contract scope')
    sources=m.get('sources')
    if not isinstance(sources,dict) or not sources:
        raise ValueError('Reference reaction source binding is missing')
    for name,digest in sources.items():
        path=Path(name) if isinstance(name,str) else None
        if (not isinstance(name,str) or not name or path.is_absolute() or '..' in path.parts or
            not isinstance(digest,str) or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest.lower())):
            raise ValueError('Invalid reference reaction source binding')
    return set(admitted)

def _validate_admitted_bank(model,systems,contracts):
    """Bind every configured admitted structure and contract to the reviewed bank."""
    bysystem={s.get('structure_id'):s for s in systems}
    bycontract={c.get('contract_id'):c for c in contracts}
    for sid,cids in model['admitted_systems'].items():
        system=bysystem.get(sid)
        if system is None:
            raise ValueError('Unknown admitted reference structure: '+sid)
        bank_cids=system.get('reference_contract_ids')
        if not isinstance(bank_cids,list) or set(bank_cids)!=set(cids):
            raise ValueError('Admitted reference contract scope changed: '+sid)
        for cid in cids:
            contract=bycontract.get(cid)
            if contract is None or contract.get('structure_id')!=sid:
                raise ValueError('Unknown admitted reference contract: '+cid)

def load_model(root):
    m=json.loads((root/'research_workflow/reference_reactions_model.json').read_text())
    _validate_model_scope(m)
    return m

def assay_bank(model,root):
    """Load gene-associated assay sets without changing any discovery/search bank."""
    distilled=Path(root)/'references/processed/assay_bank_distilled.json'
    if distilled.is_file():
        d=json.loads(distilled.read_text())
        return d['refs'],d['systems'],d['contracts']
    return [],[],[]

def sha_text(text):return hashlib.sha256(text.encode()).hexdigest()

def _validate_member_groups(groups, selected):
    """Validate an optional role partition without changing tuple evidence."""
    if groups is None:return None
    if not isinstance(groups,list) or not groups or any(not isinstance(group,list) or not group for group in groups):
        raise ValueError('Invalid local member group partition')
    roles={r['role'] for r in selected}
    flat=[role for group in groups for role in group]
    if (any(not isinstance(role,str) or role not in roles for role in flat)
            or len(flat)!=len(set(flat)) or set(flat)!=roles):
        raise ValueError('Invalid local member group partition')
    catalytic=[r['role'] for r in selected if r['kind']=='catalytic']
    if len(catalytic)!=1 or sum(catalytic[0] in group for group in groups)!=1:
        raise ValueError('Invalid catalytic local member group')
    return groups

def options(result):
    yield result
    for child in result.get('alternative_architectures',[]):yield from options(child)


def _local_core_pair(anchor, partner, thresholds):
    """Apply the existing linear gate, or a documented circular local gate."""
    if anchor.get('contig') != partner.get('contig') or anchor.get('strand') != partner.get('strand'):
        return False
    length = _circular_contig_length([anchor, partner])
    if length is not None:
        # The complete tuple, which can bridge this pair, checks actual gaps.
        return length > 0 and local_cluster([anchor, partner], dict(thresholds, local_max_intergenic_gap_bp=length))
    return abs(int(partner.get('start', 0)) - int(anchor.get('start', 0))) <= thresholds['local_max_span_bp']


def _coordinate_distance(anchor, partner):
    """Distance used only for deterministic context ranking."""
    distance = abs(int(partner.get('start', 0)) - int(anchor.get('start', 0)))
    length = _circular_contig_length([anchor, partner])
    if length is None:
        return distance
    if length == 0:
        return float('inf')
    return min(distance % length, length - (distance % length))


def _locus(gene):
    return tuple(gene[k] for k in ('contig', 'start', 'end', 'strand'))


def _coordinate_gate(pid, gene, sequences, coding_proofs, source_inputs, locus=None):
    """Accept linear coordinates, or circular coordinates with source-bound ownership."""
    try:
        start, end, strand, length = (int(gene[k]) for k in ('start', 'end', 'strand', 'contig_length'))
    except (KeyError, TypeError, ValueError):
        return False
    if strand not in {-1, 1} or start < 1 or end < start or length < 0:
        return False
    # Legacy annotation rows may omit contig length; retain their old linear
    # acceptance path, while origin-spanning rows necessarily declare it.
    if length == 0:
        return True
    if start > length:
        return False
    if end <= length:
        return True
    if (str(gene.get('topology', '')).lower() != 'circular' or
            end > start + length - 1 or pid not in sequences):
        return False
    # A circular interval cannot enter tuple enumeration on coordinates alone:
    # either the coding proof or the independently reconstructed locus must
    # bind the origin-spanning DNA to the declared source FASTA.
    return bool(locus is not None or _coding_proof_valid(pid, gene, sequences[pid], coding_proofs, source_inputs))


def _context_rank(anchor, partner):
    return (partner.get('contig') != anchor.get('contig'), partner.get('strand') != anchor.get('strand'),
            _coordinate_distance(anchor, partner))


def capability(contract):return {**{k:contract[k] for k in CAP_FIELDS},'scope':'reference_assay_conditions_only'}

def _coding_proof_valid(pid,g,seq,coding_proofs,source_inputs):
    """Return true only for a current, source-bound complete-CDS proof."""
    proof=coding_proofs.get(pid) if coding_proofs else None
    if not proof or not source_inputs:return False
    artifact=proof.get('source_artifact',{})
    path=str(artifact.get('path',''))
    if path not in source_inputs or source_inputs.get(path)!=artifact.get('sha256'):
        return False
    return bool(coding_proof_supports(proof,pid,
        expected_protein_sha256=hashlib.sha256(seq.encode()).hexdigest(),gene=g))

def _coding_locus_match(pid,ref,gene,seq,proofs,source_inputs,run):
    """Bind a registered member to DNA at the called CDS's 5' boundary."""
    if run.get('input_mode') not in DNA_INPUT_MODES or gene.get('coordinate_evidence')!='genomic_CDS' or not proofs or not source_inputs:
        return None
    reference=ref.get('sequence',''); qlen=len(seq); rlen=len(reference)
    if not _boundary_shape(reference,seq) or not isinstance(proofs.get(pid),dict): return None
    proof=proofs[pid]; artifact=proof.get('source_artifact',{}); path=str(artifact.get('path',''))
    if path not in source_inputs or source_inputs[path]!=artifact.get('sha256'): return None
    try:
        table=int(proof['translation_table']); start=int(gene['start']); end=int(gene['end']); strand=int(gene['strand'])
        contig=str(gene['contig'])
        if proof.get('translation_table_source') == 'genbank_CDS_qualifier':
            dna, _, _ = _read_genbank(path)
        else:
            with Path(path).open() as handle: records=list(SeqIO.parse(handle,'fasta'))
            if len({r.id for r in records}) != len(records): return None
            dna={r.id:str(r.seq).upper() for r in records}
        contig_len=int(gene.get('contig_length',0))
    except (KeyError,TypeError,ValueError,OSError): return None
    if contig not in dna or strand not in {-1,1} or contig_len<=0 or contig_len!=len(dna[contig]): return None
    crossing = end > contig_len
    if crossing and (str(gene.get('topology', '')).lower() != 'circular' or
                     int(gene.get('contig_length', 0) or 0) != contig_len or
                     start > contig_len or end > start + contig_len - 1): return None
    if not crossing and end > contig_len: return None
    complete=_coding_proof_valid(pid,gene,seq,proofs,source_inputs)
    edge=False
    if not complete:
        distance=start-1 if strand==1 else contig_len-end
        edge=(qlen>rlen and 0<=distance<=2 and _circular_contig_length([gene])==contig_len
              and str(gene.get('partial')).lower()=='true'
              and proof.get('rule_version')==CODING_RULE and proof.get('record_digest')==_record_digest(proof)
              and proof.get('protein_id')==pid and proof.get('protein_sha256')==sha_text(seq)
              and (proof.get('contig'),proof.get('start'),proof.get('end'),proof.get('strand'))==(contig,start,end,strand)
              and proof.get('status')=='incomplete' and proof.get('reasons')==['canonical_start_missing']
              and proof.get('source_partial') is True and sha(Path(path))==artifact.get('sha256'))
        if not edge: return None
    diff=rlen-qlen
    missing,extra=max(diff,0),max(-diff,0)
    ext_start,ext_end=(start-3*diff,end) if strand==1 else (start,end+3*diff)
    ext_crossing = ext_end > contig_len
    if ext_crossing and (str(gene.get('topology', '')).lower() != 'circular' or
                         int(gene.get('contig_length', 0) or 0) != contig_len or
                         ext_start > contig_len or ext_end > ext_start + contig_len - 1): return None
    if (ext_start<1 or (not ext_crossing and ext_end>contig_len) or ext_end<ext_start or
            ext_end-ext_start+1!=3*(rlen+1)): return None
    nt=circular_slice(dna[contig], ext_start, ext_end) if ext_crossing else dna[contig][ext_start-1:ext_end]
    if strand==-1: nt=str(Seq(nt).reverse_complement())
    table_obj=CodonTable.unambiguous_dna_by_id.get(table)
    if table_obj is None or len(nt)%3 or nt[:3] not in table_obj.start_codons: return None
    try: translated=str(Seq(nt).translate(table=table,cds=True))
    except (ValueError,KeyError,CodonTable.TranslationError): return None
    if translated!=reference: return None
    query_nt=circular_slice(dna[contig], start, end) if crossing else dna[contig][start-1:end]
    if strand==-1: query_nt=str(Seq(query_nt).reverse_complement())
    if hashlib.sha256(query_nt.encode()).hexdigest()!=proof.get('dna_sha256'): return None
    try:
        if edge:
            if len(query_nt)!=3*(qlen+1) or query_nt[:3] in table_obj.start_codons: return None
            raw_translation=str(Seq(query_nt).translate(table=table))
            if raw_translation != seq+'*': return None
            raw_query=raw_translation[:-1]
        else:
            raw_query=str(Seq(query_nt).translate(table=table,cds=True))
        # BioPython's cds=True normalizes alternate starts to M.  The coding
        # proof deliberately preserves the raw initiator residue, so restore
        # the table's ordinary codon translation for that one position.
        if proof.get('translation_match_mode')=='raw_legal_alternate_initiator':
            raw_query=table_obj.forward_table.get(query_nt[:3], raw_query[0])+raw_query[1:]
    except (ValueError,KeyError,CodonTable.TranslationError): return None
    if raw_query!=seq: return None
    mode='raw_reference_suffix' if diff>0 else ('normalized_query_initiator' if diff==0 else 'query_has_extra_n_terminal_residues')
    if edge: mode='circular_edge_exact_reference_after_trimming_partial_start'
    return {'query_protein_sha256':hashlib.sha256(seq.encode()).hexdigest(),
            'reference_protein_sha256':hashlib.sha256(reference.encode()).hexdigest(),
            'query_nt_sha256':hashlib.sha256(query_nt.encode()).hexdigest(),
            'extended_nt_sha256':hashlib.sha256(nt.encode()).hexdigest(),
            'query_translation_boundary_different':True,'query_start':start,'query_end':end,
            'boundary_translation_mode':mode,'encoded_reference_start':ext_start,'encoded_reference_end':ext_end,
            'missing_n_terminal_residues':missing,'extra_n_terminal_residues':extra,
            'source_input_path':path,'source_input_sha256':artifact.get('sha256'),'translation_table':table,
            'reference_sequence_sha256':ref.get('sequence_sha256')}

def coding_locus_candidate_ids(annotation,root=ROOT):
    """Return genome candidates whose called protein has a registered 5' boundary."""
    annotation,root=Path(annotation),Path(root)
    run=json.loads((annotation/'run_manifest.json').read_text())
    if run.get('input_mode') not in DNA_INPUT_MODES: return []
    model=load_model(root); refs,_,_=assay_bank(model,root)
    faa=annotation/'proteins.faa'; genes=table(annotation/'genes.tsv')
    with faa.open() as handle: proteins={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    return sorted({g['id'] for g in genes if g.get('coordinate_evidence')=='genomic_CDS' and g['id'] in proteins
                   for ref in refs if _boundary_shape(ref['sequence'],proteins[g['id']])})

def _boundary_shape(reference,seq):
    if not reference or not seq or reference==seq: return False
    if len(reference)==len(seq): return len(reference)>1 and reference[1:]==seq[1:] and reference[0]!=seq[0]
    d=len(reference)-len(seq)
    if d>0:
        return reference[d:]==seq or (len(seq)>1 and reference[d+1:]==seq[1:])
    return seq[-len(reference):]==reference or (len(reference)>1 and seq[-len(reference)+1:]==reference[1:])

def tuple_proof(qid,result,refs,genes,sequences,mode,thresholds,coding_proofs=None,source_inputs=None,coding_locus=None,member_groups=None):
    """Protein equivalence is a separate claim from a complete genomic CDS replay."""
    member_groups=_validate_member_groups(member_groups,refs)
    errors=[];matched=[];used=set();roles=result.get('evidence',{}).get('roles',{})
    direct=result.get('entry_route') in {'exact_full_protein_lookup','exact_coding_locus_lookup'}
    if (not direct and result.get('coupling_potential_status')!='local_coupling_potential_uncalibrated') or result.get('evidence',{}).get('resource_conflict'):
        errors.append('local_architecture_or_resource_ownership_unresolved')
    if mode not in DNA_INPUT_MODES:errors.append('local_sequence_input_scope_required')
    anchor=genes.get(qid)
    role_group={role:index for index,group in enumerate(member_groups or []) for role in group}
    for r in refs:
        candidates=roles.get(r['role'],[])
        if len(candidates)!=1:errors.append('missing_or_ambiguous_exact_reference_role_'+r['role']);continue
        pid=candidates[0]['protein_id'];g=genes.get(pid);seq=sequences.get(pid)
        if pid in used:errors.append('member_resource_reuse_'+r['role'])
        used.add(pid)
        if g is None or seq is None:raise ValueError('Reference tuple refers to an unknown query protein')
        if r['kind']=='catalytic' and pid!=qid:errors.append('catalytic_role_anchor_mismatch')
        locus=(coding_locus or {}).get(pid)
        if locus is not None and locus.get('reference_sequence_sha256')!=r.get('sequence_sha256'):
            if seq!=r['sequence']:errors.append('coding_locus_reference_role_mismatch_'+r['role'])
            locus=None
        raw_partial=str(g.get('partial')).lower() not in {'false','0'}
        proof_used=False
        if raw_partial:
            proof_used=_coding_proof_valid(pid,g,seq,coding_proofs,source_inputs)
            if not proof_used:
                edge_locus=(locus or {}).get('boundary_translation_mode')=='circular_edge_exact_reference_after_trimming_partial_start'
                if not edge_locus or locus!=_coding_locus_match(pid,r,g,seq,coding_proofs,source_inputs,{'input_mode':mode}):
                    errors.append('partial_or_unknown_member_'+r['role'])
        if (seq!=r['sequence'] and locus is None) or not set(seq)<=set('ACDEFGHIKLMNPQRSTVWY'):
            errors.append('not_exact_unambiguous_full_sequence_'+r['role'])
        start,end,strand,length=(int(g[k]) for k in ['start','end','strand','contig_length'])
        if not _coordinate_gate(pid, g, sequences, coding_proofs, source_inputs, locus=locus):
            errors.append('invalid_supplied_coordinate_'+r['role'])
        if end-start+1 not in {3*len(seq),3*(len(seq)+1)}:errors.append('protein_CDS_coordinate_length_conflict_'+r['role'])
        if anchor is None or not association_allowed(anchor,g,mode):errors.append('member_input_association_conflict_'+r['role'])
        matched.append(dict(role=r['role'],protein_id=pid,reference_id=r['reference_id'],native_reference_protein=r['protein_accession'],
            query_sha256=hashlib.sha256(seq.encode()).hexdigest(),reference_sha256=r['sequence_sha256'],
            contig=g['contig'],start=start,end=end,strand=strand,contig_length=length or 'not_supplied',
            raw_partial=raw_partial,coding_proof_used=proof_used,
            coding_locus_proof_used=bool(locus),
            coding_proof_digest=(coding_proofs[pid].get('record_digest') if proof_used or locus else None),
            **({'member_group':list(member_groups[role_group[r['role']]]) ,'group_index':role_group[r['role']]} if member_groups and r['role'] in role_group else {}),
            **(locus or {})))
    lengths={}
    for m in matched:
        if m['contig_length']!='not_supplied':
            lengths.setdefault(m['contig'],set()).add(m['contig_length'])
    if any(len(v)>1 for v in lengths.values()):errors.append('inconsistent_supplied_contig_lengths')
    if member_groups:
        for index,group in enumerate(member_groups):
            members=[genes[m['protein_id']] for m in matched if m.get('group_index')==index]
            if members and not local_cluster(members,thresholds):errors.append('reference_member_group_local_linkage_missing_'+str(index))
    elif matched and not local_cluster([genes[r['protein_id']] for r in matched],thresholds):errors.append('whole_reference_tuple_local_linkage_missing')
    return sorted(set(errors)),matched

def _component_hypothesis(qid, result, refs, roles, genes, sequences, mode,
                          thresholds, coding_proofs, source_inputs, coding_locus,
                          member_groups):
    """Choose one complete source-equivalent tuple while leaving locality unresolved."""
    linkage=lambda error: error=='whole_reference_tuple_local_linkage_missing' or error.startswith('reference_member_group_local_linkage_missing_')
    ordered_refs=sorted(refs,key=lambda ref:(len(roles.get(ref['role'],())), ref['kind']!='catalytic', ref['role']))
    anchor=genes.get(qid)
    def rank(pid):
        return _context_rank(anchor, genes[pid])
    pools={ref['role']:sorted({x['protein_id'] for x in roles.get(ref['role'],())},key=lambda pid:(rank(pid),pid)) for ref in refs}
    selected={};found=None
    def visit(index,used):
        nonlocal found
        if found is not None:return
        if index==len(ordered_refs):
            trial={name:[{'protein_id':pid}] for name,pid in selected.items()}
            trial_result={'entry_route':result.get('entry_route'),'evidence':{'roles':trial,'resource_conflict':result.get('evidence',{}).get('resource_conflict',False)}}
            errors,matched=tuple_proof(qid,trial_result,refs,genes,sequences,mode,thresholds,coding_proofs,source_inputs,coding_locus,member_groups)
            if all(linkage(error) for error in errors):found=(dict(selected),errors,matched)
            return
        ref=ordered_refs[index]
        for pid in pools.get(ref['role'],()):
            if pid in used:continue
            selected[ref['role']]=pid;visit(index+1,used|{pid})
            if found is not None:return
        selected.pop(ref['role'],None)
    visit(0,set())
    if found is None:return None
    chosen,linkage_errors,matched=found
    ambiguities=[]
    for ref in refs:
        role=ref['role'];alternatives=[]
        for pid in pools.get(role,()):
            if pid==chosen[role]:continue
            trial=dict(chosen);trial[role]=pid
            if len(set(trial.values()))!=len(trial):continue
            trial_result={'entry_route':result.get('entry_route'),'evidence':{'roles':{name:[{'protein_id':value}] for name,value in trial.items()},'resource_conflict':result.get('evidence',{}).get('resource_conflict',False)}}
            errors,_=tuple_proof(qid,trial_result,refs,genes,sequences,mode,thresholds,coding_proofs,source_inputs,coding_locus,member_groups)
            if all(linkage(error) for error in errors):alternatives.append(pid)
        if alternatives:
            tied = [pid for pid in alternatives if rank(pid)==rank(chosen[role])
                    and _locus(genes[pid])!=_locus(genes[chosen[role]])]
            ambiguities.append({'role':role,'selected_protein_id':chosen[role],'alternative_protein_ids':alternatives,
                                'selection_basis':'source_equivalent_complete_members_ranked_by_context_then_id',
                                'query_member_assay_status':'not_tested','literal_member_assignment_status':'not_unique_sequence_equivalent_loci',
                                **({'equally_supported_loci':tied} if tied else {})})
    selected_result={'entry_route':result.get('entry_route'),'evidence':{'roles':{name:[{'protein_id':value}] for name,value in chosen.items()},
                                                                           'resource_conflict':result.get('evidence',{}).get('resource_conflict',False)}}
    return selected_result,linkage_errors,matched,ambiguities

def _source_like_fdh(qid, ref, inventory, genes, sequences, mode, anchor, used, coding_proofs=None, source_inputs=None):
    """Select a source-class Fdh-like module without changing exact tuple proof."""
    source_sha=ref.get('sequence_sha256')
    classes=(inventory or {}).get('fdh_source_member_classes',{})
    if not source_sha or classes.get(source_sha)!='Fdh_molybdo': return None,[],[]
    candidates=[]
    for pid,rows in (inventory or {}).get('gene_roles',{}).items():
        if pid in used or pid not in genes or pid not in sequences: continue
        for row in rows:
            chem=row.get('fdh_chemistry',{})
            compat_rows=row.get('source_member_compatibility',[])
            if not isinstance(compat_rows,list): continue
            compat=next((x for x in compat_rows if isinstance(x,dict) and x.get('source_sequence_sha256')==source_sha),None)
            if (row.get('role_class')!='Fdh_molybdo' or row.get('preferred_role_class') not in {None,'Fdh_molybdo'} or
                row.get('status')!='supported_role_homology' or
                chem.get('status')!='formate_class_compatible' or
                compat is None or
                compat.get('status')!='source_like_module_supported' or
                not isinstance(compat.get('source_profile_proof'),list) or len(compat['source_profile_proof'])!=4): continue
            proof=compat['source_profile_proof']
            try:
                proof_ranges=[]
                for p in proof:
                    if not isinstance(p,dict) or not p.get('accession') or str(p.get('passes','')).lower()!='true':
                        raise ValueError
                    if p.get('query_id') != hashlib.sha256(sequences[pid].encode()).hexdigest() or int(p.get('query_length')) != len(sequences[pid]):
                        raise ValueError
                    sb,db,sga,dga=[float(p[k]) for k in ('sequence_bitscore','domain_bitscore','sequence_ga','domain_ga')]
                    if not all(x == x and x != float('inf') and x != float('-inf') for x in (sb,db,sga,dga)) or sb<sga or db<dga:
                        raise ValueError
                    proof_ranges.append((int(p['query_start']),int(p['query_end'])))
            except (KeyError,TypeError,ValueError):
                continue
            if (any(start < 1 or end < start or end > len(sequences[pid]) for start,end in proof_ranges) or
                proof_ranges != sorted(proof_ranges) or any(b[0] <= a[1] for a,b in zip(proof_ranges,proof_ranges[1:]))):
                continue
            seq=sequences[pid]; gene=genes[pid]
            if hashlib.sha256(seq.encode()).hexdigest()==source_sha: continue
            partial=str(gene.get('partial')).lower() not in {'false','0'}
            proof_used=_coding_proof_valid(pid,gene,seq,coding_proofs,source_inputs) if partial else False
            if partial and not proof_used: continue
            if not association_allowed(anchor,gene,mode): continue
            if gene.get('strand') not in {-1,1}: continue
            if not _coordinate_gate(pid, gene, sequences, coding_proofs, source_inputs): continue
            if gene['end']-gene['start']+1 not in {3*len(seq),3*(len(seq)+1)}: continue
            try:
                identity=float(compat.get('paired_identity',0)); rcover=float(compat.get('reference_coverage',0)); qcover=float(compat.get('query_coverage',0))
            except (TypeError,ValueError): continue
            if not 0<identity<=100 or not 0<rcover<=1 or not 0<qcover<=1: continue
            candidates.append((row, pid, compat, chem, proof_used))
    if not candidates:return None,[],[]
    def rank(item):
        row,pid,compat,_,_=item; g=genes[pid]
        return (-float(compat.get('paired_identity',0)), -float(compat.get('reference_coverage',0)),
                -float(compat.get('query_coverage',0)), g.get('contig')!=anchor.get('contig'),
                _coordinate_distance(anchor, g), pid)
    candidates.sort(key=rank); row,pid,compat,chem,proof_used=candidates[0]
    selected_partial=str(genes[pid].get('partial')).lower() not in {'false','0'}
    record=dict(role=ref['role'],protein_id=pid,reference_id=ref['reference_id'],native_reference_protein=ref['protein_accession'],
        matched_by='supported_FDH_role_homology',source_equivalent=False,query_sha256=hashlib.sha256(sequences[pid].encode()).hexdigest(),
        reference_sha256=source_sha,contig=genes[pid]['contig'],start=genes[pid]['start'],end=genes[pid]['end'],strand=genes[pid]['strand'],
        raw_partial=selected_partial,coding_proof_used=bool(proof_used),coding_locus_proof_used=False,paired_identity=compat.get('paired_identity'),
        reference_coverage=compat.get('reference_coverage'),query_coverage=compat.get('query_coverage'),source_profile_proof=compat.get('source_profile_proof'),
        fdh_chemistry=chem)
    alternatives=[{'role':ref['role'],'selected_protein_id':pid,'alternative_protein_ids':[x[1] for x in candidates[1:]],
        'selection_basis':'source_specific_identity_coverage_context_then_id','query_member_assay_status':'not_tested',
        'literal_member_assignment_status':'not_source_exact'}] if len(candidates)>1 else []
    return record,alternatives,[]


def _source_module_prior(qid, result, refs, component_roles, genes, sequences, mode,
                         thresholds, pyridine_results, caps, coding_proofs=None,
                         source_inputs=None, coding_locus=None, member_groups=None):
    """Use one independently supported module homolog as a coarse source-function prior."""
    missing=[ref for ref in refs if not component_roles.get(ref['role'])]
    if len(missing)!=1 or missing[0]['kind']!='module':return []
    ref=missing[0];backbone=[r for r in refs if r is not ref]
    if not backbone or any(not component_roles.get(r['role']) for r in backbone):return []
    groups=None if member_groups is None else [[role for role in group if role!=ref['role']] for group in member_groups]
    if groups is not None:groups=[group for group in groups if group]
    found=_component_hypothesis(qid,result,backbone,component_roles,genes,sequences,mode,
        thresholds,coding_proofs,source_inputs,coding_locus,groups)
    if found is None:return []
    _selected,_linkage,matched,_ambiguities=found
    used={m['protein_id'] for m in matched}
    participants=sorted({p for cap in caps for name in ['electron_donors','electron_acceptors'] for p in cap[name]})
    hypotheses=[]
    for current in options(pyridine_results.get(qid,{})):
        evidence=current.get('evidence',{})
        if (current.get('coupling_potential_status')!='local_coupling_potential_uncalibrated'
                or current.get('hydrogen_reaction_core_status')!='sequence_supported_candidate'
                or current.get('coupling_blockers') or not evidence.get('association_supported')
                or evidence.get('resource_conflict') or evidence.get('actual_resource_overlap')):continue
        roles=evidence.get('roles',{})
        members=[m for rows in roles.values() for m in rows if m.get('protein_id')]
        pool={m['protein_id'] for m in members}
        if not used<=pool or len(pool-used)!=1 or len(pool)!=len(refs):continue
        pid=next(iter(pool-used));query=sequences.get(pid,'')
        if not query or not set(query)<=set('ACDEFGHIKLMNPQRSTVWY') or query==ref['sequence']:continue
        module_members=[(role,m) for role,rows in roles.items() for m in rows if m.get('protein_id')==pid]
        if len(module_members)!=1:continue
        role,module=module_members[0]
        if any(key not in module for key in ['identity','query_coverage','reference_coverage']):continue
        aligner=PairwiseAligner(mode='global',match_score=2,mismatch_score=-1,open_gap_score=-10,extend_gap_score=-.5)
        alignment=aligner.align(ref['sequence'],query)[0]
        a,b=str(alignment[0]),str(alignment[1])
        paired=sum(x!='-' and y!='-' for x,y in zip(a,b))
        identical=sum(x==y and x!='-' for x,y in zip(a,b))
        identity=100*identical/len(a)
        qcov,rcov=paired/len(query),paired/len(ref['sequence'])
        if (identity<thresholds['min_identity'] or qcov<thresholds['min_query_coverage']
                or rcov<thresholds['min_reference_coverage']):continue
        selected={m['role']:[{'protein_id':m['protein_id']}] for m in matched}
        selected[ref['role']]=[{'protein_id':pid}]
        trial={'entry_route':result.get('entry_route'),'evidence':{'roles':selected,'resource_conflict':False}}
        errors,all_matched=tuple_proof(qid,trial,refs,genes,sequences,mode,thresholds,
            coding_proofs,source_inputs,coding_locus,member_groups)
        if errors!=['not_exact_unambiguous_full_sequence_'+ref['role']]:continue
        source_alignment=dict(identity=identity,query_coverage=qcov,reference_coverage=rcov,
            query_similarity=100*identical/len(query),
            query_sha256=sha_text(query),reference_sha256=ref['sequence_sha256'],
            query_length=len(query),reference_length=len(ref['sequence']),
            aligned_reference=a,aligned_query=b,score=float(alignment.score),
            method='global_pairwise_sequence_comparison_no_evalue',
            parameters=dict(match=2,mismatch=-1,open_gap=-10,extend_gap=-.5))
        # Compare on the same query denominator; different reference lengths are not a ranking gate.
        template_similarity=float(module['identity'])*float(module['query_coverage'])
        dominates=source_alignment['query_similarity']>=template_similarity
        for member in all_matched:
            member['source_equivalent']=member['protein_id']!=pid
            if member['protein_id']==pid:
                member.update(matched_by='supported_local_module_homology',source_alignment=source_alignment,
                    role_module_evidence=dict(structure_id=current['structure_id'],role=role,member=dict(module),
                        query_similarity=template_similarity),
                    source_alignment_dominates_role_template=dominates)
        hypotheses.append(dict(status='source_like_module_component_potential',structure_id=ref['structure_id'],
            scope='source_like_module_function_prior',matched_members=all_matched,homologous_roles=[ref['role']],
            source_participants=participants,source_sources=sorted({cap['primary_source'] for cap in caps}),
            prefer_source_function=dominates,confidence='medium' if dominates and current.get('confidence_hint')!='low' else 'low',
            selection_risks=list(current.get('selection_risks',[])),hard_blockers=[],
            strict_tuple_blockers=errors,capabilities=[],source_contracts_transferred=False,
            reference_reaction_details_transferred=False,query_activity='not_tested',in_vivo_direction='unknown',
            physical_assembly='unknown',native_partner_specificity='unknown',
            limitations='One nonexact noncatalytic member in a complete otherwise source-equivalent local tuple. Coarse function prior only; unvalidated homology transfer, no exact assay attribution, numerical rate, nucleotide specificity, mechanism or net direction.'))
    return hypotheses

def prepare(annotation,anchors,identities,pyridine_results,root=ROOT,coding_completeness=None,module_role_inventory=None):
    """Project already assigned contracts, or prove an admitted exact protein tuple."""
    annotation,root=Path(annotation),Path(root)
    model=load_model(root)
    refs,systems,_,_,contracts=load_bank(root)
    extra_refs,extra_systems,extra_contracts=assay_bank(model,root)
    if {s['structure_id'] for s in systems}&{s['structure_id'] for s in extra_systems} or {c['contract_id'] for c in contracts}&{c['contract_id'] for c in extra_contracts}:raise ValueError('Assay bank overwrites an existing reference')
    refs=refs+extra_refs;systems=systems+extra_systems;contracts=contracts+extra_contracts
    _validate_admitted_bank(model,systems,contracts)
    admitted=set(model['admitted_systems'])
    extended=admitted|LEGACY_EXTENDED_ONLY
    candidates={q:[o for o in options(r) if o.get('structure_id') in extended] for q,r in pyridine_results.items()}
    candidates={q:v for q,v in candidates.items() if v}
    inherited=table(annotation/'reaction_capabilities.tsv')
    run=json.loads((annotation/'run_manifest.json').read_text());bound={str(annotation/'run_manifest.json'):sha(annotation/'run_manifest.json')}
    source_inputs={str(item.get('path')):item.get('sha256') for item in run.get('input_files',[])
                   if item.get('path') and item.get('sha256')}
    coding_proofs=coding_completeness if isinstance(coding_completeness,dict) else {}
    if not anchors or (not candidates and not inherited and not any(identities.get(a['protein_id']) in ACCEPTED for a in anchors)):return {},{}
    bycontract={c['contract_id']:c for c in contracts}
    faa=annotation/'proteins.faa'
    with faa.open() as handle:records=list(SeqIO.parse(handle,'fasta'))
    sequences={r.id:str(r.seq) for r in records}
    if len(sequences)!=len(records):raise ValueError('Duplicate reference reaction protein ID')
    exact_anchors={a['protein_id']:[r['structure_id'] for r in refs if r['kind']=='catalytic' and r['structure_id'] in admitted and sequences.get(a['protein_id'])==r['sequence']] for a in anchors}
    exact_anchors={q:sids for q,sids in exact_anchors.items() if sids}
    genes_rows=table(annotation/'genes.tsv');genes={g['id']:g for g in genes_rows}
    for g in genes.values():
        for k in ['start','end','strand','contig_length']:g[k]=int(g[k])
    locus_matches={}
    assay_structures={s['structure_id'] for s in extra_systems}
    for a in anchors:
        qid=a['protein_id']; gene=genes.get(qid)
        if gene is None or identities.get(qid) not in ACCEPTED or a.get('input_mode') not in DNA_INPUT_MODES: continue
        for ref in refs:
            if ref['kind']=='catalytic' and ref['structure_id'] in admitted and ref['structure_id'] in assay_structures:
                hit=_coding_locus_match(qid,ref,gene,sequences.get(qid,''),coding_proofs,source_inputs,run)
                if hit: locus_matches.setdefault(qid,{}).setdefault(ref['structure_id'],{})[qid]=hit
    locus_anchors={q:list(sids) for q,sids in locus_matches.items()}
    if not candidates and not inherited and not exact_anchors and not locus_anchors:return {},{}
    reviews=function_reviews(contracts,root);byanchor={a['protein_id']:a for a in anchors}
    if len(byanchor)!=len(anchors):raise ValueError('Duplicate reference reaction anchor')
    for name in ['genes.tsv','catalytic_anchors.tsv','complex_members.tsv','reaction_capabilities.tsv']:
        path=annotation/name
        if not path.exists():
            if name=='reaction_capabilities.tsv' and not inherited:continue
            raise ValueError('Missing reference reaction input: '+name)
        if run.get('output_sha256',{}).get(name)!=sha(path):raise ValueError('Reference reaction annotation output changed: '+name)
        bound[str(path)]=sha(path)
    if run.get('normalized_proteins_sha256')!=sha(faa):raise ValueError('Reference reaction proteins changed')
    bound[str(faa)]=sha(faa)
    bound.update({str(root/n):sha(root/n) for n in ['research_workflow/reference_reactions.py','research_workflow/reference_reactions_model.json']})
    if len(sequences)!=len(records) or len(genes)!=len(genes_rows) or set(sequences)!=set(genes):raise ValueError('Reference reaction protein/gene ownership mismatch')
    if coding_proofs:bound[str(root/'research_workflow/coding_completeness.py')]=sha(root/'research_workflow/coding_completeness.py')
    for pid,g in genes.items():
        if _coding_proof_valid(pid,g,sequences[pid],coding_proofs,source_inputs):
            artifact=coding_proofs[pid]['source_artifact'];bound[artifact['path']]=artifact['sha256']
    sequence_ids={}
    for pid,seq in sequences.items():sequence_ids.setdefault(seq,[]).append(pid)
    bysystem={s['structure_id']:s for s in systems}
    for qid,sids in {q:sorted(set(exact_anchors.get(q,[]))|set(locus_anchors.get(q,[]))) for q in set(exact_anchors)|set(locus_anchors)}.items():
        anchor=genes[qid];mode=byanchor[qid]['input_mode']
        for sid in sids:
            roles={}; local_members=dict(locus_matches.get(qid,{}).get(sid,{}))
            component_roles={}; component_locus=dict(local_members)
            member_groups=bysystem[sid].get('local_member_groups')
            rs=[r for r in refs if r['structure_id']==sid]
            catalytic_role=next((r['role'] for r in rs if r['kind']=='catalytic'),None)
            core_group_index=next((i for i,g in enumerate(member_groups or []) if catalytic_role in g),None)
            for ref in rs:
                ids=[qid] if ref['kind']=='catalytic' else sequence_ids.get(ref['sequence'],[])
                role_is_core=core_group_index is None or ref['role'] in member_groups[core_group_index]
                ids=[pid for pid in ids if association_allowed(anchor,genes[pid],mode)
                     and (not role_is_core or _local_core_pair(anchor, genes[pid], run['thresholds']))]
                role_locus={}
                if ref['kind']=='module':
                    for pid,g in genes.items():
                        if pid==qid or pid in ids or (role_is_core and not _local_core_pair(anchor, g, run['thresholds'])):
                            continue
                        if not association_allowed(anchor,g,mode): continue
                        hit=_coding_locus_match(pid,ref,g,sequences[pid],coding_proofs,source_inputs,run)
                        if hit: role_locus[pid]=hit
                    local_members.update(role_locus)
                    ids=ids+list(role_locus)
                # Recompute after adding only this role's verified boundary
                # members; a locus hit for another role must not make this
                # role ambiguous.
                valid= [pid for pid in ids if pid in genes and pid in sequences and
                        _coordinate_gate(pid, genes[pid], sequences, coding_proofs, source_inputs,
                                         locus=role_locus.get(pid)) and
                        (pid in role_locus or str(genes[pid]['partial']).lower() in {'false','0'} or
                         _coding_proof_valid(pid,genes[pid],sequences[pid],coding_proofs,source_inputs))]
                roles[ref['role']]=[dict(protein_id=pid) for pid in dict.fromkeys(valid or ids)]
                # Component hypotheses deliberately omit the strict local
                # contig/strand/span prefilter, but retain association,
                # complete sequence and coding-proof gates.
                component_ids=[qid] if ref['kind']=='catalytic' else sequence_ids.get(ref['sequence'],[])
                component_ids=[pid for pid in component_ids if pid in genes and association_allowed(anchor,genes[pid],mode)]
                component_role_locus={}
                if ref['kind']=='module':
                    for pid,g in genes.items():
                        if pid==qid or not association_allowed(anchor,g,mode):continue
                        hit=_coding_locus_match(pid,ref,g,sequences[pid],coding_proofs,source_inputs,run)
                        if hit:component_role_locus[pid]=hit
                    component_locus.update(component_role_locus)
                    component_ids.extend(component_role_locus)
                component_valid=[pid for pid in dict.fromkeys(component_ids) if
                    _coordinate_gate(pid, genes[pid], sequences, coding_proofs, source_inputs,
                                     locus=component_role_locus.get(pid)) and
                    (pid in component_role_locus or str(genes[pid]['partial']).lower() in {'false','0'} or
                     _coding_proof_valid(pid,genes[pid],sequences[pid],coding_proofs,source_inputs))]
                component_roles[ref['role']]=[dict(protein_id=pid) for pid in component_valid]
            locus_matches.setdefault(qid,{})[sid]=local_members
            conflicts=[c for option in candidates.get(qid,[]) for c in option.get('evidence',{}).get('cross_anchor_conflicts',[]) if identities.get(c['rival_anchor']) in ACCEPTED]
            candidates[qid]=[r for r in candidates.get(qid,[]) if r['structure_id']!=sid]
            is_locus=bool(locus_matches.get(qid,{}).get(sid)) and sid not in exact_anchors.get(qid,[])
            selection_ambiguities=[];selection_rejections=[]
            route='exact_coding_locus_lookup' if is_locus else 'exact_full_protein_lookup'
            if route in {'exact_full_protein_lookup','exact_coding_locus_lookup'} and member_groups:
                singleton_roles=[group[0] for group in member_groups if len(group)==1]
                ambiguous=[role for role in singleton_roles if next((r for r in rs if r['role']==role and r['kind']!='catalytic'),None)
                           and len(roles.get(role,[]))>1]
                if len(ambiguous)==1 and all(len(roles.get(ref['role'],[]))==1 for ref in rs if ref['role']!=ambiguous[0]):
                    role=ambiguous[0];qualified=[]
                    for candidate in roles.get(role,[]):
                        test_roles={name:list(values) for name,values in roles.items()}
                        test_roles[role]=[candidate]
                        test_result={'entry_route':route,'evidence':{'roles':test_roles,'resource_conflict':bool(conflicts)}}
                        candidate_locus=dict(local_members) if local_members else None
                        candidate_errors,_=tuple_proof(qid,test_result,rs,genes,sequences,mode,run['thresholds'],coding_proofs,source_inputs,
                                                        coding_locus=candidate_locus,member_groups=member_groups)
                        if candidate_errors:
                            selection_rejections.append({'role':role,'protein_id':candidate['protein_id'],'blockers':candidate_errors})
                        else: qualified.append(candidate['protein_id'])
                    if qualified:
                        anchor_gene=genes[qid]
                        def rank(pid):
                            return _context_rank(anchor_gene, genes[pid])
                        ordered=sorted(qualified,key=lambda pid:(rank(pid),pid))
                        first=ordered[0]
                        tied=[pid for pid in ordered[1:] if rank(pid)==rank(first)
                              and _locus(genes[pid])!=_locus(genes[first])]
                        roles[role]=[{'protein_id':pid} for pid in [first,*tied]]
                        if len(ordered)>1 and not tied:
                            selection_ambiguities=[{'role':role,'selected_protein_id':ordered[0],
                                'alternative_protein_ids':ordered[1:],
                                'selection_basis':'source_equivalent_complete_members_ranked_by_context_then_id',
                                'query_member_assay_status':'not_tested',
                                'literal_member_assignment_status':'not_unique_sequence_equivalent_loci'}]
            candidates[qid].insert(0,dict(structure_id=sid,entry_route='exact_coding_locus_lookup' if is_locus else 'exact_full_protein_lookup',
                coding_locus=locus_matches.get(qid,{}).get(sid),evidence=dict(roles=roles,resource_conflict=bool(conflicts),
                **({'member_selection_ambiguities':selection_ambiguities} if selection_ambiguities else {}),
                **({'member_selection_rejections':selection_rejections} if selection_rejections else {}),
                component_roles=component_roles,component_locus=component_locus)))
    original_rows=table(annotation/'catalytic_anchors.tsv');original={a['protein_id']:a for a in original_rows}
    if len(original)!=len(original_rows):raise ValueError('Duplicate original reference reaction anchor')
    results={}
    def entry(qid):
        if qid not in byanchor:raise ValueError('Reference reaction belongs to an unknown anchor')
        return results.setdefault(qid,dict(status='not_supported',capabilities=[],tuple_matches=[],blockers=[],inherited_assignments=[],tuple_hypotheses=[],
            sequence_component_hypotheses=[],source_like_component_hypotheses=[],query_activity='not_tested',in_vivo_direction='unknown',scope='reference_assay_conditions_only',rule_version=VERSION))
    for row in inherited:
        qid=row['complex_anchor']
        if qid not in original:raise ValueError('Assigned reaction has no original anchor')
        if qid not in byanchor:continue
        e=entry(qid);a=original[qid];cid=row['function_id'];c=bycontract.get(cid)
        if c is None or cid not in set(filter(None,a.get('potential_function_ids','').split(';'))):raise ValueError('Unassigned or unknown reference reaction contract')
        expected=reviews[cid]
        if (row.get('status') not in {'homology_supported_potential','reference_sequence_supported_potential'} or row['status']!=a.get('function_status') or
            row.get('structure_id')!=a['structure_id'] or c['structure_id']!=a['structure_id'] or row.get('source')!=c['primary_source'] or
            any(row.get(k)!=expected[k] for k in ['molecular_reaction','reference_evidence_scope','capability_scope']) or
            row.get('query_net_direction') not in {'unresolved','unknown'}):raise ValueError('Reference reaction record differs from its bound assignment/review')
        if any(r['contract_id']==cid for r in e['inherited_assignments']):raise ValueError('Repeated assigned reference reaction contract')
        e['inherited_assignments'].append(dict(contract_id=cid,status=row['status'],reference_direction=row.get('reference_direction','')))
        if identities.get(qid) in ACCEPTED:
            e['status']='inherited_reference_potential';e['capabilities'].append(capability(c))
        else:e['blockers'].append('catalytic_identity_support_insufficient')
    for qid,alternatives in candidates.items():
        e=entry(qid)
        if identities.get(qid) not in ACCEPTED:e['blockers'].append('catalytic_identity_support_insufficient');continue
        for result in alternatives:
            sid=result['structure_id']
            if sid not in admitted:
                e['tuple_hypotheses'].append(dict(structure_id=sid,blockers=['assay_to_sequence_binding_not_admitted_'+sid]));continue
            mode=byanchor[qid]['input_mode']
            if mode!=run['input_mode']:raise ValueError('Reference tuple input mode changed')
            rs=[r for r in refs if r['structure_id']==sid]
            member_groups=bysystem[sid].get('local_member_groups')
            local_locus=locus_matches.get(qid,{}).get(sid) or None
            if result.get('entry_route')=='exact_coding_locus_lookup' and local_locus is None:
                e['tuple_hypotheses'].append(dict(structure_id=sid,blockers=['unbound_coding_locus_evidence'],matched_members=[]))
                continue
            errors,matched=tuple_proof(qid,result,rs,genes,sequences,mode,run['thresholds'],coding_proofs,source_inputs,
                coding_locus=local_locus,member_groups=member_groups)
            e['tuple_hypotheses'].append(dict(structure_id=sid,blockers=errors,matched_members=matched,
                                               **({'declared_member_groups':member_groups} if member_groups else {}),
                                               **({'member_selection_ambiguities':result.get('evidence',{}).get('member_selection_ambiguities')} if result.get('evidence',{}).get('member_selection_ambiguities') else {}),
                                               **({'member_selection_rejections':result.get('evidence',{}).get('member_selection_rejections')} if result.get('evidence',{}).get('member_selection_rejections') else {})))
            if errors:
                component_roles=result.get('evidence',{}).get('component_roles',{})
                component_locus=result.get('evidence',{}).get('component_locus',{})
                component=None
                if component_roles and all(component_roles.get(ref['role']) for ref in rs):
                    component=_component_hypothesis(qid,result,rs,component_roles,genes,sequences,mode,run['thresholds'],
                                                    coding_proofs,source_inputs,component_locus,member_groups)
                if component:
                    _selected_result,linkage_blockers,component_matched,component_ambiguities=component
                    hard=[]
                    caps=[bycontract[cid] for cid in model['admitted_systems'][sid]]
                    e['sequence_component_hypotheses'].append(dict(
                        status='source_equivalent_components_potential',structure_id=sid,
                        matched_members=component_matched,
                        capabilities=[dict(capability(c),scope='sequence_component_function_hypothesis') for c in caps],
                        linkage_blockers=linkage_blockers,strict_tuple_blockers=errors,hard_blockers=hard,
                        query_activity='not_tested',in_vivo_direction='unknown',physical_assembly='unknown',
                        claim_scope='source_equivalent_components_function_hypothesis',
                        member_selection_ambiguities=component_ambiguities))
                else:
                    e['sequence_component_hypotheses'].append(dict(status='not_supported',structure_id=sid,
                        matched_members=[],capabilities=[],linkage_blockers=[],hard_blockers=errors,
                        member_selection_ambiguities=[]))
                e['source_like_component_hypotheses'].extend(_source_module_prior(qid,result,rs,component_roles,
                    genes,sequences,mode,run['thresholds'],pyridine_results,
                    [bycontract[cid] for cid in model['admitted_systems'][sid]],coding_proofs,source_inputs,
                    component_locus,member_groups))
                # A source-class Fdh homolog can support only a coarse formate
                # role when the non-Fdh source backbone is otherwise proven.
                if module_role_inventory:
                    fdh_classes=module_role_inventory.get('fdh_source_member_classes',{}) if isinstance(module_role_inventory,dict) else {}
                    fdh_refs=[r for r in rs if r['kind']=='module' and fdh_classes.get(r.get('sequence_sha256'))=='Fdh_molybdo']
                    if len(fdh_refs)==1:
                        fdh_ref=fdh_refs[0]; backbone_refs=[r for r in rs if r is not fdh_ref]
                        backbone_groups=None if member_groups is None else [[role for role in group if role!=fdh_ref['role']] for group in member_groups]
                        if backbone_groups is not None: backbone_groups=[group for group in backbone_groups if group]
                        backbone_roles={role:component_roles.get(role,[]) for r in backbone_refs for role in [r['role']]}
                        backbone_result={'entry_route':result.get('entry_route'),'evidence':{'roles':backbone_roles,'resource_conflict':result.get('evidence',{}).get('resource_conflict',False)}}
                        backbone=_component_hypothesis(qid,backbone_result,backbone_refs,backbone_roles,genes,sequences,mode,run['thresholds'],coding_proofs,source_inputs,component_locus,backbone_groups)
                        if backbone:
                            _,backbone_linkage,backbone_matched,backbone_ambiguities=backbone
                            used={m['protein_id'] for m in backbone_matched}
                            fdh,fdh_ambiguities,_=_source_like_fdh(qid,fdh_ref,module_role_inventory,genes,sequences,mode,genes[qid],used,coding_proofs,source_inputs)
                            if fdh:
                                caps=[bycontract[cid] for cid in model['admitted_systems'][sid]]
                                confidence='low' if fdh.get('fdh_chemistry',{}).get('confidence')=='low' else 'medium'
                                e['source_like_component_hypotheses'].append(dict(status='source_like_formate_component_potential',structure_id=sid,
                                    matched_members=backbone_matched+[fdh],homologous_roles=[fdh_ref['role']],supported_functions=['candidate_H2_formate_coupling'],
                                    hard_blockers=[],linkage_risks=backbone_linkage,member_selection_ambiguities=backbone_ambiguities+fdh_ambiguities,
                                    confidence=confidence,scope='source_like_FDH_role_combination',query_activity='not_tested',in_vivo_direction='unknown',
                                    reference_reaction_details_transferred=False,role_proof=fdh,capabilities=[]))
                continue
            caps=[bycontract[cid] for cid in model['admitted_systems'][sid]]
            if set(model['admitted_systems'][sid])!=set(bysystem[sid]['reference_contract_ids']) or any(c['structure_id']!=sid for c in caps):raise ValueError('Reference tuple contract scope changed')
            e.update(status='exact_reference_coding_locus_potential' if local_locus is not None else 'exact_reference_sequence_potential',capabilities=[capability(c) for c in caps],tuple_matches=matched,matched_structure_id=sid,
                      claim_scope='gene_associated_encoded_reference_CDS_potential' if local_locus is not None else 'reference_assay_conditions_only',
                      **({'declared_member_groups':member_groups} if member_groups else {}),
                      **({'member_selection_ambiguities':result.get('evidence',{}).get('member_selection_ambiguities')} if result.get('evidence',{}).get('member_selection_ambiguities') else {}))
    exact=[(q,e) for q,e in results.items() if e['status'] in TUPLE_POTENTIAL_STATUSES]
    for q,e in exact:
        owned={m['protein_id'] for m in e['tuple_matches']}
        rivals=[other for other,r in exact if other!=q and owned & {m['protein_id'] for m in r['tuple_matches']}]
        if rivals:
            inherited_caps=[capability(bycontract[r['contract_id']]) for r in e['inherited_assignments']]
            e.update(status='inherited_reference_potential' if inherited_caps else 'not_supported',capabilities=inherited_caps)
            e['blockers'].append('exact_tuple_resources_shared_with_other_anchor')
    def component_resources(h):
        return {m['protein_id'] for m in h.get('matched_members',[])} | {
            pid for a in h.get('member_selection_ambiguities',[]) for pid in a.get('equally_supported_loci',[])
        }
    witnesses=[]
    for q,e in results.items():
        if e['status'] in TUPLE_POTENTIAL_STATUSES:
            witnesses.append((q,{m['protein_id'] for m in e['tuple_matches']}))
        for h in e.get('sequence_component_hypotheses',[]):
            if h.get('status')=='source_equivalent_components_potential':
                witnesses.append((q,component_resources(h)))
        for h in e.get('source_like_component_hypotheses',[]):
            if h.get('status')=='source_like_formate_component_potential':
                witnesses.append((q,component_resources(h)))
    for q,e in results.items():
        for h in e.get('sequence_component_hypotheses',[]):
            if h.get('status')!='source_equivalent_components_potential':continue
            owned=component_resources(h)
            if any(other!=q and owned & other_owned for other,other_owned in witnesses):
                h.update(status='not_supported',capabilities=[],hard_blockers=sorted(set(h.get('hard_blockers',[])+['component_resources_shared_with_other_anchor'])))
        for h in e.get('source_like_component_hypotheses',[]):
            if h.get('status')!='source_like_formate_component_potential':continue
            owned=component_resources(h)
            if any(other!=q and owned & other_owned for other,other_owned in witnesses):
                h.update(status='not_supported',matched_members=[],hard_blockers=sorted(set(h.get('hard_blockers',[])+['component_resources_shared_with_other_anchor'])))
    for e in results.values():
        if e['status']=='not_supported':e['blockers'].extend(b for h in e['tuple_hypotheses'] for b in h['blockers'])
        e['blockers']=sorted(set(e['blockers']))
    return results,bound
