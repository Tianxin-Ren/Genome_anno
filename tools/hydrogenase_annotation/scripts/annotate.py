#!/usr/bin/env python3
"""Conservative, sequence-based hydrogenase seed annotation with module evidence."""
import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from circular_coordinates import circular_gaps, interval_segments
from circular_gene_calls import call_prodigal, protein_call, recover_circular_calls
from discovery_evidence import FIELDS as DISCOVERY_FIELDS, search_discovery, stage_fingerprint
from input_translation import cds_coordinates, evidence as translation_evidence, normalized_translation, verified_origin_translation
from nucleotide_input import fasta_records, materialize_fasta

ROOT = Path(__file__).resolve().parents[1]
GENE_FIELDS = ["id","contig","start","end","strand","gene","aliases","protein_accession","partial","contig_length","bin_id","coordinate_evidence","topology"]
BLAST_FIELDS = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qseq sseq".split()
COMPLEX_FIELDS=('protein_id contig start end strand structure_id catalytic_family structure_confidence function_confidence '
    'potential_function_ids catalytic_status structure_status function_status missing_roles net_direction flags input_mode bin_id '
    'reference_clade clade_method reference_hyddb_label reference_ko_label hyddb_prediction ko_prediction diagnostic_positions_conserved '
    'diagnostic_residue_map hmm_support catalytic_identity catalytic_query_coverage catalytic_reference_coverage local_cluster_complete '
    'missing_context_roles small_subunit_topology potential_direct_partner reference_mechanism reference_energy_coupling '
    'in_vivo_physiology operon_status evidence_url coordinate_evidence active_center_status active_center_template '
    'active_center_span_status catalytic_domain_integrity module_graph_search_status module_graph_assignment_states localization_status '
    'module_organization_status missing_optional_roles structure_completeness_scope cross_anchor_pairing_status '
    'cross_anchor_competing_anchors cross_anchor_partner_proteins cross_anchor_function_ambiguity_ids physical_interaction_status '
    'split_center_status split_center_component_proteins identity_confirmation_status catalytic_identity_confidence identity_confirmation_family').split()
MODULE_FIELDS=('anchor structure_id role protein_id scope bin_id identity query_coverage reference_coverage rule_query_coverage '
    'coverage_scope query_domain_start query_domain_end reference_source_region_start reference_source_region_end domain_evidence_url '
    'reference_cysteines_conserved reference_protein contig start end coordinate_evidence module_node_id module_category '
    'node_reference_start node_reference_end node_coverage node_aligned_pair_identity node_realization required_for_structure localization_status '
    'cross_anchor_pairing_status cross_anchor_competing_anchors').split()
REACTION_FIELDS='complex_anchor function_id structure_id status reference_reaction reference_direction query_net_direction source molecular_reaction reference_evidence_scope capability_scope reference_context query_physiology limitations association_pairing_status physical_interaction_status'.split()


def software_versions():
    versions=dict(python=dict(version=sys.version,executable=sys.executable),biopython=dict(version=Bio.__version__))
    for name,argument in [('blastp','-version'),('hmmsearch','-h'),('prodigal','-v'),('mafft','--version'),('iqtree3','--version')]:
        executable=shutil.which(name)
        if executable:
            result=subprocess.run([executable,argument],capture_output=True,text=True,timeout=10)
            versions[name]=dict(executable=executable,version_output='\n'.join((result.stdout+result.stderr).strip().splitlines()[:6]),exit_code=result.returncode)
        else:
            versions[name]=dict(executable=None,version_output='not_available_not_used_by_this_run')
    return versions


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else [])
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def extract_genbank(path, topology_override=None):
    genes, contigs = [], {}
    with Path(path).open() as handle:
        for record in SeqIO.parse(handle, "genbank"):
            contigs[record.id] = len(record)
            topology = topology_override or record.annotations.get("topology")
            for feature in record.features:
                q = feature.qualifiers
                if feature.type != "CDS" or "translation" not in q or "pseudo" in q or "pseudogene" in q:
                    continue
                gid = q.get("locus_tag", q.get("protein_id", [""]))[0]
                if not gid:
                    raise ValueError("Translated CDS requires a locus_tag or protein_id")
                sequence, _sec_edits = normalized_translation(record, feature)
                coordinate_evidence = "genomic_CDS"
                try:
                    start, end, strand, partial = cds_coordinates(record, feature, topology)
                    if end > len(record) and not verified_origin_translation(record, feature, sequence):
                        raise ValueError("Origin-spanning CDS translation is not verified")
                except ValueError:
                    # Preserve the supplied protein without treating a gapped or unverified join as a genomic span.
                    start = end = strand = 0
                    partial = False
                    coordinate_evidence = "unavailable_noncontiguous_or_unverified_CDS"
                gene = dict(id=gid, contig=record.id, start=start,
                    end=end, strand=strand,
                    gene=q.get("gene",[""])[0], aliases=";".join(q.get("old_locus_tag",[])),
                    protein_accession=q.get("protein_id",[""])[0], partial=partial,
                    contig_length=len(record), sequence=sequence,
                    topology=topology if topology in {"linear", "circular"} else None)
                if coordinate_evidence != "genomic_CDS":
                    gene["coordinate_evidence"] = coordinate_evidence
                genes.append(gene)
    if not contigs:
        raise ValueError(f"No GenBank records in {path}")
    counts=Counter(g["id"] for g in genes)
    for gene in genes:
        if counts[gene["id"]]>1:
            original=gene["id"]
            if not gene["protein_accession"]:
                raise ValueError(f"Repeated GenBank locus {original} has no protein accession")
            gene["aliases"]=";".join(filter(None,[gene["aliases"],original]))
            gene["id"]=original+"|"+gene["protein_accession"]
    return genes, contigs


def read_faa_gff(faa, gff, topology_override=None):
    sequences, aliases = {}, {}
    with Path(faa).open() as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in sequences:
                raise ValueError(f"Duplicate protein ID: {record.id}")
            sequences[record.id] = str(record.seq).rstrip("*")
            internal=re.search(r"(?:^|[; ])ID=([^; ]+)",record.description)
            for key in [record.id]+([internal.group(1)] if internal else []):
                if key in aliases and aliases[key]!=record.id:
                    raise ValueError(f"Ambiguous FASTA/GFF ID alias: {key}")
                aliases[key]=record.id
    coords, lengths, topologies = {}, {}, {}
    with Path(gff).open() as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                _, contig, start, end = line.split()
                lengths[contig] = int(end)
            if line.lower().startswith("##topology="):
                value = line.split("=", 1)[1].strip().lower()
                if value in {"linear", "circular"}:
                    topologies["*"] = value
            if line.lower().startswith("##topology "):
                _, contig, value = line.split(None, 2)
                value = value.strip().lower()
                if value in {"linear", "circular"}:
                    topologies[contig] = value
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields)!=9:
                raise ValueError("GFF requires nine tab-separated fields")
            attrs = dict(x.split("=",1) for x in fields[8].split(";") if "=" in x)
            if fields[2] == "region" and fields[3] == "1" and attrs.get("Is_circular") == "true":
                lengths.setdefault(fields[0], int(fields[4]))
                topologies.setdefault(fields[0], "circular")
            if fields[2]!="CDS":
                continue
            keys = [unquote(attrs.get(k,"")) for k in ["protein_id","ID","Parent","locus_tag"]]
            matched = list(dict.fromkeys(aliases[k] for k in keys if k in aliases))
            if len(matched)>1:
                raise ValueError(f"Ambiguous GFF/protein mapping: {matched}")
            if not matched:
                continue
            gid = matched[0]
            if gid in coords:
                raise ValueError(f"Multiple CDS segments for {gid}; supply a prokaryotic, unsplit CDS GFF")
            start,end = int(fields[3]),int(fields[4])
            if start<1 or end<start or fields[6] not in ["+","-"]:
                raise ValueError(f"Invalid coordinates/strand for {gid}")
            partial=attrs.get("partial","00")!="00" or attrs.get("start_range","")!="" or attrs.get("end_range","")!=""
            coords[gid] = dict(id=gid,contig=fields[0],start=start,end=end,strand=1 if fields[6]=="+" else -1,
                gene="",aliases="",protein_accession="",partial=partial,contig_length=lengths.get(fields[0],0),sequence=sequences[gid],
                topology=topologies.get(fields[0], topologies.get("*")))
    missing = set(sequences)-set(coords)
    if missing:
        raise ValueError(f"No unambiguous GFF CDS for {len(missing)} proteins, e.g. {sorted(missing)[:5]}")
    for gene in coords.values():
        gene["topology"] = topology_override or topologies.get(gene["contig"], topologies.get("*"))
        gene["contig_length"] = lengths.get(gene["contig"], 0)
        length = gene["contig_length"]
        if length and gene["end"] > length:
            if gene["topology"] != "circular":
                raise ValueError("CDS exceeds GFF sequence-region length")
            interval_segments(gene["start"], gene["end"], length)
    return list(coords.values()), lengths


def validate_genes(genes):
    if not genes:
        raise ValueError("No translated CDS proteins found")
    if len({g["id"] for g in genes})!=len(genes):
        raise ValueError("Protein IDs must be unique across the input genome/MAG")
    for gene in genes:
        seq=gene["sequence"].upper().rstrip("*")
        if not seq or re.search(r"[^ACDEFGHIKLMNPQRSTVWYUBZOXJ]",seq):
            raise ValueError(f"Invalid or internally stopped protein: {gene['id']}")
        gene["sequence"]=seq


def repeat_indel_ambiguity(hit):
    """Reference positions with equally scoring shifts of an exact repeat indel."""
    query,reference=hit['qseq'],hit['sseq']
    if len(query)!=len(reference):raise ValueError('Alignment rows have different lengths')
    positions=[];position=int(hit['sstart'])-1
    for residue in reference:
        if residue!='-':position+=1
        positions.append(position if residue!='-' else None)
    ambiguous=set()
    for gapped,opposite in [(query,reference),(reference,query)]:
        index=0
        while index<len(gapped):
            if gapped[index]!='-':index+=1;continue
            start=index
            while index<len(gapped) and gapped[index]=='-':index+=1
            end=index;left,right=start,end
            while (left>0 and gapped[left-1]!='-' and opposite[left-1]!='-' and
                   opposite[right-1]!='-' and opposite[left-1]==opposite[right-1]):left-=1;right-=1
            earliest=left;left,right=start,end
            while (right<len(gapped) and gapped[right]!='-' and opposite[left]!='-' and
                   opposite[right]!='-' and opposite[left]==opposite[right]):left+=1;right+=1
            if earliest<start or right>end:
                ambiguous.update(p for p in positions[earliest:right] if p is not None)
    return ambiguous


def repeat_mapping_ambiguity(mapped, reference, query):
    """Exact homopolymer indels between observed flanks, including HSP joins."""
    ambiguous=set();start=0
    while start<len(reference):
        end=start+1
        while end<len(reference) and reference[end]==reference[start]:end+=1
        left,_=mapped.get(start,(None,''));right,_=mapped.get(end+1,(None,''))
        if left is not None and right is not None and left<right:
            if not 1<=left<right<=len(query):raise ValueError('Mapped flanks outside query sequence')
            middle=query[left:right-1]
            if middle and len(middle)!=end-start and set(middle)=={reference[start]}:
                ambiguous.update(range(start+1,end+1))
        start=end
    return ambiguous


def residue_states(hit, positions, expected="C"):
    """Map observed residues separately from absent coverage and unknown letters."""
    rpos,qpos=int(hit["sstart"])-1,int(hit["qstart"])-1
    mapped={}
    for q,r in zip(hit["qseq"],hit["sseq"]):
        if r!="-":
            rpos+=1
        if q!="-":
            qpos+=1
        if r!="-" and rpos in positions:
            mapped[rpos]=(qpos if q!="-" else None,q)
    ambiguous_positions=repeat_indel_ambiguity(hit)
    states=[]
    for pos in positions:
        allowed=expected[str(pos)] if isinstance(expected,dict) else expected
        query_pos,observed=mapped.get(pos,(None,""))
        status=("not_covered" if observed in ["", "-"] else "ambiguous" if observed in "XBZJO" else
                "present" if observed in allowed else "substituted")
        if pos in ambiguous_positions:status="ambiguous";query_pos=None
        states.append(dict(reference_position=pos,expected_residue=allowed,query_position=query_pos,
            observed_residue=observed,status=status,coverage_reason="repeat_indel_alignment_ambiguity" if pos in ambiguous_positions else "alignment_deletion" if observed=="-" else
            "outside_alignment" if observed=="" else "aligned"))
    return states,mapped


def conserved_positions(hit, positions, expected="C"):
    states,mapped=residue_states(hit,positions,expected)
    return all(s["status"]=="present" for s in states),mapped


def read_bin_map(path, genes):
    """Require explicit contig-to-bin ownership; '.' means unbinned, never shared."""
    with Path(path).open() as handle:
        reader=csv.DictReader(handle,delimiter="\t")
        if reader.fieldnames != ["contig", "bin_id"]:
            raise ValueError("Bin map must have exactly contig and bin_id TSV columns")
        mapping={}
        for row in reader:
            contig,bin_id=row["contig"],row["bin_id"]
            if not contig or contig in mapping or bin_id is None:
                raise ValueError(f"Missing or duplicate contig in bin map: {contig}")
            mapping[contig]="" if bin_id in ["", "."] else bin_id
    contigs={g["contig"] for g in genes}
    if contigs-set(mapping):
        raise ValueError(f"Contigs missing from bin map: {sorted(contigs-set(mapping))[:5]}")
    for gene in genes:
        gene["bin_id"]=mapping[gene["contig"]]
    return mapping


def association_allowed(anchor, partner, mode):
    if any(g.get("coordinate_evidence") == "unavailable_noncontiguous_or_unverified_CDS" for g in (anchor, partner)):
        return anchor["id"] == partner["id"]
    if mode == "genome":
        return True
    if mode == "bin":
        if anchor.get("bin_id") or partner.get("bin_id"):
            return bool(anchor.get("bin_id")) and anchor.get("bin_id")==partner.get("bin_id")
        return anchor["contig"]==partner["contig"]
    if mode == "assembly":
        return anchor["contig"]==partner["contig"]
    return anchor["id"]==partner["id"]


def load_blast(path, thresholds, references=None):
    refs={r["reference_id"]:r for r in references or []}
    hits=[]
    for line in Path(path).read_text().splitlines():
        h=dict(zip(BLAST_FIELDS,line.split("\t")))
        for key in ["pident","evalue","bitscore"]:
            h[key]=float(h[key])
        for key in ["length","qstart","qend","sstart","send","qlen","slen"]:
            h[key]=int(h[key])
        h["qcov"]=(h["qend"]-h["qstart"]+1)/h["qlen"]
        h["scov"]=(h["send"]-h["sstart"]+1)/h["slen"]
        h["coverage_scope"]=refs.get(h["sseqid"],{}).get("coverage_scope","full_reference_protein")
        # A reviewed domain can occur in a longer protein; retain raw whole-query coverage separately.
        h["rule_query_coverage"]=(min(1,(h["qend"]-h["qstart"]+1)/min(h["qlen"],h["slen"]))
            if h["coverage_scope"]=="reviewed_reference_domain" else h["qcov"])
        h["passes_thresholds"]=(h["evalue"]<=thresholds["blast_evalue"] and h["pident"]>=thresholds["min_identity"] and
            h["rule_query_coverage"]>=thresholds["min_query_coverage"] and h["scov"]>=thresholds["min_reference_coverage"])
        hits.append(h)
    return hits


def _circular_contig_length(genes):
    """Return a documented circular contig length, or ``None``.

    Older tables do not have a topology column, but the source contig IDs used
    by the natural-assembly importer carry the explicit ``circular-yes``
    marker.  Accept the newer explicit fields too, while refusing to infer
    circularity from a length alone.
    """
    if not genes or len({g.get("contig") for g in genes}) != 1:
        return None
    values = []
    saw_explicit = False
    saw_linear = False
    for gene in genes:
        topology = gene.get("topology") or gene.get("contig_topology")
        circular = gene.get("circular", gene.get("is_circular"))
        contig = str(gene.get("contig", ""))
        # An explicit linear declaration is authoritative over a legacy
        # circular marker embedded in the contig name.
        if topology is not None:
            normalized = str(topology).lower()
            if normalized in {"linear", "line", "no", "false", "0"}:
                saw_linear = True
                continue
            explicit = normalized in {"circular", "circle", "yes", "true", "1"}
            if not explicit:
                return 0
        elif circular is not None:
            normalized = str(circular).lower()
            if normalized in {"false", "no", "linear", "0"}:
                saw_linear = True
                continue
            explicit = circular is True or normalized in {"circular", "circle", "yes", "true", "1"}
            if not explicit:
                return 0
        else:
            explicit = bool(re.search(r"(?:^|[_-])circular-(?:yes|true|1)(?:[_-]|$)", contig, re.I))
        if explicit:
            saw_explicit = True
            try:
                length = int(gene.get("contig_length") or 0)
            except (TypeError, ValueError):
                return 0
            if length <= 0:
                return 0
            values.append(length)
    if saw_linear and saw_explicit:
        return 0
    if saw_explicit and (len(values) != len(genes) or len(set(values)) != 1):
        return 0
    return values[0] if saw_explicit else None


def local_cluster(genes, thresholds):
    if any(g.get("coordinate_evidence") == "unavailable_noncontiguous_or_unverified_CDS" for g in genes):
        return False
    if len({g["contig"] for g in genes})!=1 or len({g["strand"] for g in genes})!=1:
        return False
    genes=sorted(genes,key=lambda g:g["start"])
    circular_length = _circular_contig_length(genes)
    if circular_length is not None:
        if circular_length <= 0:
            return False
        # Coordinates are 1-based inclusive.  The largest unoccupied gap is
        # the cut point; the remaining arc is the shortest circular cluster.
        try:
            gaps = circular_gaps([(g["start"], g["end"]) for g in genes], circular_length)
        except ValueError:
            return False
        cut_index = max(range(len(gaps)), key=gaps.__getitem__)
        included_gaps = [gap for index, gap in enumerate(gaps) if index != cut_index]
        occupied_span = circular_length - gaps[cut_index]
        return (occupied_span <= thresholds["local_max_span_bp"] and
                all(gap <= thresholds["local_max_intergenic_gap_bp"] for gap in included_gaps))
    return (genes[-1]["end"]-genes[0]["start"]+1<=thresholds["local_max_span_bp"] and
        all(b["start"]-a["end"]-1<=thresholds["local_max_intergenic_gap_bp"] for a,b in zip(genes,genes[1:])))


def classify(genes, hits, hmm_hits, references, systems, thresholds, mode="assembly", graph_policy=None):
    by_id={g["id"]:g for g in genes}
    refs={r["reference_id"]:r for r in references}
    by_query=defaultdict(list)
    for h in hits:
        if h["passes_thresholds"]:
            by_query[h["qseqid"]].append(h)
    for rows in by_query.values():
        rows.sort(key=lambda h:h["bitscore"],reverse=True)
    from module_graph import load_policy, resolve, export_graph
    from identity_truth import adjudications
    identity_review=adjudications(ROOT)
    if graph_policy is None:graph_policy=load_policy(ROOT,references)
    annotations,modules,rejected,graph_nodes,graph_edges=[],[],[],[],[]
    hmm_by_query=defaultdict(list)
    for h in hmm_hits:
        hmm_by_query[h["protein_id"]].append(h)
    for values in hmm_by_query.values():
        values.sort(key=lambda h:h["independent_evalue"])
    anchors=set(hmm_by_query) | {q for q,rows in by_query.items() if any(refs[h["sseqid"]]["kind"] in ["catalytic","unresolved","decoy"] for h in rows)}
    for qid in sorted(anchors):
        retained_hmd=False
        rows=by_query.get(qid,[])
        cat=[h for h in rows if refs[h["sseqid"]]["kind"]=="catalytic"]
        competing=[h for h in rows if refs[h["sseqid"]]["kind"] in ["unresolved","decoy"]]
        if competing and (not cat or competing[0]["bitscore"]*thresholds["catalytic_winner_ratio"]>=cat[0]["bitscore"]):
            best=refs[competing[0]["sseqid"]]
            if best["kind"]=="decoy":
                review=identity_review.get(best['sequence_sha256'],{})
                # Hmd paralog similarity cannot exclude H2 chemistry with a different substrate.
                retained_hmd=(review.get('family')=='Fe' and review.get('identity_label')=='indeterminate' and
                    any(h['family']=='Fe' for h in hmm_by_query[qid]))
                if not retained_hmd:
                    rejected.append(dict(protein_id=qid,reason="non_catalytic_homolog_wins",reference=best["reference_id"]))
                    continue
            cat=[]
        if not cat:
            if hmm_by_query[qid] or (competing and refs[competing[0]["sseqid"]]["kind"]=="unresolved"):
                family=('Fe' if retained_hmd else refs[competing[0]["sseqid"]]["center"] if competing else hmm_by_query[qid][0]["family"])
                gene=by_id[qid]
                annotations.append(dict(protein_id=qid,contig=gene["contig"],start=gene["start"],end=gene["end"],strand=gene["strand"],
                    structure_id="unresolved",catalytic_family=family,
                    structure_confidence="low",function_confidence="unassigned",potential_function_ids="",
                    catalytic_status="candidate",structure_status="unresolved",function_status="unassigned",
                    hmm_support=';'.join(sorted({h.get('profile',h['family']) for h in hmm_by_query[qid]})),
                    evidence_url=refs[competing[0]['sseqid']].get('evidence_url','') if competing else str(ROOT/'references/processed/profile_manifest.json'),
                    missing_roles="no_reviewed_structure_rule",net_direction="unresolved",flags="family_candidate_only"+
                    (';Hmd_paralog_competitor_blocks_function_transfer_not_discovery' if retained_hmd else '')))
            continue
        best=cat[0]
        ref=refs[best["sseqid"]]
        sid=ref["structure_id"]
        system=next(s for s in systems if s["structure_id"]==sid)
        alternate=cat[1]["bitscore"] if len(cat)>1 else 0
        distinct=not alternate or best["bitscore"]>=alternate*thresholds["catalytic_winner_ratio"]
        motifs, mapping=conserved_positions(best,ref["diagnostic_positions"],expected=ref.get("diagnostic_residues","C"))
        anchor=by_id[qid]
        context=system.get("context",{})
        graph=resolve(anchor,best,system,hits,by_id,thresholds,mode,graph_policy,association_allowed,local_cluster)
        node_rows,edge_rows=export_graph(anchor,sid,graph,by_id);graph_nodes.extend(node_rows);graph_edges.extend(edge_rows)
        chosen=graph['chosen']
        missing=sorted({n['role'] for n in graph['nodes'] if n['node_id'] in graph['missing_function']})
        ambiguous=['module_graph_alternative_assignment'] if graph['ambiguous'] else []
        local_ok=not any(n['required_for_structure'] and n['linkage']=='local' and n['node_id'] in graph['missing_structure'] for n in graph['nodes'])
        complete=local_ok and not graph['missing_structure'] and not graph['structural_ambiguous'] and not graph['search_truncated'] and distinct and motifs
        structural_nodes=[graph['selected'][n['node_id']] for n in graph['nodes'] if n['required_for_structure'] and n['node_id'] in graph['selected']]
        partial=any(by_id[h['qseqid']]['partial'] for n in graph['nodes'] if n['required_for_structure'] and n['node_id'] in graph['selected'] for h in [p['hit'] for p in graph['selected'][n['node_id']]['parts']])
        high_similarity=all(o['identity']>=thresholds['high_structure_identity'] and o['coverage']>=thresholds['high_structure_coverage'] for o in structural_nodes)
        cys_ok=True;structure_cys_ok=True;sec_failures=[]
        for role,hs in chosen.items():
            r=refs[hs[0]['sseqid']]
            pieces=[p for n in graph['nodes'] if n['role']==role and n['node_id'] in graph['selected'] for p in graph['selected'][n['node_id']]['parts']]
            def conserved_in_assigned_resources(positions,expected='C'):
                states=[s for p in pieces for s in residue_states(p['hit'],positions,expected)[0] if s['reference_position'] in p['reference_positions']]
                return all(any(s['reference_position']==pos and s['status']=='present' for s in states) for pos in positions)
            keep=conserved_in_assigned_resources(r['cysteine_positions'])
            if role not in context:cys_ok &= keep
            if any(n['role']==role and n['required_for_structure'] for n in graph['nodes']):structure_cys_ok &= keep
            sec_ok=conserved_in_assigned_resources(r.get('selenocysteine_positions',[]),'U')
            if not sec_ok:sec_failures.append(role)
            for n in graph['nodes']:
                if n['role']!=role or n['node_id'] not in graph['selected']:continue
                option=graph['selected'][n['node_id']]
                for piece in option['parts']:
                    h=piece['hit'];g=by_id[h['qseqid']]
                    modules.append(dict(anchor=qid,structure_id=sid,role=role,protein_id=h['qseqid'],
                        scope="local" if role in system["local"] else "same_genome_context" if mode=="genome" and role in context else
                            "same_genome_remote" if mode=="genome" else "same_bin_context" if mode=="bin" and anchor.get("bin_id") and role in context else
                            "same_bin_remote" if mode=="bin" and anchor.get("bin_id") else "same_contig",
                        bin_id=g.get('bin_id',''),identity=h['pident'],query_coverage=h['qcov'],reference_coverage=h['scov'],
                        rule_query_coverage=h.get('rule_query_coverage',h['qcov']),coverage_scope='module_graph_reference_segment',
                        query_domain_start=min(piece['query_positions']),query_domain_end=max(piece['query_positions']),
                        reference_source_region_start=r.get('protein_region_start',1),reference_source_region_end=r.get('protein_region_end',h['slen']),
                        domain_evidence_url=r.get('domain_evidence_url',''),reference_cysteines_conserved=keep,reference_protein=r['protein_accession'],
                        contig=g['contig'],start=g['start'],end=g['end'],module_node_id=n['node_id'],module_category=n['module_category'],
                        node_reference_start=n['reference_start'],node_reference_end=n['reference_end'],node_coverage=option['coverage'],
                        node_aligned_pair_identity=piece['identity'],
                        node_realization='collinear_CDS_fragments' if option['fragmented'] else 'single_polypeptide_segment',
                        required_for_structure=n['required_for_structure'],localization_status=graph['localization_status']))
        flags=[]
        if not distinct:flags.append("catalytic_reference_ambiguous")
        if not motifs:flags.append("diagnostic_evidence_incomplete")
        diagnostic_states,_=residue_states(best,ref["diagnostic_positions"],ref.get("diagnostic_residues","C"))
        if any(s["status"]=="not_covered" for s in diagnostic_states):flags.append("diagnostic_site_not_covered")
        if any(s["status"]=="substituted" for s in diagnostic_states):flags.append("diagnostic_site_substituted")
        if any(s["status"]=="ambiguous" for s in diagnostic_states):flags.append("diagnostic_site_ambiguous")
        if not cys_ok:flags.append("module_cysteine_topology_not_conserved")
        if sec_failures:flags.append("selenocysteine_or_recoding_requires_review:"+",".join(sec_failures))
        if not local_ok:flags.append("local_gene_cluster_incomplete_or_disrupted")
        if partial:flags.append("partial_CDS")
        if ambiguous:flags.append("ambiguous_modules:"+",".join(ambiguous))
        if system["remote"]:flags.append("remote_association_inferred_not_operonic")
        if context:flags.append("physiological_context_is_not_an_enzyme_subunit")
        below_transfer=[n['node_id'] for n in graph['nodes'] if n['required_for_function'] and n['node_id'] in graph['selected'] and
            (graph['selected'][n['node_id']]['identity']<thresholds['function_transfer_identity'] or graph['selected'][n['node_id']]['coverage']<.8)]
        if graph['fragmented']:flags.append('module_CDS_fragments_natural_fission_or_assembly_error_unresolved')
        if graph['fused_queries']:flags.append('fusion_relative_to_reference_partition_function_unvalidated')
        if graph['search_truncated']:flags.append('module_graph_search_limit_no_confident_assignment')
        if below_transfer:flags.append("module_below_function_transfer_threshold:"+",".join(below_transfer))
        flags.append("net_direction_and_kinetics_not_measured")
        structure="high" if complete and high_similarity and structure_cys_ok and not partial and not graph['organization_changed'] else "medium" if local_ok and distinct and motifs and not partial else "low"
        # shortcut: a small seed set cannot establish novel mechanisms; require close similarity across all modules and expand with independent biochemical references before relaxing.
        transferable=complete and not missing and not ambiguous and not graph['organization_changed'] and cys_ok and not sec_failures and not partial and not below_transfer and not any(by_id[h["qseqid"]]["partial"] for hs in chosen.values() for h in hs)
        missing_accessory=[n['node_id'] for n in graph['nodes'] if n['module_category']=='assembly_accessory' and n['node_id'] not in graph['selected']]
        annotations.append(dict(protein_id=qid,contig=anchor["contig"],start=anchor["start"],end=anchor["end"],strand=anchor["strand"],
            structure_id=sid,catalytic_family=system["center"],structure_confidence=structure,
            function_confidence="medium_homology_potential" if transferable else "unassigned",
            catalytic_status="supported" if distinct and motifs and not anchor["partial"] else "candidate",
            structure_status="composition_supported_reorganized" if complete and graph['organization_changed'] else "core_complete_accessory_missing" if complete and missing_accessory else "complete_supported" if complete and not partial else "incomplete_or_ambiguous",
            function_status="homology_supported_potential" if transferable else "unassigned",
            input_mode=mode,bin_id=anchor.get("bin_id",""),
            potential_function_ids=";".join(system["function_ids"]) if transferable else "",
            reference_clade=system["reference_clade"],clade_method="nearest_reference; not a phylogenetic placement",
            reference_hyddb_label=";".join(sorted({r["HydDB对照（非定义）"] for r in system["reference_function_records"]})),
            reference_ko_label=";".join(sorted({r["KO已独立核实（未核实留空）"] for r in system["reference_function_records"] if r["KO已独立核实（未核实留空）"]})),
            hyddb_prediction="not_run",ko_prediction="not_run",diagnostic_positions_conserved=motifs,
            diagnostic_residue_map=json.dumps(mapping),hmm_support=";".join(sorted({x.get("profile",x["family"]) for x in hmm_by_query[qid]})),
            catalytic_identity=best["pident"],catalytic_query_coverage=best["qcov"],catalytic_reference_coverage=best["scov"],
            local_cluster_complete=local_ok,missing_roles=";".join(missing),
            missing_context_roles=";".join(r for r in missing if r in context),
            small_subunit_topology=("homology_prediction: "+system["reference_small_subunit_topology"]) if structure=="high" else "unresolved",
            potential_direct_partner=system["direct_partner"] if transferable else "unresolved",
            reference_mechanism=system["mechanism"],reference_energy_coupling=system["reference_energy_coupling"],
            net_direction="unresolved",in_vivo_physiology="unresolved",operon_status="predicted_gene_cluster; transcription_not_tested",
            module_graph_search_status="truncated" if graph["search_truncated"] else "complete",module_graph_assignment_states=graph["assignment_states"],
            module_organization_status='changed_from_reference_unvalidated' if graph['organization_changed'] else 'reference_partition_compatible',
            missing_optional_roles=';'.join(missing_accessory),structure_completeness_scope='configured_required_modules; physical interaction and oligomer assembly unresolved',
            localization_status=graph["localization_status"],evidence_url=system["evidence_url"],flags=";".join(flags)))
    return annotations,modules,rejected,graph_nodes,graph_edges


def annotate(genes, out, threads=2, input_files=(), mode="assembly", translation_table=None):
    if mode not in ["genome", "bin", "assembly", "transcript"]:
        raise ValueError(f"Unsupported input mode: {mode}")
    if mode=="bin" and any("bin_id" not in g for g in genes):
        raise ValueError("Bin mode requires explicit bin_id ownership for every protein")
    validate_genes(genes)
    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    (out/"run_manifest.json").write_text(json.dumps(dict(status="running"))+"\n")
    reference_manifest=json.loads((ROOT/"references/processed/build_manifest.json").read_text())
    if reference_manifest.get("status")!="complete":
        raise ValueError("Reference build is incomplete; rerun build_reference.py before annotation")
    from audit_reference_domains import verify_applied_review
    reference_role_audit=verify_applied_review(ROOT,reference_manifest)
    if reference_manifest.get("reference_domains_sha256")!=hashlib.sha256((ROOT/"config/reference_domains.json").read_bytes()).hexdigest():
        raise ValueError("Reference domain review changed; rebuild reference before annotation")
    if (reference_manifest.get('branch_profile_config_sha256')!=hashlib.sha256((ROOT/'config/discovery_profiles.json').read_bytes()).hexdigest() or
        reference_manifest.get('profile_manifest_sha256')!=hashlib.sha256((ROOT/'references/processed/profile_manifest.json').read_bytes()).hexdigest()):
        raise ValueError('Discovery profile configuration or manifest changed; rebuild reference before annotation')
    fes_model_path=ROOT/'references/processed/fes/manifest.json'
    if not fes_model_path.exists() or reference_manifest.get('fes_model_manifest_sha256')!=hashlib.sha256(fes_model_path.read_bytes()).hexdigest():
        raise ValueError('FeS reference build is missing or changed; rebuild reference before annotation')
    thresholds=json.loads((ROOT/"config/thresholds.json").read_text())
    references=json.loads((ROOT/"references/processed/references.json").read_text())
    refs_by_id={r["reference_id"]:r for r in references}
    systems=json.loads((ROOT/"config/systems.json").read_text())
    for gene in genes:
        gene.setdefault("coordinate_evidence","genomic_CDS")
    faa=out/"proteins.faa"
    SeqIO.write([SeqRecord(Seq(g["sequence"]),id=g["id"],description="") for g in genes],faa,"fasta")
    write_tsv(out/"genes.tsv",genes,GENE_FIELDS)
    blast=out/"blast.tsv"
    with (out/"blast.log").open("w") as log:
        subprocess.run(["blastp","-query",str(faa),"-db",str(ROOT/"references/processed/seed"),"-evalue","1e-5",
            "-max_target_seqs","1000","-max_hsps","1","-num_threads",str(threads),"-outfmt","6 "+" ".join(BLAST_FIELDS),"-out",str(blast)],stderr=log,check=True)
    hmm_hits,profile_manifest_path,discovery_policy_path=search_discovery(faa,out,ROOT,threads,thresholds)
    hits=load_blast(blast,thresholds,references)
    annotations,modules,rejected,graph_nodes,graph_edges=classify(genes,hits,hmm_hits,references,systems,thresholds,mode)
    from module_graph import NODE_FIELDS, EDGE_FIELDS, PAIR_FIELDS, apply_cross_anchor_guard
    genes_by_id={g["id"]:g for g in genes}
    for a in annotations:
        a.update(input_mode=mode,bin_id=genes_by_id[a["protein_id"]].get("bin_id",""))
        for key,value in dict(module_graph_search_status='not_applicable_no_reference_system',module_graph_assignment_states=0,
            module_organization_status='unresolved',missing_optional_roles='',structure_completeness_scope='not_established',
            localization_status='not_predicted; gene colocalization is not subcellular localization').items():a.setdefault(key,value)
        a.setdefault('in_vivo_physiology','unresolved')
        a.setdefault('coordinate_evidence',genes_by_id[a['protein_id']]['coordinate_evidence'])
        if mode=="transcript":
            a.update(start="",end="",strand="",operon_status="not_applicable_translated_transcript",
                coordinate_evidence="unavailable_translated_transcript")
    if mode=="transcript":
        for m in modules:
            m.update(start="",end="",strand="",scope="same_polypeptide",coordinate_evidence="unavailable_translated_transcript")
    from active_centers import annotate_active_centers, apply_site_guard
    active_rows,active_manifest=annotate_active_centers(faa,out,threads,anchors={a['protein_id'] for a in annotations})
    from split_centers import annotate_split_centers
    active_rows,split_manifest=annotate_split_centers(faa,out,genes,annotations,active_rows,mode,threads,ROOT)
    apply_site_guard(annotations,active_rows)
    from catalytic_confirmation import annotate_confirmation
    confirmation_manifest=annotate_confirmation(out,annotations,threads,ROOT)
    # BEGIN reviewed post-confirmation L4 selection
    from copy import deepcopy
    from complex_selection import refine, FIELDS as HYPOTHESIS_FIELDS
    selection_policy_path=ROOT/'config/complex_selection.json'
    selection_policy=json.loads(selection_policy_path.read_text())
    protected_names=['proteins.faa','genes.tsv','active_center_summary.tsv','active_center_blast.tsv','catalytic_sites.tsv','catalytic_confirmation.tsv','confirmation/features.json']
    protected_hashes={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in protected_names}
    # BEGIN reviewed auxiliary noncanonical literature evidence
    from noncanonical_center_evidence import annotate_evidence
    noncanonical_manifest=annotate_evidence(out/'proteins.faa',out/'noncanonical_centers',
        frozenset(a['protein_id'] for a in annotations if a['catalytic_family']=='FeFe'),
        frozenset(g['id'] for g in genes if g['partial']),threads=threads,root=ROOT)
    if protected_hashes!={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in protected_names}:
        raise ValueError('Auxiliary noncanonical evidence changed protected identity/site data')
    # END reviewed auxiliary noncanonical literature evidence
    from complex_extensions import search_and_merge
    references,systems,hits,selection_graph_policy,extension_manifest,reference_contracts=search_and_merge(deepcopy(genes),deepcopy(hits),deepcopy(references),deepcopy(systems),deepcopy(thresholds),mode,deepcopy(annotations),selection_policy,out/'structure_extensions',threads)
    refs_by_id={r['reference_id']:r for r in references}
    patches,replacements,hypotheses=refine(deepcopy(genes),deepcopy(hits),deepcopy(references),deepcopy(systems),deepcopy(thresholds),mode,deepcopy(annotations),selection_policy,graph_policy=selection_graph_policy,
        confirmation_features=json.loads((out/'confirmation/features.json').read_text()),
        center_candidates=list(csv.DictReader((out/'noncanonical_centers/alternatives.tsv').read_text().splitlines(),delimiter='\t')))
    allowed_fields=set(('structure_id structure_confidence function_confidence potential_function_ids structure_status function_status '
        'missing_roles reference_clade clade_method reference_hyddb_label reference_ko_label diagnostic_positions_conserved diagnostic_residue_map '
        'catalytic_identity catalytic_query_coverage catalytic_reference_coverage local_cluster_complete missing_context_roles small_subunit_topology '
        'potential_direct_partner reference_mechanism reference_energy_coupling operon_status module_graph_search_status module_graph_assignment_states '
        'module_organization_status missing_optional_roles structure_completeness_scope localization_status evidence_url flags').split())
    if set(patches)!=set(replacements) or set(patches)-{a['protein_id'] for a in annotations}:
        raise ValueError('L4 selection cannot add or remove catalytic anchors')
    for a in annotations:
        qid=a['protein_id']
        if qid not in patches:continue
        if set(patches[qid])-allowed_fields:raise ValueError('L4 selection cannot alter protected identity/site decisions')
        if patches[qid].get('potential_function_ids') or patches[qid].get('structure_confidence')=='high':
            raise ValueError('Alternative L4 selection cannot promote function or high confidence')
        new_modules,new_nodes,new_edges=replacements[qid]
        if any(r['anchor']!=qid or r['structure_id']!=patches[qid]['structure_id'] for rows in replacements[qid] for r in rows):
            raise ValueError('L4 replacement has incorrect anchor/structure ownership')
        a.update(patches[qid])
        modules[:]=[r for r in modules if r['anchor']!=qid]+new_modules
        graph_nodes[:]=[r for r in graph_nodes if r['anchor']!=qid]+new_nodes
        graph_edges[:]=[r for r in graph_edges if r['anchor']!=qid]+new_edges
    write_tsv(out/'complex_hypotheses.tsv',[{k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in h.items()} for h in hypotheses],HYPOTHESIS_FIELDS)
    selection_manifest=dict(policy_sha256=hashlib.sha256(selection_policy_path.read_bytes()).hexdigest(),script_sha256=hashlib.sha256((ROOT/'scripts/complex_selection.py').read_bytes()).hexdigest(),rescued_anchors=sorted(patches),scope='L4-only post-identity provisional core selection; no function or high-confidence promotion')
    if protected_hashes!={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in protected_names}:
        raise ValueError('Post-identity selection changed protected query/identity/site evidence')
    selection_manifest['protected_before_and_after_sha256']=protected_hashes
    selection_manifest['extensions']=extension_manifest
    selection_manifest['module_architecture']=dict(policy_sha256=hashlib.sha256((ROOT/'config/module_architecture.json').read_bytes()).hexdigest(),script_sha256=hashlib.sha256((ROOT/'scripts/module_architecture.py').read_bytes()).hexdigest(),scope='Existing FeS coordinate resources and source-reviewed compatible core competition only; no identity or function promotion')
    # END reviewed post-confirmation L4 selection
    cross_pairs=apply_cross_anchor_guard(annotations,graph_nodes,graph_edges,modules,genes_by_id,mode,
        json.loads((ROOT/'config/module_graph.json').read_text()),association_allowed)
    # BEGIN reviewed post-pairing exact-reference L5 potential
    from complex_extensions import exact_reference_patches
    exact_policy_path=ROOT/'config/exact_reference_functions.json'
    exact_policy=json.loads(exact_policy_path.read_text())
    exact_patches,exact_evidence=exact_reference_patches(deepcopy(genes),deepcopy(references),deepcopy(systems),deepcopy(annotations),deepcopy(modules),deepcopy(reference_contracts),exact_policy,deepcopy(thresholds))
    exact_allowed={'function_status','function_confidence','potential_function_ids','potential_direct_partner','flags'}
    if set(exact_patches)-{a['protein_id'] for a in annotations} or any(set(p)-exact_allowed for p in exact_patches.values()):
        raise ValueError('Exact-reference L5 rule cannot alter catalytic identity, sites, structure or anchors')
    for a in annotations:
        if a['protein_id'] in exact_patches:a.update(exact_patches[a['protein_id']])
    if protected_hashes!={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in protected_names}:
        raise ValueError('Exact-reference L5 rule changed protected identity/site evidence')
    exact_function_manifest=dict(policy_sha256=hashlib.sha256(exact_policy_path.read_bytes()).hexdigest(),script_sha256=hashlib.sha256((ROOT/'scripts/complex_extensions.py').read_bytes()).hexdigest(),assigned_anchors=sorted(exact_patches),evidence=exact_evidence,protected_before_and_after_sha256=protected_hashes,scope=exact_policy['output_scope'])
    # END reviewed post-pairing exact-reference L5 potential
    from complex_extensions import write_contract_candidates
    contract_manifest=write_contract_candidates(out,annotations,systems,reference_contracts)
    write_tsv(out/'module_graph_nodes.tsv',graph_nodes,NODE_FIELDS)
    write_tsv(out/'module_graph_edges.tsv',graph_edges,EDGE_FIELDS)
    write_tsv(out/'cross_anchor_associations.tsv',cross_pairs,PAIR_FIELDS)
    fields=list(dict.fromkeys(COMPLEX_FIELDS+[k for row in annotations for k in row]))
    write_tsv(out/"complexes.tsv",annotations,fields)
    write_tsv(out/"modules.tsv",modules,MODULE_FIELDS)
    write_tsv(out/"complex_members.tsv",modules,MODULE_FIELDS)
    sites=[]
    for a in annotations:
        best=next((h for h in hits if h["qseqid"]==a["protein_id"] and
                   refs_by_id[h["sseqid"]]["kind"]=="catalytic" and
                   refs_by_id[h["sseqid"]]["structure_id"]==a["structure_id"]),None)
        if best is None:
            continue
        ref=refs_by_id[best["sseqid"]]
        states,_=residue_states(best,ref["diagnostic_positions"],ref.get("diagnostic_residues","C"))
        for s in states:
            sites.append(dict(protein_id=a["protein_id"],reference_id=ref["reference_id"],family=ref["center"],
                **s,evidence_url=ref.get("evidence_url",""),site_model="seed_diagnostic_subset_not_complete_active_site"))
    write_tsv(out/"seed_diagnostic_sites.tsv",sites,["protein_id","reference_id","family","reference_position","expected_residue",
        "query_position","observed_residue","status","coverage_reason","evidence_url","site_model"])
    write_tsv(out/"catalytic_anchors.tsv",annotations,fields)
    write_tsv(out/"rejected.tsv",rejected,["protein_id","reason","reference"])
    write_tsv(out/"family_candidates.tsv",hmm_hits,DISCOVERY_FIELDS)
    write_tsv(out/"sequence_evidence.tsv",[{k:v for k,v in h.items() if k not in ["qseq","sseq"]} for h in hits],
              [k for k in BLAST_FIELDS if k not in ["qseq","sseq"]]+["qcov","scov","rule_query_coverage","coverage_scope","passes_thresholds"])
    from audit_maturation_domains import scan_profile_domains, architecture_supported
    domain_policy_path=ROOT/"config/maturation_domain_policy.json"
    domain_policy=json.loads(domain_policy_path.read_text())
    domain_audit_dir=ROOT/"references/processed/maturation_domain_audit"
    domain_reference=json.loads((domain_audit_dir/"summary.json").read_text())
    domain_model=domain_audit_dir/"maturation_domains.hmm"
    if (hashlib.sha256(domain_model.read_bytes()).hexdigest()!=domain_reference["combined_model_sha256"] or
            hashlib.sha256(domain_policy_path.read_bytes()).hexdigest()!=domain_reference["policy_sha256"]):
        raise ValueError("Maturation model/policy changed; rebuild reference before annotation")
    query_domains=scan_profile_domains(domain_model,faa,domain_reference["profile_models"],
        domain_reference["min_model_coverage"],out/"maturation_domains",threads)
    domain_rows=[dict(protein_id=qid,**d) for qid,ds in query_domains.items() for d in ds]
    write_tsv(out/"domain_hits.tsv",domain_rows,list(domain_rows[0]) if domain_rows else
        ["protein_id","accession","profile","domain_bitscore","query_start","query_end","passes"])
    maturation=[]
    refs={r["reference_id"]:r for r in references}
    for h in hits:
        r=refs[h["sseqid"]]
        if h["passes_thresholds"] and r["kind"]=="maturation":
            role=r["role"].split("_")[-1]
            accessions={d["accession"].split('.')[0] for d in query_domains[h["qseqid"]] if d["passes"]}
            supported=architecture_supported(accessions,domain_policy["rules"][role])
            maturation.append(dict(protein_id=h["qseqid"] ,maturation_gene=role,identity=h["pident"],reference=r["reference_id"],
                evidence_scope=r.get("evidence_scope","unreviewed_legacy_reference"),source=r.get("evidence_url",""),
                reference_domain_validation=r.get("domain_validation","not_run"),
                query_domain_validation="profile_architecture_supported" if supported else "homology_candidate_only",
                query_domain_accessions=";".join(sorted(accessions)),
                domain_specificity_scope=r.get("domain_specificity_scope",""),cofactor_assembly="unknown"))
    write_tsv(out/"maturation.tsv",maturation,["protein_id","maturation_gene","identity","reference","evidence_scope","source",
        "reference_domain_validation","query_domain_validation","query_domain_accessions","domain_specificity_scope","cofactor_assembly"])
    functions=[]
    smap={s["structure_id"]:s for s in systems}
    function_review={r['function_id']:r for r in json.loads((ROOT/'config/reference_function_review.json').read_text())['records']}
    for a in annotations:
        if not a.get("potential_function_ids"):continue
        for record in smap[a["structure_id"]]["reference_function_records"]:
            if record['功能ID'] not in a['potential_function_ids'].split(';'):continue
            functions.append(dict(complex_anchor=a["protein_id"],function_id=record["功能ID"],structure_id=a["structure_id"],
                status="homology_supported_potential",reference_reaction=record["核心反应（简化；非配平模型）"],
                reference_direction=record["H₂角色/方向"],query_net_direction="unresolved",source=smap[a["structure_id"]]["evidence_url"],
                **{k:function_review[record['功能ID']][k] for k in ['molecular_reaction','reference_evidence_scope','capability_scope','reference_context','limitations']},
                query_physiology='unresolved',association_pairing_status=a['cross_anchor_pairing_status'],physical_interaction_status=a['physical_interaction_status']))
    from complex_extensions import exact_reaction_rows
    functions.extend(exact_reaction_rows(annotations,reference_contracts))
    write_tsv(out/"potential_functions.tsv",functions,REACTION_FIELDS)
    write_tsv(out/"reaction_capabilities.tsv",functions,REACTION_FIELDS)
    from fes_sites import annotate_fes
    fes_rows,fes_manifest=annotate_fes(faa,out,modules,threads,anchors=[a['protein_id'] for a in annotations])
    evidence=[]
    anchors_by_id={a['protein_id']:a for a in annotations}
    for hit in hmm_hits:
        evidence.append(dict(protein_id=hit['protein_id'],complex_anchor=hit['protein_id'] if hit['protein_id'] in anchors_by_id else '',layer='L1',
            criterion='catalytic_profile_discovery',status='candidate_evidence',reference_id=hit['profile'],source=str(profile_manifest_path),details=json.dumps(hit),
            limitation='Discovery only. Collinear-chain E is the maximum component E, not a combined significance; no catalytic identity, full-domain integrity or function established.'))
    with (out/'catalytic_confirmation.tsv').open() as handle:
        confirmation_rows=list(csv.DictReader(handle,delimiter='\t'))
    for row in confirmation_rows:
        if row['protein_id'] not in anchors_by_id:continue
        evidence.append(dict(protein_id=row['protein_id'],complex_anchor=row['protein_id'],layer='L1',criterion='separate_catalytic_identity_confirmation',
            status=row['confirmation_status'],reference_id=row['positive_prototype'],source=str(ROOT/'references/processed/catalytic_confirmation/manifest.json'),details=json.dumps(row),
            limitation='Source-separated identity rule and explicit query-domain/ligand evidence; public-profile training overlap cannot be fully reconstructed. No molecular reaction, actual folding, occupancy or activity inferred.'))
    for site in fes_rows:
        evidence.append(dict(protein_id=site['protein_id'],complex_anchor=';'.join(site['complex_anchors']),layer='L3',
            criterion='structural_FeS_ligand_mapping',status=site['status'],reference_id=site['template_site_id'],source=site['source'],
            details=json.dumps(site),limitation='Provisional structural homology; occupancy unknown; not independent hydrogenase or function identity'))
    for s in sites:
        evidence.append(dict(protein_id=s["protein_id"],complex_anchor=s["protein_id"],layer="L2",criterion="seed_diagnostic_residue",
            status=s["status"],reference_id=s["reference_id"],source=s["evidence_url"],details=json.dumps(s),
            limitation="Diagnostic subset only, not a complete active-site model"))
    for row in active_rows:
        evidence.append(dict(protein_id=row['protein_id'],complex_anchor=row['protein_id'],layer='L2',criterion='complete_template_protein_ligand_mapping',
            status=row['active_center_status'],reference_id=row['reference_id'],source=row['sites'][0]['evidence_url'],details=json.dumps(row),
            limitation='Required sites use separately assessed local sequence context; whole-template span coverage is retained and independently guards function transfer. Site mapping does not establish a complete catalytic domain, cofactor occupancy or activity.'))
    for row in graph_nodes:
        evidence.append(dict(protein_id=row['protein_ids'],complex_anchor=row['anchor'],layer='L4',criterion='coordinate_module_assignment',
            status=row['status'],reference_id=row['reference_id'],source=row['source'],details=json.dumps(row),
            limitation='Reference-segment homology and ownership-compatible linkage only; alternative pairings, fragmentation and changed protein partitions are explicit. Localization, physical assembly and novel-organization function remain unvalidated.'))
    for row in cross_pairs:
        for side in ['a','b']:
            evidence.append(dict(protein_id=row['protein_id'],complex_anchor=row['anchor_'+side],layer='L4',criterion='cross_anchor_partner_hypotheses',
                status=row['pairing_status'],reference_id=row['structure_id_'+side]+':'+row['node_'+side],source=row['source_'+side],
                details=json.dumps(row),limitation=row['limitation']))
    for d in domain_rows:
        evidence.append(dict(protein_id=d["protein_id"],complex_anchor="",layer="L3",criterion="independent_profile_domain",
            status="supported" if d["passes"] else "below_domain_coverage",reference_id=d["accession"],source=d["source"],
            details=json.dumps(d),limitation="Shared domains may support multiple families; cofactor occupancy unknown"))
    for a in annotations:
        source=a.get('evidence_url') or str(ROOT/'references/processed/profile_manifest.json')
        for layer,criterion,status,detail,limitation in [
            ('L2','active_site_model_coverage',a['active_center_status'],
                dict(template=a['active_center_template'],positions=sum(len(r['sites']) for r in active_rows if r['protein_id']==a['protein_id'])),
                'Reviewed structural templates cover complete protein ligand sets and selected neighbors; local context quality and whole-template span quality are separate. No mapping does not prove ligand absence; independent calibration remains pending'),
            ('L3','FeS_mapping_coverage','mapped' if any(a['protein_id'] in s['complex_anchors'] for s in fes_rows) else 'not_evaluable',
                dict(site_count=sum(a['protein_id'] in s['complex_anchors'] for s in fes_rows)),
                'No reviewed match does not demonstrate cofactor absence; Hmd has no required classical FeS small subunit'),
            ('L6','in_vivo_measurement','not_measured',dict(net_direction=a['net_direction'],physiology=a['in_vivo_physiology']),
                'Sequence annotation cannot establish actual activity, cofactor occupancy or physiological flux')]:
            evidence.append(dict(protein_id=a['protein_id'],complex_anchor=a['protein_id'],layer=layer,criterion=criterion,status=status,
                reference_id=a['structure_id'],source=source,details=json.dumps(detail),limitation=limitation))
        for layer,field in [("L1","catalytic_status"),("L4","structure_status"),("L5","function_status")]:
            evidence.append(dict(protein_id=a["protein_id"],complex_anchor=a["protein_id"],layer=layer,criterion=field,status=a[field],
                reference_id=a["structure_id"],source=source,
                details=json.dumps({k:a.get(k,"") for k in ["hmm_support","missing_roles","flags","potential_function_ids"]}),
                limitation="Uncalibrated seed rules; in vivo activity and net direction unresolved"))
    write_tsv(out/"annotation_evidence.tsv",evidence,["protein_id","complex_anchor","layer","criterion","status","reference_id","source","details","limitation"])
    translation=translation_evidence(input_files, mode=mode, requested=translation_table,
        generated_gff=(out/'prodigal.gff') if (out/'prodigal.gff').exists() else None)
    manifest=dict(status="complete",input_files=[dict(path=str(Path(p).resolve()),sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest()) for p in input_files],
        protein_count=len(genes),complex_candidates=len(annotations),thresholds=thresholds,input_mode=mode,
        threads=threads,translation_table=translation_table,
        translation_table_requested=translation_table,
        translation_table_legacy_semantics='requested CLI/default value; not evidence of the effective input translation model',
        translation_table_used=translation['used'],
        translation_table_evidence=translation['table_evidence'],
        translation_exception_evidence=translation['exception_evidence'],
        software=software_versions(),
        output_sha256={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in [
            'complexes.tsv','modules.tsv','potential_functions.tsv','catalytic_anchors.tsv','catalytic_sites.tsv','complex_members.tsv',
            'reaction_capabilities.tsv','maturation.tsv','domain_hits.tsv','annotation_evidence.tsv','fes_sites.tsv','fes_ligands.tsv',
            'fes_template_evidence.tsv','family_candidates.tsv','sequence_evidence.tsv','rejected.tsv','genes.tsv',
            'active_center_summary.tsv','active_center_template_evidence.tsv','seed_diagnostic_sites.tsv','module_graph_nodes.tsv','module_graph_edges.tsv','cross_anchor_associations.tsv',
            'catalytic_center_candidates.tsv','catalytic_center_components.tsv','catalytic_confirmation.tsv','complex_hypotheses.tsv','reaction_contract_candidates.tsv',
            'structure_extensions/search_manifest.json','structure_extensions/queries.faa','structure_extensions/blast.tsv',
            'family.domtbl','family.tbl','blast.tsv']},
        bin_assignments={g["id"]:g.get("bin_id","") for g in genes} if mode=="bin" else {},
        reference_manifest=reference_manifest,
        reference_role_audit_scope=reference_role_audit['scope'],
        reference_fasta_sha256=hashlib.sha256((ROOT/"references/processed/references.faa").read_bytes()).hexdigest(),
        systems_config_sha256=hashlib.sha256((ROOT/"config/systems.json").read_bytes()).hexdigest(),
        annotation_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256((ROOT/"references/processed/catalytic_families.hmm").read_bytes()).hexdigest(),
        discovery_profile_manifest_sha256=hashlib.sha256(profile_manifest_path.read_bytes()).hexdigest(),
        discovery_evidence_policy_sha256=hashlib.sha256(discovery_policy_path.read_bytes()).hexdigest(),
        discovery_evidence_script_sha256=hashlib.sha256((ROOT/'scripts/discovery_evidence.py').read_bytes()).hexdigest(),
        discovery_stage_fingerprint=stage_fingerprint(ROOT),
        maturation_domain_model_sha256=domain_reference["combined_model_sha256"],
        maturation_domain_policy_sha256=domain_reference["policy_sha256"],
        fes_model_manifest_sha256=hashlib.sha256((ROOT/'references/processed/fes/manifest.json').read_bytes()).hexdigest(),
        fes_site_policy_sha256=fes_manifest['policy_sha256'],
        active_center_model_manifest_sha256=hashlib.sha256((ROOT/'references/processed/active_centers/manifest.json').read_bytes()).hexdigest(),
        active_center_policy_sha256=active_manifest['policy_sha256'],
        active_center_script_sha256=hashlib.sha256((ROOT/'scripts/active_centers.py').read_bytes()).hexdigest(),
        split_center_policy_sha256=split_manifest['policy_sha256'],
        split_center_model_manifest_sha256=hashlib.sha256((ROOT/'references/processed/split_centers/manifest.json').read_bytes()).hexdigest(),
        split_center_script_sha256=hashlib.sha256((ROOT/'scripts/split_centers.py').read_bytes()).hexdigest(),
        catalytic_confirmation=confirmation_manifest,
        complex_selection=selection_manifest,
        reaction_contract_candidates=contract_manifest,
        exact_reference_functions=exact_function_manifest,
        noncanonical_center_evidence=noncanonical_manifest,
        normalized_proteins_sha256=hashlib.sha256(faa.read_bytes()).hexdigest(),
        module_graph_policy_sha256=hashlib.sha256((ROOT/"config/module_graph.json").read_bytes()).hexdigest(),
        module_graph_script_sha256=hashlib.sha256((ROOT/"scripts/module_graph.py").read_bytes()).hexdigest(),
        references_json_sha256=hashlib.sha256((ROOT/"references/processed/references.json").read_bytes()).hexdigest(),
        limitations="Input mode controls partner linkage. Reference-system retrieval is not independent accuracy validation. Function never labelled experimentally proven for a query. Complete template protein ligand sets and selected neighbors are mapped independently of legacy diagnostics; broad confirmation calibration and full domain architecture remain pending. Per-site FeS structural homology is provisional, occupancy unknown.")
    origin = out / "circular_origin.json"
    if origin.is_file():
        manifest["output_sha256"][origin.name] = hashlib.sha256(origin.read_bytes()).hexdigest()
        manifest["circular_origin_evidence"] = origin.name
    (out/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n")
    return annotations


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--genbank",type=Path)
    group.add_argument("--proteins",type=Path)
    parser.add_argument("--genome",type=Path)
    parser.add_argument("--gff",type=Path)
    parser.add_argument("--mode",choices=["assembly","genome","bin","transcript"],default="assembly",
        help="Default assembly forbids cross-contig partners; genome asserts one trusted genome; bin requires --bin-map; transcript accepts translated proteins only")
    parser.add_argument("--bin-map",type=Path,help="TSV with contig and bin_id columns; required for bin mode")
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--threads",type=int,default=2)
    parser.add_argument("--translation-table",type=int,default=11)
    parser.add_argument("--topology",choices=["auto","linear","circular"],default="auto",
        help="Explicit topology for genome/GenBank input; auto reads GenBank annotations or FASTA header declarations")
    parser.add_argument("--phylogeny",action="store_true",help="Build review trees with MAFFT + IQ-TREE; slower, no automatic function upgrade")
    args=parser.parse_args()
    if not any([args.genome, args.genbank, args.proteins]):
        parser.error("one of --genome, --genbank, or --proteins is required")
    if args.genbank and args.genome:
        parser.error("--genbank cannot be combined with --genome")
    if args.threads<1:parser.error("--threads must be positive")
    if args.mode=="bin" and not args.bin_map:parser.error("--mode bin requires --bin-map")
    if args.bin_map and args.mode!="bin":parser.error("--bin-map requires --mode bin")
    if args.mode=="transcript" and (not args.proteins or args.gff):
        parser.error("transcript mode requires translated --proteins without GFF; eukaryotic genomic DNA is not supported by Prodigal")
    if args.genbank:
        genes,_=extract_genbank(args.genbank, None if args.topology=="auto" else args.topology)
        inputs=[args.genbank]
    elif args.proteins:
        if args.mode=="transcript":
            genes=[dict(id=r.id,contig=r.id,start=1,end=len(r.seq),strand=0,gene="",aliases="",protein_accession="",
                partial=False,contig_length=0,sequence=str(r.seq),coordinate_evidence="unavailable_translated_transcript")
                for r in SeqIO.parse(args.proteins,"fasta")]
            inputs=[args.proteins]
        else:
            if not args.gff:parser.error("--proteins requires a matching --gff")
            genes,_=read_faa_gff(args.proteins,args.gff,None if args.topology=="auto" else args.topology)
            inputs=[args.proteins,args.gff] + ([args.genome] if args.genome else [])
    else:
        args.out.mkdir(parents=True,exist_ok=True)
        # Prodigal uses an internal GFF ID differing from its FASTA ID; normalize explicitly.
        faa,gff=args.out/"prodigal.faa",args.out/"prodigal.gff"
        genes=[]
        tool_genome = materialize_fasta(args.genome, args.out)
        records=list(fasta_records(tool_genome))
        lengths={r.id:len(r) for r in records}
        if len(lengths) != len(records):
            raise ValueError("Duplicate nucleotide replicon IDs")
        topologies={r.id:(args.topology if args.topology!="auto" else
            ("circular" if (re.search(r"(?:^|[_-])circular-(?:yes|true|1)(?:[_-]|$)", r.id, re.I) or
                              re.search(r"(?:^|[; ])topology=circular(?:[; ]|$)", r.description, re.I)) else
             ("linear" if re.search(r"(?:^|[; ])topology=linear(?:[; ]|$)", r.description, re.I) else None)))
                     for r in records}
        if "circular" in topologies.values():
            recover_circular_calls(tool_genome, records, topologies, args.out, args.translation_table)
        else:
            call_prodigal(tool_genome, faa, gff, args.out/"prodigal.log", args.translation_table)
        for r in SeqIO.parse(faa,"fasta"):
            start,end,strand,attrs=protein_call(r)
            contig=r.id.rsplit("_",1)[0]
            genes.append(dict(id=r.id,contig=contig,start=start,end=end,strand=strand,
                gene="",aliases="",protein_accession="",partial=attrs.get("partial","00")!="00",
                contig_length=lengths[contig],sequence=str(r.seq).rstrip("*"),topology=topologies.get(contig)))
        inputs=[args.genome]
    if args.bin_map:
        read_bin_map(args.bin_map,genes)
        inputs.append(args.bin_map)
    result=annotate(genes,args.out,args.threads,inputs,args.mode,translation_table=args.translation_table if args.genome else None)
    if args.phylogeny:
        from phylogeny import build_trees
        build_trees(args.out,args.threads)
    print(f"{len(genes)} proteins, {len(result)} complex candidates: {args.out/'complexes.tsv'}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
