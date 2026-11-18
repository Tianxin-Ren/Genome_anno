"""Build occupancy-trimmed branch HMMs and retain every training sequence ID."""
import hashlib
import subprocess
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def profile_training(root,config,records):
    """Only explicitly registered development identities may expand discovery training."""
    pool=dict(records);provenance={}
    if config.get('development_seed_manifest'):
        spec=config['development_seed_manifest'];path=(root/spec['path']).resolve()
        manifest=json.loads(path.read_text())
        fasta=(root/manifest['fasta_file']).resolve()
        with fasta.open() as handle:extra=list(SeqIO.parse(handle,'fasta'))
        for seq in extra:
            if seq.id in pool:raise ValueError('Discovery training ID overlaps an existing reference')
            pool[seq.id]=seq
        provenance=dict(development_seed_manifest=spec,training_corpus_unit_ids=[r['unit_id'] for r in manifest['records']],training_truth_scope=manifest.get('truth_scope'))
    if any(rid not in pool for rid in config['reference_ids']):
        raise ValueError('Unknown branch training reference')
    return [pool[rid] for rid in config['reference_ids']],provenance


def core_columns(msa,occupancy,reference_core=None):
    """Mask from training occupancy and an optional sequence-pinned structural interval."""
    if not msa or not 0 < occupancy <= 1 or len({len(s) for s in msa}) != 1:
        raise ValueError('Invalid training alignment or occupancy')
    keep=[i for i in range(len(msa[0])) if sum(str(s.seq)[i]!="-" for s in msa)/len(msa)>=occupancy]
    if reference_core:
        anchors=[s for s in msa if s.id==reference_core['reference_id']]
        if len(anchors)!=1:
            raise ValueError('Expected exactly one reviewed core anchor in training')
        seq=str(anchors[0].seq).replace('-','')
        if hashlib.sha256(seq.encode()).hexdigest()!=reference_core['sequence_sha256']:
            raise ValueError('Reviewed core anchor sequence changed')
        start,end=reference_core['start'],reference_core['end']
        if not 1 <= start <= end <= len(seq):
            raise ValueError('Reviewed core interval outside anchor')
        if any(not start<=int(pos)<=end or seq[int(pos)-1]!=aa for pos,aa in reference_core['required_positions'].items()):
            raise ValueError('Reviewed core ligand coordinates or identities changed')
        positions=[i for i,a in enumerate(str(anchors[0].seq)) if a!='-']
        keep=[i for i in keep if positions[start-1]<=i<=positions[end-1]]
    return keep


def build_profile(records,directory,name,occupancy=0.6,reference_core=None):
    directory.mkdir(parents=True,exist_ok=True)
    raw,aligned,core,hmm=[directory/f"{name}.{suffix}" for suffix in ["faa","aligned.faa","core.faa","hmm"]]
    SeqIO.write([SeqRecord(Seq(str(s.seq).replace("U","X")),id=s.id,description="") for s in records],raw,"fasta")
    with aligned.open("w") as handle:
        subprocess.run(["mafft","--quiet","--auto",str(raw)],stdout=handle,check=True)
    msa=list(SeqIO.parse(aligned,"fasta"))
    keep=core_columns(msa,occupancy,reference_core)
    if len(keep)<100:
        raise ValueError(f"Insufficient catalytic core columns: {name}: {len(keep)}")
    SeqIO.write([SeqRecord(Seq("".join(str(s.seq)[i] for i in keep)),id=s.id,description="") for s in msa],core,"fasta")
    with (directory/f"{name}.hmmbuild.log").open("w") as handle:
        subprocess.run(["hmmbuild","--amino","-n",name,str(hmm),str(core)],stdout=handle,check=True)
    return hmm,dict(profile=name,sequences=[s.id for s in records],alignment_columns=len(msa[0]),retained_columns=len(keep),occupancy=occupancy,
        reference_core=reference_core,retained_alignment_columns_1based=[i+1 for i in keep],
        model_sha256=hashlib.sha256(hmm.read_bytes()).hexdigest(),
        training_sequence_sha256={s.id:hashlib.sha256(str(s.seq).encode()).hexdigest() for s in records})
