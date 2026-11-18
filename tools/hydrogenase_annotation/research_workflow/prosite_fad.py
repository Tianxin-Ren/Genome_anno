"""Strict local scan of the official PROSITE PS51384 FAD_FR profile."""
import csv
import hashlib
import json
import math
import re
import subprocess
import string
from pathlib import Path

from Bio import SeqIO

AA = set(string.ascii_uppercase)
HIT = re.compile(
    r"^>(?P<sid>[^/\s]+)/(?P<start>\d+)-(?P<end>\d+) "
    r"motif=(?P<motif>[^ ]+) norm_score=(?P<norm>[-+0-9.eE]+) "
    r"raw_score=(?P<raw>[-+0-9.eE]+) level=(?P<level>-?\d+) "
    r"seq_end=(?P<seq_end>-?\d+) motif_start=(?P<mstart>-?\d+) motif_end=(?P<mend>-?\d+)\s*$"
)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _folder_sha(path):
    digest = hashlib.sha256()
    for item in sorted(Path(path).glob("pfsearch*")):
        if item.is_file():
            digest.update(item.name.encode()); digest.update(item.read_bytes())
    return digest.hexdigest()


def _profile_bounds(path):
    text = Path(path).read_text()
    m = re.search(r"/CUT_OFF: LEVEL=0; SCORE=(\d+);.*?N_SCORE=([0-9.]+)", text)
    if not m or int(m.group(1)) != 1658 or float(m.group(2)) != 8.5:
        raise ValueError("PS51384 profile has no official LEVEL-0 cutoff 1658/8.5")
    lm = re.search(r"LENGTH=(\d+)", text)
    if "ID   FAD_FR; MATRIX." not in text or "AC   PS51384;" not in text or not lm:
        raise ValueError("profile is not the official PS51384 FAD_FR record")
    return int(m.group(1)), float(m.group(2)), int(lm.group(1))


def _read_fasta(path):
    with Path(path).open() as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if not records:
        raise ValueError("annotation proteins.faa is empty")
    seen = set(); out = {}
    for record in records:
        if record.id in seen:
            raise ValueError("duplicate FASTA protein id: " + record.id)
        seen.add(record.id)
        seq = str(record.seq).upper()
        if not seq or set(seq) - AA:
            raise ValueError("invalid or empty protein sequence: " + record.id)
        out[record.id] = seq
    return out


def _parse_output(text, id_map, sequences, profile_sha, scanner_sha, cutoff):
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        m = HIT.match(line)
        if not m or line[m.end():].strip() and not line[m.end():].strip().startswith("seq_end="):
            raise ValueError("malformed pfsearch output: " + line)
        sid = m["sid"]
        if sid not in id_map:
            raise ValueError("unknown scanner FASTA id: " + sid)
        start, end = int(m["start"]), int(m["end"])
        raw, norm, level = float(m["raw"]), float(m["norm"]), int(m["level"])
        qid, seq = id_map[sid], sequences[id_map[sid]]
        seq_end, mstart, mend = int(m["seq_end"]), int(m["mstart"]), int(m["mend"])
        pstart = mstart if mstart > 0 else cutoff[2] + mstart + 1
        pend = mend if mend > 0 else cutoff[2] + mend + 1
        resolved_seq_end = seq_end if seq_end > 0 else len(seq) + seq_end + 1
        if not all(math.isfinite(v) for v in (raw, norm)) or start < 1 or end < start or end > len(seq) or resolved_seq_end != end or not (1 <= pstart <= pend <= cutoff[2]):
            raise ValueError("invalid PS51384 hit coordinates or scores: " + line)
        if m["motif"] != "PS51384|FAD_FR" or level != 0:
            raise ValueError("unexpected PS51384 motif/level: " + line)
        if norm < cutoff[1] or raw < cutoff[0]:
            raise ValueError("LEVEL-0 scanner output is below official cutoff: " + line)
        rows.append(dict(protein_id=qid, accession="PS51384", query_start=start, query_end=end,
            passes=True, raw_score=raw, norm_score=norm, level=level,
            source_sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
            profile_sha256=profile_sha, scanner_sha256=scanner_sha, query_length=len(seq),
            profile_start=pstart, profile_end=pend))
    return rows


def scan(annotation: Path, out: Path, model, root: Path):
    """Run one strict PS51384 scan and return rows plus all bound input hashes."""
    root, annotation, out = Path(root), Path(annotation), Path(out)
    cfg = model.get("prosite_fad")
    if not cfg or cfg.get("normalization_mode") != 1 or cfg.get("min_normalized_score") != 8.5 or cfg.get("min_raw_score") != 1658:
        raise ValueError("model lacks the locked PS51384 LEVEL-0 rule")
    import shutil
    profile = root / cfg["profile"]
    scanner = root / cfg["scanner"]
    if not scanner.exists():
        found = shutil.which(Path(cfg["scanner"]).name) or shutil.which("pfsearch")
        if found:
            scanner = Path(found)
    for path in (profile, scanner, annotation / "proteins.faa"):
        if not path.exists():
            hint = "; build the local runtime per experiments/round65/prosite_runtime/README.md" if path == scanner else ""
            raise FileNotFoundError("required PROSITE input missing: " + str(path) + hint)
    cutoff = _profile_bounds(profile)
    bound = {str(p): _sha(p) for p in (annotation / 'proteins.faa', profile, scanner)}
    sequences = _read_fasta(annotation / "proteins.faa")
    out.mkdir(parents=True, exist_ok=False)
    id_map = {f"p{i:08d}": qid for i, qid in enumerate(sequences, 1)}
    query = out / "queries.faa"
    # The reference Fortran reader has a fixed line buffer; wrap long proteins.
    query.write_text(''.join('>' + sid + '\n' + ''.join(seq[i:i+80] + '\n' for i in range(0, len(seq), 80))
                            for sid, qid in id_map.items() for seq in [sequences[qid]]))
    map_path = out / "faa_id_map.tsv"
    with map_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["scanner_id", "protein_id", "query_length", "sequence_sha256"])
        writer.writerows((sid, qid, len(sequences[qid]), hashlib.sha256(sequences[qid].encode()).hexdigest()) for sid, qid in id_map.items())
    command = [str(scanner), "-f", "-M1", "-C8.5", "-k", "-l", "-z", str(profile), str(query)]
    result = subprocess.run(command, cwd=str(root), text=True, capture_output=True, check=False)
    (out / "stdout.txt").write_text(result.stdout); (out / "stderr.txt").write_text(result.stderr)
    (out / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    if result.returncode:
        raise RuntimeError(f"pfsearch failed ({result.returncode}); see {out / 'stderr.txt'}")
    profile_sha, scanner_sha = _sha(profile), _sha(scanner)
    rows = _parse_output(result.stdout, id_map, sequences, profile_sha, scanner_sha, cutoff)
    row_path = out / "rows.tsv"
    fields = list(rows[0]) if rows else ["protein_id", "accession", "query_start", "query_end", "passes", "raw_score", "norm_score", "level", "source_sequence_sha256", "profile_sha256", "scanner_sha256", "query_length", "profile_start", "profile_end"]
    with row_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    if any(_sha(path) != digest for path, digest in bound.items()):
        raise ValueError('PROSITE source changed during scan')
    artifacts = {name: _sha(out / name) for name in ("queries.faa", "faa_id_map.tsv", "stdout.txt", "stderr.txt", "command.json", "rows.tsv")}
    manifest = {"profile": "PS51384", "cutoff": {"raw": cutoff[0], "normalized": cutoff[1], "level": 0},
        "status": "complete", "sequence_count": len(sequences), "hit_count": len(rows), "command": command,
        "returncode": result.returncode, "query_sha256": _sha(query), "artifacts_sha256": artifacts,
        "scanner_family_sha256": _folder_sha(scanner.parent), "bound_inputs": bound}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return rows, bound
