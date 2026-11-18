#!/usr/bin/env python3
"""Export runtime code, models, and their required validation evidence."""
import argparse, hashlib, json, re, shutil
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''): h.update(block)
    return h.hexdigest()


PATH_PREFIXES = ('config/', 'experiments/', 'references/', 'research_workflow/',
                 'scripts/', 'tools/', 'results/', 'revisions/', 'reports/')


def existing_path(root, value):
    """Recognize repository paths before inspecting large sequence/provenance JSON."""
    if not isinstance(value, str) or not value.startswith(PATH_PREFIXES):
        return None
    if len(value) > 300 or any(c.isspace() for c in value):
        return None
    p = (root / value).resolve()
    try:
        rel = p.relative_to(root.resolve())
    except ValueError:
        return None
    if any(x in rel.parts for x in ('.env', '__pycache__')) or p.suffix in ('.pyc', '.pyo'):
        return None
    return rel.as_posix() if p.is_file() or p.is_dir() else None


def existing_relative(root, value):
    rel = existing_path(root, value)
    return rel if rel and (root / rel).is_file() else None


def walk_strings(value):
    if isinstance(value, dict):
        for k, v in value.items():
            yield None, k
            yield k, v
            yield from walk_strings(v)
    elif isinstance(value, list):
        for v in value:
            if isinstance(v, str): yield None, v
            else: yield from walk_strings(v)


def collect(root):
    files = {}
    path_cache = {}
    def lookup(value):
        if not isinstance(value, str) or not value.startswith(PATH_PREFIXES):
            return None
        if value not in path_cache:
            path_cache[value] = existing_path(root, value)
        return path_cache[value]
    # Validation locks sometimes name an annotation directory while the
    # validator consumes only this fixed artifact contract.  Do not turn that
    # shorthand into a multi-gigabyte experiment archive.
    directory_contract = {
        'proteins.faa', 'confirmation/features.json',
        'catalytic_confirmation.tsv', 'active_center_summary.tsv',
        'catalytic_sites.tsv', 'catalytic_anchors.tsv',
    }
    def add(rel, role, source):
        rel = Path(rel).as_posix(); p = root / rel
        if p.is_file() and '.env' not in Path(rel).parts and '__pycache__' not in Path(rel).parts and not p.name.endswith(('.pyc', '.pyo')):
            if rel not in files:
                files[rel] = dict(path=rel, sha256=sha(p), roles=[], sources=[])
            if role not in files[rel]['roles']: files[rel]['roles'].append(role)
            if source not in files[rel]['sources']: files[rel]['sources'].append(source)
    for directory, role in [('scripts','stage_script'),('research_workflow','workflow_code'),('config','config'),('references/processed','processed_reference')]:
        for p in (root / directory).rglob('*'):
            if p.is_file(): add(p.relative_to(root), role, directory)
    for rel in ('run_hydrogenase.py', 'run_batch.py', 'scripts/install_prodigal.py',
                'tools/prodigal/manifest.json', 'docs/batch.md', 'docs/install_frozen.md'):
        add(rel, 'runtime_entrypoint' if rel.endswith(('.py', '.sh')) else 'pinned_tool_source', 'explicit_runtime_closure')
    for p in (root / 'tools/prodigal/source').rglob('*'):
        if p.is_file(): add(p.relative_to(root), 'pinned_tool_source', 'tools/prodigal/manifest.json')
    for p in (root / 'tools/pfsearch_runtime').rglob('*'):
        if p.is_file(): add(p.relative_to(root), 'pinned_native_runtime', 'tools/pfsearch_runtime/manifest.json')
    domain_policy = root / 'config/reference_domain_policy.json'
    if domain_policy.is_file():
        for model in json.loads(domain_policy.read_text())['models']:
            rel = 'references/raw/reference_domain_review/' + model['accession'] + '.hmm'
            if not (root / rel).is_file():
                raise FileNotFoundError(rel)
            add(rel, 'curated_domain_profile', str(domain_policy.relative_to(root)))

    # Follow source manifests and active validation locks transitively.  This is
    # Hash-map keys are paths too; absolute provenance strings are not package
    # lookups, while relative files/directories inside the repository are.
    pending = []
    seen_json = set()
    def expandable_manifest(rel):
        rel = str(rel)
        return (rel.startswith('config/') or rel.startswith('references/validation_locks/') or
                rel.startswith('references/processed/') or rel.startswith('references/validation/') or
                rel.startswith('research_workflow/') or
                rel.endswith('/topology_transition_v1.json'))
    def queue(rel, role, source):
        rel = Path(rel).as_posix()
        source = str(source)
        p = root / rel
        if p.is_dir():
            for child in p.rglob('*'):
                child_rel = child.relative_to(p).as_posix()
                if child.is_file() and child_rel in directory_contract:
                    add(child.relative_to(root), role, source)
        elif p.is_file():
            add(rel, role, source)
            if p.suffix.lower() == '.json' and rel not in seen_json:
                pending.append((rel, role)); seen_json.add(rel)
    for p in sorted((root/'research_workflow').rglob('*.json')):
        queue(p.relative_to(root), 'workflow_model_or_review', str(p.relative_to(root)))
    for p in sorted((root/'references/processed').rglob('*.json')):
        queue(p.relative_to(root), 'processed_manifest', str(p.relative_to(root)))
    # Runtime policies (not only *_model.json files) can name a validation
    # corpus directory whose freeze manifest is loaded by ordinary stages.
    for p in sorted((root/'config').glob('*.json')):
        queue(p.relative_to(root), 'runtime_policy', str(p.relative_to(root)))
        try: data=json.loads(p.read_text())
        except (OSError, json.JSONDecodeError): continue
        for key, value in walk_strings(data):
            if isinstance(value, str) and value.startswith('references/validation/'):
                target=root/value
                if target.is_dir():
                    for child in target.rglob('*'):
                        if child.is_file(): queue(child.relative_to(root), 'declared_config_source', str(p.relative_to(root)))
    for p in sorted((root/'references/validation_locks').glob('*.json')):
        try: data=json.loads(p.read_text())
        except (OSError, json.JSONDecodeError): continue
        queue(p.relative_to(root), 'validation_lock', p.name)
        for key, value in walk_strings(data):
            rel=lookup(value)
            if rel: queue(rel, 'lock_declared_source', p.name)

    # A few ordinary stages name reference files in code rather than a model.
    # Capture only existing, repository-relative literals, never globs or dirs.
    literal = re.compile(r"(?:config|experiments|references|research_workflow|scripts|tools|reports|revisions)/(?:[A-Za-z0-9_.-]+/)*[A-Za-z0-9_.-]+")
    for base in (root/'scripts', root/'research_workflow'):
        for p in base.rglob('*.py'):
            if '__pycache__' in p.parts: continue
            try: text = p.read_text()
            except OSError: continue
            for value in literal.findall(text):
                # Code literals are normally files.  A small validation corpus
                # is intentionally named as a directory and is safe to close
                # over; never do this for experiment directories.
                rel = existing_relative(root, value)
                if rel:
                    queue(rel, 'code_literal_source', str(p.relative_to(root)))
                elif value.startswith('references/validation/'):
                    rel = lookup(value)
                    if rel:
                        target = root / rel
                        if target.is_dir():
                            for child in target.rglob('*'):
                                if child.is_file(): queue(child.relative_to(root), 'code_literal_source', str(p.relative_to(root)))
                        else:
                            queue(rel, 'code_literal_source', str(p.relative_to(root)))

    # Resolve policies and locks after all entry points have been collected.
    while pending:
        rel, _ = pending.pop(0)
        p = root / rel
        if not expandable_manifest(rel): continue
        try: data = json.loads(p.read_text())
        except (OSError, UnicodeDecodeError, json.JSONDecodeError): continue
        for key, value in walk_strings(data):
            child = lookup(value)
            if child: queue(child, 'recursive_declared_source', rel)
    return [files[k] for k in sorted(files)]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--map',type=Path,required=True)
    ap.add_argument('--export',type=Path)
    args=ap.parse_args(); root=args.root.resolve(); rows=collect(root)
    args.map.parent.mkdir(parents=True,exist_ok=True)
    args.map.write_text(json.dumps(dict(status='candidate_closure',root=str(root),file_count=len(rows),files=rows),ensure_ascii=False,indent=2)+'\n')
    if args.export:
        dest=args.export.resolve()
        if dest.exists(): raise FileExistsError(dest)
        for row in rows:
            target=dest/row['path']; target.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(root/row['path'],target)
            if sha(target)!=row['sha256']: raise ValueError('copy hash mismatch: '+row['path'])
        (dest/'package_manifest.json').write_text(json.dumps(dict(status='exported',file_count=len(rows),files=rows),ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status='exported' if args.export else 'mapped',file_count=len(rows),map=str(args.map),export=str(args.export) if args.export else None)))


if __name__=='__main__': main()
