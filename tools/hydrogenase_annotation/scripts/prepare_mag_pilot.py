"""Freeze a metadata-stratified original-MAG pilot without inventing sequence access."""
import argparse
import collections
import hashlib
import json
import math
import re
from pathlib import Path
from datetime import datetime, timezone
import csv

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def profile(row, policy):
    ranked = dict(part.split('__', 1) for part in row['classification'].split(';') if '__' in part)
    completeness, contamination = float(row['completeness']), float(row['contamination'])
    n50, size = int(row['N50']), int(row['size'])
    if not (math.isfinite(completeness) and 0 <= completeness <= 100 and math.isfinite(contamination)
            and contamination >= 0 and size > 0 and 0 < n50 <= size):
        raise ValueError('Invalid quality/length metadata: ' + row['bin'])
    subtypes = sorted({label for label, count in re.findall(
        r'(\[(?:NiFe|FeFe|Fe)\](?: Group [\w.]+)?):(\d+)', row['target_summary']) if int(count) > 0})
    groups = sorted({('NiFe_Group4' if label.startswith('[NiFe] Group 4') else 'NiFe_other'
                      if label.startswith('[NiFe]') else 'FeFe' if label.startswith('[FeFe]') else 'Fe')
                     for label in subtypes})
    c1, c2 = policy['completeness_boundaries_percent']; n1, n2 = policy['n50_boundaries_bp']
    record = dict(bin=row['bin'], classification=row['classification'], domain=ranked.get('d') or 'unknown',
                  phylum=ranked.get('p') or 'unknown', completeness_percent=completeness,
                  contamination_percent=contamination, source_quality_category=row['quality_category'],
                  completeness_band=f'lt{c1}' if completeness < c1 else f'{c1}_to_lt{c2}' if completeness < c2 else f'ge{c2}',
                  n50_bp=n50, genome_size_bp=size,
                  n50_band=f'lt{n1}' if n50 < n1 else f'{n1}_to_lt{n2}' if n50 < n2 else f'ge{n2}',
                  old_subtypes=subtypes, potential_group='+'.join(groups) or 'no_old_hydrogenase_label',
                  historical_annotation_tsv=row['resolved_file'], sequence_input=None, input_status='sequence_not_resolved')
    record['joint_stratum'] = '|'.join(record[k] for k in ['domain', 'completeness_band', 'n50_band', 'potential_group'])
    record['selection_hash'] = hashlib.sha256((policy['seed'] + ':' + row['bin']).encode()).hexdigest()
    return record


def select(records, policy):
    ids = [r['bin'] for r in records]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate bin identifiers cannot define a sampling population')
    size = policy['sample_size']
    if not 0 < size <= len(records):
        raise ValueError('Sample size must be positive and no larger than population')
    ranked = sorted(records, key=lambda r: (r['selection_hash'], r['bin']))
    levels, strata = collections.defaultdict(list), collections.defaultdict(list)
    for r in ranked:
        strata[r['joint_stratum']].append(r)
        for dimension in policy['coverage_dimensions']:
            for value in (r[dimension] if dimension == 'old_subtypes' else [r[dimension]]):
                levels[(dimension, value)].append(r)
    selected, reasons = {}, collections.defaultdict(list)
    for (dimension, value), members in sorted(levels.items()):
        row = members[0]; selected[row['bin']] = row
        reasons[row['bin']].append(dimension + ':' + value)
    if len(selected) > size:
        raise ValueError('Sample too small for preregistered coverage; do not silently drop strata')
    counts = collections.Counter(r['joint_stratum'] for r in selected.values())
    offsets = collections.Counter()
    while len(selected) < size:
        choices = [s for s in sorted(strata) if counts[s] < len(strata[s])]
        s = max(choices, key=lambda s: math.sqrt(len(strata[s])) / (counts[s] + 0.5))
        while strata[s][offsets[s]]['bin'] in selected:
            offsets[s] += 1
        row = strata[s][offsets[s]]; offsets[s] += 1
        selected[row['bin']] = row; counts[s] += 1
        reasons[row['bin']].append('sqrt_stratum_allocation')
    return [dict(selected[key], selection_reasons=reasons[key]) for key in sorted(selected)]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    policy_path = ROOT / 'config/mag_pilot_sampling.json'
    policy = json.loads(policy_path.read_text()); source = (ROOT / policy['population_file']).resolve()
    with source.open() as handle:
        records = [profile(r, policy) for r in csv.DictReader(handle, delimiter='\t')]
    if len(records) != policy['population_rows']:
        raise ValueError('Original population size changed; review protocol before sampling')
    chosen = select(records, policy)
    coverage = {}
    for dimension in policy['coverage_dimensions']:
        def count(rows):
            return dict(sorted(collections.Counter(v for r in rows for v in
                (r[dimension] if dimension == 'old_subtypes' else [r[dimension]])).items()))
        population, sample = count(records), count(chosen)
        if set(population) != set(sample):
            raise ValueError('Observed metadata level omitted: ' + dimension)
        coverage[dimension] = dict(population=population, selected=sample)
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    (out / 'selected_samples.json').write_text(json.dumps(chosen, indent=2) + '\n')
    freeze = dict(status='metadata_sample_frozen_sequences_unresolved_not_run', created_utc=datetime.now(timezone.utc).isoformat(),
                  source_file=str(source), source_sha256=sha(source), policy_sha256=sha(policy_path),
                  script_sha256=sha(Path(__file__)), population=len(records), selected=len(chosen),
                  artifacts={'selected_samples.json': sha(out / 'selected_samples.json')},
                  coverage=coverage, annotation_runs=0, accessible_sequence_inputs=0,
                  limitations=policy['rules'])
    (out / 'freeze.json').write_text(json.dumps(freeze, indent=2) + '\n')
    print(json.dumps({k: freeze[k] for k in ['status', 'population', 'selected', 'annotation_runs']}))


if __name__ == '__main__':
    main()
