"""Apply reviewed identity corrections without rewriting historical source labels."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def adjudications(root=ROOT):
    root=Path(root);path=root/'config/catalytic_identity_adjudications.json'
    review=json.loads(path.read_text());result={}
    for row in review['records']:
        digest=row['sequence_sha256']
        if row['identity_label'] not in ['positive','negative','indeterminate'] or digest in result:
            raise ValueError('Identity adjudication correspondence changed: '+row['unit_id'])
        result[digest]=row
    return result


def load_units(corpus,root=ROOT):
    rows=json.loads((Path(corpus)/'units.json').read_text());review=adjudications(root)
    for row in rows:
        correction=review.get(row['sequence_sha256'])
        if correction:
            row['frozen_label']=row['label'];row['identity_adjudication']=correction
            row['label']=correction['identity_label']
            row['categories']=sorted(set(row['categories'])|{correction['family']})
    return rows
