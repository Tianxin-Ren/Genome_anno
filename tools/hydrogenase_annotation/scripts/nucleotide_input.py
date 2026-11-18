"""Small shared reader for plain and gzip-compressed nucleotide FASTA files."""
from contextlib import contextmanager
import gzip
import shutil
from pathlib import Path

from Bio import SeqIO


def is_gzip(path):
    return Path(path).name.lower().endswith(".gz")


@contextmanager
def open_text(path):
    handle = gzip.open(path, "rt") if is_gzip(path) else Path(path).open()
    try:
        yield handle
    finally:
        handle.close()


def fasta_records(path):
    with open_text(path) as handle:
        yield from SeqIO.parse(handle, "fasta")


def materialize_fasta(path, directory):
    """Give command-line tools a plain FASTA while preserving the source path/hash."""
    source = Path(path)
    if not is_gzip(source):
        return source
    output = Path(directory) / ".input_genome.fna"
    with gzip.open(source, "rb") as handle, output.open("wb") as target:
        shutil.copyfileobj(handle, target)
    return output
