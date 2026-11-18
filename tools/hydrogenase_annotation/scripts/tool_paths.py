"""Resolve external bioinformatics binaries with local .env preference and PATH fallback."""
from pathlib import Path
import shutil


def resolve_tool(root: Path, name: str) -> Path:
    local = Path(root) / ".env/bin" / name
    if local.is_file():
        return local
    found = shutil.which(name)
    if found:
        return Path(found).resolve()
    return local
