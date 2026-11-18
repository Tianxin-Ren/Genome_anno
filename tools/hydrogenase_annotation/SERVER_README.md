# Hydrogenase Annotation Server Deployment & Parallel Runner Guide

## 1. Server Environment Setup (Linux / macOS)

Using `mamba` or `conda` with the pinned specification in `config/environment.yml`:

```bash
cd hydrogenase_annotation
mamba env create --prefix ./.env --file config/environment.yml
# Or with conda:
# conda env create --prefix ./.env --file config/environment.yml
conda activate ./.env
```

Verify core binaries (`prodigal`, `blastp`, `makeblastdb`, `tblastn`, `hmmsearch`, `mafft`, `iqtree3`):

```bash
./.env/bin/python -c "import Bio, openpyxl; print('Python deps OK')"
prodigal -v && blastp -version && hmmsearch -h | head -n 2
```

*(Optional)* If running PROSITE `PS51384` FAD_FR profile scans (`research_workflow/prosite_fad.py`), install `pftools` (`conda install -c bioconda pftools`) and ensure `pfsearch2.3` (or `pfsearch`) is in `./.env/bin/`.

---

## 2. High-Concurrency Parallel Execution (e.g., 40 Cores)

Use `run_parallel.py` to annotate a directory of genomes (`.fna`, `.fa`, `.fasta`, or `.gz` compressed) concurrently and automatically aggregate results across the cohort:

```bash
./.env/bin/python -u run_parallel.py \
  --input-dir /path/to/genomes \
  --out /path/to/results_batch1 \
  --workers 40 \
  --threads-per-worker 1 \
  --mode genome \
  --resume
```

Key options:
- `--workers 40 --threads-per-worker 1`: runs 40 genomes simultaneously, pinning `OMP_NUM_THREADS=1` / `OPENBLAS_NUM_THREADS=1` per worker to avoid thread oversubscription.
- `--manifest <manifest.tsv|manifest.json>`: alternative to `--input-dir` when using explicit `sample` and `path` columns.
- `--resume`: skips genomes that already completed (`workflow.json` status `complete`) and archives/retries any incomplete or failed runs.
- `--aggregate-only`: re-generates cohort summary tables from an existing `--out` directory without re-running predictions.

### Cohort Outputs in `--out`

- `samples/<sample>/`: full per-genome 7-stage output (`workflow.json`, `candidate_function_catalog.tsv`, `functions/functional_annotations.tsv`, `fes_refinement/fes_sites.tsv`, etc.).
- `logs/<sample>.log`: per-genome stdout/stderr log.
- `cohort_summary.tsv`: 1 row per genome (status, elapsed time, candidate count, core anchors, predicted functional classes, unresolved core count).
- `cohort_functional_annotations.tsv`: merged table of all non-trivial / functionally assigned hydrogenase candidates across all genomes.
- `cohort_class_matrix.tsv`: sample × functional-class count matrix (columns are dynamically read from `config/function_rules_34.json`).
- `cohort_manifest.json`: run summary, rule versions, and counts.

---

## 3. Extending or Modifying Functional Rules

The functional classification layer is decoupled from sequence discovery:

1. **Edit or rebuild rules**:
   - Rules live in `config/function_rules_34.json` and the member sequence bank in `references/processed/function_rules/members.faa`.
   - You can add new functional classes or update module/exclusivity rules either by editing `scripts/build_function_rules.py` and re-running it, or by updating `config/function_rules_34.json` + `references/processed/function_rules/members.faa` (updating `fasta_sha256` and per-reference `sha256`).
2. **Validate before deploying**:
   ```bash
   ./.env/bin/python -m unittest discover -s tests -v
   ./.env/bin/python -m unittest discover -s research_workflow -p "test_*.py" -v
   ./.env/bin/python scripts/validate_modular_function_rules.py
   ```
3. **Re-run or re-aggregate**:
   - `cohort_class_matrix.tsv` in `run_parallel.py` automatically picks up any new class IDs from `config/function_rules_34.json`.
