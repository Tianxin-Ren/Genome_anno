# Standalone Protein HMM Annotation Workflow (config_hmm)

This configuration directory is customized to run **parallelized HMM annotation on standalone protein files (FAA)** without running gene prediction (Bakta) or circular assembly analysis.

---

## 📂 Configuration Files

1. **`pipeline.yaml`**: The master module registry. Only the `annotation` module is enabled.
2. **`annotation.yaml`**: Configures HMM options. Sample-specific HMM/database paths are defined in `samples.yaml`.
3. **`samples.yaml`**: Defines input protein files using the `type: "protein"` parameter.

---

## 🚀 How to Run

1. **Define your inputs in `samples.yaml`**:
   Specify a unique name for your run (e.g. `my_proteins`) and provide the path to your protein file:
   ```yaml
   external_genomes:
     my_proteins:
       path: "/absolute/path/to/your_proteins.faa"
       type: "protein"
   ```

2. **Pair each input with its HMM database in `samples.yaml`**:
   Set `path` and `hmm_db` in the same sample entry. `hmm_db` may be a folder containing `.hmm` files, a single concatenated `.hmm` file, or a list:
   ```yaml
   external_genomes:
     my_genomes:
       path: "/path/to/genomes"
       type: "dna_dir"
       hmm_db:
         - "/path/to/target.hmm"
       ko_list: "/path/to/ko_list"
   ```
   Platform database paths remain fallback defaults for samples that do not specify their own paths.

3. **Launch Snakemake**:
   Run the annotation and its downstream integration together on `slurm_pi`:
   ```bash
   snakemake --workflow-profile profiles/slurm_pi \
     --config conf=projects/config_hmm platform=config/platforms/slurm_pi.yaml
   ```

---

## 📊 Outputs

The results will be written to `results/annotation/{sample_name}/`:
* `advanced_hmm/unified_all_raw.tbl`: The gathered raw HMMER tblout file.
* `{sample_name}_advanced_hmm_resolved.tsv`: The resolved annotations with best-hit scoring (if `resolve_hmm` is set to `True`).

For batch tasks, `postprocess.yaml` additionally connects the batch functional
summary to the historical data-integration scripts. The current
`bathyarchaeia_gtdb_all` task produces:

* `results/postprocess/bathyarchaeia_gtdb_all/gene_counts.tsv`
* `results/postprocess/bathyarchaeia_gtdb_all/oxphos_evaluation.tsv`

## Large-scale protein redundancy strategy

When `annotation.protein_dedup.enabled` is enabled, the batch workflow first
creates a lossless SHA256 protein catalog. HMM target searches run against all
100%-unique proteins, and the result is expanded through `protein_map.tsv`
before per-genome resolution. This preserves rare target-HMM hits and avoids
running HMMER repeatedly on identical proteins.

In parallel, MMseqs2 creates a conservative 99% identity / 95% coverage cluster
under `results/external_annotation/<batch>_protein_dedup/cluster_99/`. Those
representatives are intended for future general-purpose annotation; they are
not used to suppress direct target-HMM searches.

The deduplication directory also contains `catalog.sqlite`, which is the
disk-backed index used for memory-bounded and audit-friendly reprocessing.
