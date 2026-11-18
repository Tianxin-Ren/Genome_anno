# Batch entry point

`run_batch.py` consumes a JSON list of `{sample, path, sha256}` records. The existing `run_hydrogenase.py` is invoked once per sample with the full seven-stage workflow. The compatible historical Hadar manifest is `results/hadar_pilot_round49_Fes_expansion/input_manifest.json`.

Example:

```bash
.env/bin/python run_batch.py manifest.json --out results/batch --mode assembly --threads 2
.env/bin/python run_batch.py manifest.json --out results/batch --resume --mode assembly --threads 2
```

Each sample is isolated under `out/<sample>/attempt_001/`, then `attempt_002/` on a retry; the runner log is beside the attempt directory. Sample names are restricted to one safe path component and SHA-256 is mandatory. Failures do not stop later samples. `batch_state.json`, `summary.tsv`, and `summary.json` are written atomically after each sample. A resumed batch validates the manifest path/content, every input SHA-256, arguments, entrypoint hash, and the runtime source/lock fingerprint. A sample is skipped only when its workflow is complete, is a fresh non-reused seven-stage run, its annotation input hash matches the manifest, and every recorded output hash still verifies; status alone is never sufficient.

SIGINT terminates the current child process group, records `interrupted`, and exits 130. POSIX `flock` prevents two writers from using one output directory. The runner intentionally has no stage-level resume: unfinished samples get a fresh attempt and rerun from the input.

For a completed sample, `claim_assertions` links to the report, `claim_status` is
`valid`, and `claim_count` gives the number of assertions. Completion and resume
verify the report's hash, schema version, candidate IDs and count, and all six
dimension statuses. These summary fields are empty for pending or failed samples.

Each completed sample also includes `maturation_observation.json`. It records
the maturation search over the supplied proteins, distinguishes homology-only
hits from profile-supported architectures, and reports genomic observability
separately. No hit means no detection within that search scope; it does not
establish biological absence or remove a supported catalytic core. Maturation
factors may be remote or shared. The report is required and hash-checked before
a batch sample is considered complete or skipped on resume.

`translation_hypotheses.json` also records bounded split-component coding
hypotheses. A seed with incomplete local coverage may be checked against a
complete start-to-stop ORF using a global alignment over the same reviewed
structural window. The global mapping must agree with the seed's exact TGA/C
codons; ambiguous optimal mappings are rejected. This additional method is
uncalibrated and is marked on each affected coding alternative. Its `X` residues
are conditional TGA placeholders, carry no identity-match credit, and do not
establish Sec incorporation, a called protein, complex membership or function.
