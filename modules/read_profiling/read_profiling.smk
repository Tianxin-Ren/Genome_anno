# ==============================================================================
# Metagenomic Read Profiling Module (Sylph Taxonomic Profiling)
# ==============================================================================
import os
import glob
import re
import sys
import csv
import tempfile

# 1. Global Path Logic (Namespaced)
RP_C = config.get("read_profiling", {}) or config.get("sylph", {})
RP_ROOT = config.get("project_root", "") or config.get("paths", {}).get("project_root", "")
if RP_ROOT:
    RP_ROOT = os.path.expanduser(os.path.expandvars(str(RP_ROOT)))
    if not os.path.isabs(RP_ROOT):
        RP_ROOT = os.path.abspath(RP_ROOT)
else:
    RP_ROOT = os.getcwd()


def rp_full_path(path):
    if not path:
        return ""
    path_str = os.path.expanduser(os.path.expandvars(str(path)))
    if os.path.isabs(path_str):
        return os.path.normpath(path_str)
    if RP_ROOT:
        return os.path.normpath(os.path.join(RP_ROOT, path_str))
    return os.path.normpath(path_str)


def rp_resolve_path(path):
    """Use an existing path, then try configured and cluster storage aliases."""
    if not path:
        return ""
    p = rp_full_path(path)

    # This sample sheet mixes /dssg and /union locations. Prefer the path
    # recorded for each sample when its target exists, including symlink targets.
    if os.path.exists(p):
        return p

    prefix_map = config.get("path_prefix_replace") or RP_C.get("path_prefix_replace", {})
    if isinstance(prefix_map, dict):
        for src_pre, dst_pre in prefix_map.items():
            src_pre = str(src_pre).rstrip("/")
            dst_pre = str(dst_pre).rstrip("/")
            if src_pre and (p == src_pre or p.startswith(src_pre + "/")):
                alt = dst_pre + p[len(src_pre):]
                if os.path.exists(alt):
                    return alt
            elif dst_pre and (p == dst_pre or p.startswith(dst_pre + "/")):
                alt = src_pre + p[len(dst_pre):]
                if os.path.exists(alt):
                    return alt

    # Automatic bidirectional fallback between /dssg/home and /union/home.
    if "/dssg/home/" in p:
        alt = p.replace("/dssg/home/", "/union/home/", 1)
        if os.path.exists(alt):
            return alt
    elif "/union/home/" in p:
        alt = p.replace("/union/home/", "/dssg/home/", 1)
        if os.path.exists(alt):
            return alt

    # Preserve the TSV path in a missing-input error instead of inventing a
    # location when neither namespace contains the file.
    return p


def rp_resolve_write_store(path):
    """Resolve a sample store to the configured writable storage namespace."""
    if not path:
        return ""
    p = rp_full_path(path)
    prefix_map = RP_C.get("write_path_prefix_replace", {})
    if isinstance(prefix_map, dict):
        for src_pre, dst_pre in prefix_map.items():
            src_pre = str(src_pre).rstrip("/")
            dst_pre = str(dst_pre).rstrip("/")
            if not src_pre or not dst_pre:
                continue
            if p == src_pre or p.startswith(src_pre + "/"):
                alt = dst_pre + p[len(src_pre):]
                # This mapping explicitly selects the writable namespace;
                # the sample directory may not exist before fastp creates it.
                return os.path.normpath(alt)
    return p


def rp_env(name):
    env_name = config.get("environments", {}).get(name) or RP_C.get("conda_envs", {}).get(name, name)
    env_dir = config.get("conda_env_dir", "") or config.get("paths", {}).get("conda_envs", "")
    if not env_dir or os.path.isabs(env_name):
        return env_name
    candidate = os.path.join(env_dir, env_name)
    if os.path.exists(candidate):
        return candidate
    return env_name


RP_RES = rp_full_path(RP_C.get("results_dir", "results/MetaG-reads/sylph"))
RP_LOGS = rp_full_path(RP_C.get("logs_dir", "logs/read_profiling"))
RP_SKIP_MISSING = RP_C.get("skip_missing_inputs", False)
RP_INPUT_AUDIT = rp_full_path(RP_C.get("input_audit_file", os.path.join(RP_LOGS, "input_status.tsv")))
RP_CORRUPT_REPORT = rp_full_path(RP_C.get("corrupted_report_file", os.path.join(RP_LOGS, "corrupted_samples.tsv")))

# 2. Datasets and Samples Registration
RP_DATASETS = dict(config.get("datasets") or {})
RP_SAMPLES = dict(config.get("samples") or {})

# Auto-parse TSV sample sheet (e.g. RUN_id \t Store_address) if provided
_data_samples = config.get("data", {}).get("samples")
_is_tsv = isinstance(_data_samples, str) and (
    _data_samples.endswith(".tsv") or _data_samples.endswith(".txt")
)

SAMPLE_TSV_PATH = (
    (_data_samples if _is_tsv else None)
    or config.get("sample_tsv")
    or config.get("samples_tsv")
    or config.get("sample_table")
    or config.get("data", {}).get("sample_tsv")
    or config.get("data", {}).get("samples_tsv")
    or RP_C.get("sample_tsv")
    or RP_C.get("samples_tsv")
)
if SAMPLE_TSV_PATH:
    # The mapping describes read storage, not the location of the sample sheet.
    resolved_tsv = rp_full_path(SAMPLE_TSV_PATH)
    if not os.path.exists(resolved_tsv):
        for candidate_tsv in [
            os.path.join(CONF_DIR, str(SAMPLE_TSV_PATH)) if "CONF_DIR" in globals() else "",
            os.path.join(WORKFLOW_ROOT, str(SAMPLE_TSV_PATH)) if "WORKFLOW_ROOT" in globals() else "",
            os.path.join(RP_ROOT, str(SAMPLE_TSV_PATH)) if RP_ROOT else "",
            os.path.join(os.getcwd(), str(SAMPLE_TSV_PATH)),
            os.path.abspath(str(SAMPLE_TSV_PATH)),
            resolved_tsv.replace("/dssg/home/", "/union/home/", 1) if "/dssg/home/" in resolved_tsv else "",
            resolved_tsv.replace("/union/home/", "/dssg/home/", 1) if "/union/home/" in resolved_tsv else "",
        ]:
            if candidate_tsv and os.path.exists(candidate_tsv):
                resolved_tsv = candidate_tsv
                break
    if os.path.exists(resolved_tsv):
        with open(resolved_tsv, "r", encoding="utf-8") as f:
            lines = [l.strip().split("\t") for l in f if l.strip() and not l.startswith("#")]
        if lines:
            header = [h.strip() for h in lines[0]]
            run_col = 0
            store_col = 1
            has_header = False
            for idx, h in enumerate(header):
                if h.lower() in ("run_id", "sample", "sample_id", "id", "sampleid"):
                    run_col = idx
                    has_header = True
                elif h.lower() in ("store_address", "store_addr", "address", "path", "dir", "directory"):
                    store_col = idx
                    has_header = True
            data_rows = lines[1:] if has_header else lines
            for row in data_rows:
                if len(row) > max(run_col, store_col):
                    rid = row[run_col].strip()
                    saddr = row[store_col].strip()
                    if rid and saddr:
                        saddr_resolved = rp_resolve_write_store(saddr)
                        r1_candidate = rp_resolve_path(os.path.join(saddr_resolved, "00.QC", f"{rid}.clean.1.fq.gz"))
                        r2_candidate = rp_resolve_path(os.path.join(saddr_resolved, "00.QC", f"{rid}.clean.2.fq.gz"))
                        if rid not in RP_SAMPLES:
                            RP_SAMPLES[rid] = {
                                "short_read": {
                                    "r1": r1_candidate,
                                    "r2": r2_candidate,
                                },
                                "store_address": saddr_resolved,
                                "tsv_store_address": saddr,
                            }
    elif RP_SKIP_MISSING:
        raise FileNotFoundError(f"Read profiling sample TSV not found: {SAMPLE_TSV_PATH}")

# If user provided a dataset directly in read_profiling config or legacy format
if not RP_DATASETS and RP_C.get("input_dir"):
    batch_name = RP_C.get("dataset_name", "default_dataset")
    RP_DATASETS = {
        batch_name: {
            "input_dir": RP_C.get("input_dir"),
            "r1_pattern": RP_C.get("r1_pattern", "*/00.QC/*.clean.1.fq.gz"),
            "r2_pattern": RP_C.get("r2_pattern", "*/00.QC/*.clean.2.fq.gz"),
            "output_file": RP_C.get("output_file", f"profiling.{batch_name}.tsv"),
        }
    }

READ_PROFILING_TARGETS = []

# Dataset / Batch Targets (Matches collective profiling like OSMA_3nd_4and5_20241024)
for dataset_name, d_info in RP_DATASETS.items():
    out_file = d_info.get("output_file") or f"profiling.{dataset_name}.tsv"
    READ_PROFILING_TARGETS.append(os.path.join(RP_RES, out_file))

def find_sample_raw_reads(store_dir, sample_id):
    """
    Locates paired raw read files (_1/_2 or _R1/_R2) for a given sample_id
    within store_dir, accounting for delimiter differences (e.g. '_' vs '-').
    """
    if not store_dir or not sample_id:
        return None, None
    norm_sample = re.sub(r"[-_]", "", sample_id).lower()
    stores = []
    for candidate in (store_dir, rp_resolve_path(store_dir)):
        if candidate:
            candidate = os.path.normpath(rp_full_path(candidate))
            if candidate not in stores:
                stores.append(candidate)

    # Search both namespaces even when the writable store directory already
    # exists but its raw reads are still only present in the read-only store.
    for prefix in ("/dssg/home/", "/union/home/"):
        for other in ("/union/home/", "/dssg/home/"):
            for candidate in list(stores):
                if candidate.startswith(prefix):
                    alias = candidate.replace(prefix, other, 1)
                    if alias not in stores:
                        stores.append(alias)

    for resolved_store in stores:
        if not os.path.isdir(resolved_store):
            continue
        search_dirs = [resolved_store]
        try:
            sub_entries = os.listdir(resolved_store)
        except OSError:
            sub_entries = []

        for sub in sub_entries:
            sub_path = os.path.join(resolved_store, sub)
            if (
                os.path.isdir(sub_path)
                and not sub.startswith(".")
                and sub not in ("00.QC", "05.Kraken2", ".snakemake", "logs", "results")
            ):
                search_dirs.append(sub_path)

        for d in search_dirs:
            try:
                entries = os.listdir(d)
            except OSError:
                continue
            candidates = [
                f for f in entries
                if f.endswith((".fq.gz", ".fastq.gz"))
                and not f.startswith(".")
                and ".clean." not in f
            ]
            for f in candidates:
                m = re.search(r"^(.*?)(?:[._-](R?1))([._-].*|\.f(?:ast)?q\.gz)$", f, re.IGNORECASE)
                if not m:
                    continue
                prefix, r1_tag, suffix = m.groups()
                norm_prefix = re.sub(r"[-_]", "", prefix).lower()
                if norm_prefix == norm_sample:
                    r2_tag = r1_tag.replace("1", "2").replace("r", "r").replace("R", "R")
                    r2_file = prefix + f[len(prefix):].replace(r1_tag, r2_tag, 1)
                    r1_path = os.path.join(d, f)
                    r2_path = os.path.join(d, r2_file)
                    if os.path.isfile(r1_path) and os.path.isfile(r2_path):
                        return r1_path, r2_path

    return None, None


# Individual Sample Targets (if any)
VALID_SAMPLE_KEYS = []
RP_INPUT_AUDIT_ROWS = []
if isinstance(RP_SAMPLES, dict):
    for s_name, s_info in RP_SAMPLES.items():
        if isinstance(s_info, dict) and ("short_read" in s_info or "r1" in s_info):
            sr = s_info.get("short_read", s_info)
            r1_path = rp_resolve_path(sr.get("r1") or sr.get("external_clean_r1") or "")
            r2_path = rp_resolve_path(sr.get("r2") or sr.get("external_clean_r2") or "")
            r1_ok = os.path.isfile(r1_path)
            r2_ok = os.path.isfile(r2_path)
            store = s_info.get("store_address", "")
            qc_dir = os.path.join(store, "00.QC") if store else ""
            qc_target = os.readlink(qc_dir) if qc_dir and os.path.islink(qc_dir) else ""

            raw_r1, raw_r2 = (None, None)
            raw_r1_path, raw_r2_path = ("", "")
            corrupt_marker = os.path.join(qc_dir, f"{s_name}.corrupted") if qc_dir else ""
            corrupt_txt = os.path.join(qc_dir, "CORRUPTED.txt") if qc_dir else ""
            is_corrupt_qc = bool(
                (corrupt_marker and os.path.isfile(corrupt_marker))
                or (corrupt_txt and os.path.isfile(corrupt_txt))
            )

            if r1_ok and r2_ok:
                if is_corrupt_qc and (os.path.getsize(r1_path) <= 30 or os.path.getsize(r2_path) <= 30):
                    status = "corrupted_qc"
                else:
                    status = "ready_clean"
            else:
                raw_r1, raw_r2 = find_sample_raw_reads(store, s_name)
                if raw_r1 and raw_r2:
                    raw_r1_path = rp_resolve_path(raw_r1)
                    raw_r2_path = rp_resolve_path(raw_r2)
                    if os.path.isfile(raw_r1_path) and os.path.isfile(raw_r2_path):
                        if is_corrupt_qc:
                            status = "corrupted_raw"
                        else:
                            status = "ready_needs_qc"
                            s_info["raw_reads"] = {"r1": raw_r1_path, "r2": raw_r2_path}
                            # Keep the downstream Sylph dependency on the same
                            # writable QC targets produced by fastp. The input
                            # paths above may resolve to an existing /union copy,
                            # but that namespace is read-only on compute nodes.
                            qc_r1 = os.path.join(store, "00.QC", f"{s_name}.clean.1.fq.gz")
                            qc_r2 = os.path.join(store, "00.QC", f"{s_name}.clean.2.fq.gz")
                            if "short_read" in s_info:
                                s_info["short_read"]["r1"] = qc_r1
                                s_info["short_read"]["r2"] = qc_r2
                            else:
                                s_info["r1"] = qc_r1
                                s_info["r2"] = qc_r2
                    else:
                        status = "missing_raw"
                else:
                    status = (
                        "missing_both"
                        if not r1_ok and not r2_ok
                        else ("missing_r1" if not r1_ok else "missing_r2")
                    )

            RP_INPUT_AUDIT_ROWS.append({
                "RUN_id": s_name,
                "Store_address": s_info.get("tsv_store_address", store),
                "resolved_store_address": store,
                "qc_dir": qc_dir,
                "qc_link_target": qc_target,
                "r1_path": r1_path,
                "r1_exists": r1_ok,
                "r2_path": r2_path,
                "r2_exists": r2_ok,
                "raw_r1_path": raw_r1_path,
                "raw_r2_path": raw_r2_path,
                "status": status,
            })

            if RP_SKIP_MISSING and status not in ("ready", "ready_clean", "ready_needs_qc"):
                continue

            VALID_SAMPLE_KEYS.append(s_name)
            READ_PROFILING_TARGETS.append(os.path.join(RP_RES, s_name, f"{s_name}.sylph.tsv"))

if RP_SKIP_MISSING:
    # Written while constructing the DAG so a dry-run also records every choice.
    os.makedirs(os.path.dirname(RP_INPUT_AUDIT), exist_ok=True)
    audit_columns = [
        "RUN_id", "Store_address", "resolved_store_address", "qc_dir",
        "qc_link_target", "r1_path", "r1_exists", "r2_path",
        "r2_exists", "raw_r1_path", "raw_r2_path", "status",
    ]
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=os.path.dirname(RP_INPUT_AUDIT),
        prefix=".input_status.", suffix=".tsv", delete=False,
    ) as audit_handle:
        writer = csv.DictWriter(audit_handle, fieldnames=audit_columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(RP_INPUT_AUDIT_ROWS)
        audit_tmp = audit_handle.name
    os.replace(audit_tmp, RP_INPUT_AUDIT)
    skipped = len(RP_INPUT_AUDIT_ROWS) - len(VALID_SAMPLE_KEYS)
    ready_clean_count = sum(1 for r in RP_INPUT_AUDIT_ROWS if r.get("status") in ("ready", "ready_clean"))
    ready_qc_count = sum(1 for r in RP_INPUT_AUDIT_ROWS if r.get("status") == "ready_needs_qc")
    corrupt_count = sum(1 for r in RP_INPUT_AUDIT_ROWS if str(r.get("status", "")).startswith("corrupt"))
    corrupt_msg = f", {corrupt_count} corrupted" if corrupt_count > 0 else ""
    print(
        f"Read profiling inputs: {len(VALID_SAMPLE_KEYS)} active "
        f"({ready_clean_count} clean reads ready, {ready_qc_count} scheduled for fastp QC{corrupt_msg}), "
        f"{skipped} skipped; audit: {RP_INPUT_AUDIT}"
    )

if len(VALID_SAMPLE_KEYS) > 1:
    READ_PROFILING_TARGETS.append(os.path.join(RP_RES, "merged_sylph_profiling.tsv"))

# Sylph-tax taxonomic profiling targets (sylph-tax taxprof / merge)
RP_RUN_TAXONOMY = bool(RP_C.get("run_taxonomy", True) and RP_C.get("taxonomy"))
RP_TAX_PREFIX = str(RP_C.get("taxonomy_prefix", "gtdb_r226")).strip()

if RP_RUN_TAXONOMY:
    for s_name in VALID_SAMPLE_KEYS:
        READ_PROFILING_TARGETS.append(
            os.path.join(RP_RES, s_name, f"{RP_TAX_PREFIX}_{s_name}.sylphmpa")
        )
    if len(VALID_SAMPLE_KEYS) > 1:
        READ_PROFILING_TARGETS.append(
            os.path.join(RP_RES, f"merged_{RP_TAX_PREFIX}_relative_abundance.tsv")
        )
        READ_PROFILING_TARGETS.append(
            os.path.join(RP_RES, f"merged_{RP_TAX_PREFIX}_sequence_abundance.tsv")
        )
    for dataset_name in RP_DATASETS.keys():
        READ_PROFILING_TARGETS.append(
            os.path.join(RP_RES, f"{RP_TAX_PREFIX}_{dataset_name}.sylphmpa")
        )


# Helpers for Rules
def get_batch_input_dir(dataset):
    d_info = RP_DATASETS.get(dataset, {})
    return rp_resolve_path(d_info.get("input_dir", ""))


def get_batch_r1_pattern(dataset):
    d_info = RP_DATASETS.get(dataset, {})
    return d_info.get("r1_pattern", "*/00.QC/*.clean.1.fq.gz")


def get_batch_r2_pattern(dataset):
    d_info = RP_DATASETS.get(dataset, {})
    return d_info.get("r2_pattern", "*/00.QC/*.clean.2.fq.gz")


def get_batch_inputs(wildcards):
    d = wildcards.dataset
    in_dir = get_batch_input_dir(d)
    r1_pat = get_batch_r1_pattern(d)
    if in_dir and os.path.isdir(in_dir):
        files = glob.glob(os.path.join(in_dir, r1_pat))
        if files:
            return sorted(files)
    return []


def get_sample_r1(sample):
    s_info = RP_SAMPLES.get(sample, {})
    sr = s_info.get("short_read", {}) if "short_read" in s_info else s_info
    r1 = sr.get("r1") or sr.get("external_clean_r1") or ""
    return rp_resolve_path(r1)


def get_sample_r2(sample):
    s_info = RP_SAMPLES.get(sample, {})
    sr = s_info.get("short_read", {}) if "short_read" in s_info else s_info
    r2 = sr.get("r2") or sr.get("external_clean_r2") or ""
    return rp_resolve_path(r2)


def get_sample_raw_r1(wildcards):
    s_name = getattr(wildcards, "sample", str(wildcards))
    s_info = RP_SAMPLES.get(s_name, {})
    raw = s_info.get("raw_reads", {})
    if isinstance(raw, dict) and raw.get("r1"):
        return rp_resolve_path(raw["r1"])
    elif isinstance(raw, (list, tuple)) and len(raw) >= 1 and raw[0]:
        return rp_resolve_path(raw[0])
    store = s_info.get("store_address", "") or getattr(wildcards, "store", "")
    r1, _ = find_sample_raw_reads(store, s_name)
    return rp_resolve_path(r1) if r1 else []


def get_sample_raw_r2(wildcards):
    s_name = getattr(wildcards, "sample", str(wildcards))
    s_info = RP_SAMPLES.get(s_name, {})
    raw = s_info.get("raw_reads", {})
    if isinstance(raw, dict) and raw.get("r2"):
        return rp_resolve_path(raw["r2"])
    elif isinstance(raw, (list, tuple)) and len(raw) >= 2 and raw[1]:
        return rp_resolve_path(raw[1])
    store = s_info.get("store_address", "") or getattr(wildcards, "store", "")
    _, r2 = find_sample_raw_reads(store, s_name)
    return rp_resolve_path(r2) if r2 else []


# ------------------------------------------------------------------------------
# 3. Rules
# ------------------------------------------------------------------------------

rule read_profiling_sylph_batch:
    input:
        reads = get_batch_inputs
    output:
        tsv = os.path.join(RP_RES, "profiling.{dataset}.tsv")
    log:
        os.path.join(RP_LOGS, "batch/{dataset}.log")
    threads: get_threads(RP_C.get("batch_threads", RP_C.get("threads", 60)))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = rp_env("sylph"),
        db = rp_full_path(RP_C.get("database", "")),
        input_dir = lambda w: get_batch_input_dir(w.dataset),
        r1_pattern = lambda w: get_batch_r1_pattern(w.dataset),
        r2_pattern = lambda w: get_batch_r2_pattern(w.dataset),
        extra = RP_C.get("params", "")
    wildcard_constraints:
        dataset = r"[^/]+"
    shell:
        """
        set -euo pipefail

        mkdir -p "$(dirname "{output.tsv}")"
        mkdir -p "$(dirname "{log}")"

        exec > "{log}" 2>&1
        set -x

        ulimit -s unlimited || true
        ulimit -l unlimited || true

        if [ -x "{params.env}/bin/sylph" ]; then
            SYLPH_CMD="{params.env}/bin/sylph"
        elif command -v sylph >/dev/null 2>&1; then
            SYLPH_CMD="$(command -v sylph)"
        elif [ -d "{params.env}" ]; then
            SYLPH_CMD="conda run -p {params.env} sylph"
        else
            SYLPH_CMD="conda run -n {params.env} sylph"
        fi

        cd "{params.input_dir}"

        echo "=== Starting Sylph Batch Taxonomic Profiling ==="
        echo "Dataset: {wildcards.dataset}"
        echo "Input Directory: {params.input_dir}"
        echo "Database: {params.db}"
        echo "Threads: {threads}"
        echo "Date: $(date)"
        echo "Sylph Command: $SYLPH_CMD"

        $SYLPH_CMD profile \
            "{params.db}" \
            -1 {params.r1_pattern} \
            -2 {params.r2_pattern} \
            -t {threads} \
            {params.extra} > "{output.tsv}"

        echo "=== Completed Sylph Batch Profiling at $(date) ==="
        """


rule read_profiling_qc_fastp:
    input:
        r1_raw = lambda w: get_sample_raw_r1(w),
        r2_raw = lambda w: get_sample_raw_r2(w)
    output:
        r1_clean = "{store}/00.QC/{sample}.clean.1.fq.gz",
        r2_clean = "{store}/00.QC/{sample}.clean.2.fq.gz"
    log:
        os.path.join(RP_LOGS, "fastp{store}/{sample}.log")
    threads: get_threads(RP_C.get("qc_threads", 16))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = rp_env("qc"),
        adapter1 = RP_C.get("qc_adapter1", "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"),
        adapter2 = RP_C.get("qc_adapter2", "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"),
        length_required = RP_C.get("qc_length_required", 100),
        extra = RP_C.get("qc_params", "--cut_front --cut_right --cut_window_size 4 --cut_mean_quality 20"),
        path_out = lambda w: os.path.join(w.store, "00.QC"),
        fq_merged = lambda w: os.path.join(w.store, "00.QC", f"{w.sample}.merged.fq.gz"),
        fq1_unmerged = lambda w: os.path.join(w.store, "00.QC", f"{w.sample}.unmerged.1.fq.gz"),
        fq2_unmerged = lambda w: os.path.join(w.store, "00.QC", f"{w.sample}.unmerged.2.fq.gz"),
        json_m = lambda w: os.path.join(w.store, "00.QC", f"{w.sample}.fastp_m.json"),
        html_m = lambda w: os.path.join(w.store, "00.QC", f"{w.sample}.fastp_m.html"),
        json_p = lambda w: os.path.join(w.store, "00.QC", f"{w.sample}.fastp_p.json"),
        html_p = lambda w: os.path.join(w.store, "00.QC", f"{w.sample}.fastp_p.html"),
        corrupt_report = RP_CORRUPT_REPORT,
        precheck_gzip = bool(RP_C.get("qc_precheck_gzip", False))
    wildcard_constraints:
        store = r".+",
        sample = r"[^/]+"
    shell:
        """
        set -euo pipefail

        mkdir -p "$(dirname "{output.r1_clean}")"
        mkdir -p "$(dirname "{log}")"

        exec > "{log}" 2>&1
        set -x

        ulimit -s unlimited || true
        ulimit -l unlimited || true

        for conda_sh in \
            /dssg/home/acct-schhui/gses0/Soft/Conda/etc/profile.d/conda.sh \
            "$HOME/miniconda3/etc/profile.d/conda.sh" \
            "$HOME/anaconda3/etc/profile.d/conda.sh"; do
            if [ -f "$conda_sh" ]; then
                source "$conda_sh"
                break
            fi
        done

        if [ -d "{params.env}/bin" ]; then
            export PATH="{params.env}/bin:$PATH"
        fi

        if command -v conda >/dev/null 2>&1; then
            conda activate "{params.env}" 2>/dev/null || conda activate trench 2>/dev/null || true
        fi

        if [ -x "{params.env}/bin/fastp" ]; then
            FASTP_CMD="{params.env}/bin/fastp"
        elif [ -x "/dssg/home/acct-schhui/gses0/Soft/Conda/envs/trench/bin/fastp" ]; then
            FASTP_CMD="/dssg/home/acct-schhui/gses0/Soft/Conda/envs/trench/bin/fastp"
        elif command -v fastp >/dev/null 2>&1; then
            FASTP_CMD="$(command -v fastp)"
        elif [ -d "{params.env}" ]; then
            FASTP_CMD="conda run -p {params.env} fastp"
        else
            FASTP_CMD="conda run -n {params.env} fastp"
        fi

        echo "=== Starting fastp QC: {wildcards.sample} ==="
        echo "Date: $(date)"
        echo "Host: $(hostname)"
        echo "Fastp Command: $FASTP_CMD"
        echo "Raw R1: {input.r1_raw}"
        echo "Raw R2: {input.r2_raw}"
        echo "Output Dir: {params.path_out}"
        echo "Threads: {threads}"

        CORRUPT_MARKER="{params.path_out}/{wildcards.sample}.corrupted"
        CORRUPT_TXT="{params.path_out}/CORRUPTED.txt"

        # Check if already completed and marked as corrupted previously
        if [ -f "$CORRUPT_MARKER" ] || [ -f "$CORRUPT_TXT" ]; then
            echo "Notice: Sample {wildcards.sample} was previously detected as corrupted. Skipping fastp execution."
            if [ ! -s "{output.r1_clean}" ]; then
                gzip -n -c < /dev/null > "{output.r1_clean}"
            fi
            if [ ! -s "{output.r2_clean}" ]; then
                gzip -n -c < /dev/null > "{output.r2_clean}"
            fi
            exit 0
        fi

        if [ -s "{output.r1_clean}" ] && [ -s "{output.r2_clean}" ] && [ -f "{params.path_out}/done" ]; then
            echo "Clean reads already exist and non-empty. Skipping fastp execution."
            exit 0
        fi

        # Graceful corruption reporting handler
        handle_corrupted_sample() {{
            local err_type="$1"
            local err_msg="$2"

            echo "================================================================================"
            echo "WARNING: SAMPLE GZ CORRUPTION DETECTED FOR {wildcards.sample}"
            echo "Error Type:   $err_type"
            echo "Error Detail: $err_msg"
            echo "Action:       Recording corruption to report; workflow will NOT be interrupted."
            echo "================================================================================"

            # Clean up partial/broken files
            rm -f "{output.r1_clean}" "{output.r2_clean}" "{params.fq_merged}" "{params.fq1_unmerged}" "{params.fq2_unmerged}"

            # Create valid empty gzip placeholders (20 bytes) so Snakemake does not fail MissingOutputException
            gzip -n -c < /dev/null > "{output.r1_clean}"
            gzip -n -c < /dev/null > "{output.r2_clean}"

            # Create dummy JSON and HTML reports
            printf '{{\n  "error": "corrupted_input",\n  "sample": "%s",\n  "error_type": "%s",\n  "error_message": "%s"\n}}\n' \
                "{wildcards.sample}" "$err_type" "$err_msg" > "{params.json_p}"
            printf '<html><head><title>Corrupted Sample %s</title></head><body><h1>Corrupted Sample: %s</h1><p><b>Error Type:</b> %s</p><p><b>Detail:</b> %s</p></body></html>\n' \
                "{wildcards.sample}" "{wildcards.sample}" "$err_type" "$err_msg" > "{params.html_p}"
            [ -f "{params.json_m}" ] || cp -f "{params.json_p}" "{params.json_m}"
            [ -f "{params.html_m}" ] || cp -f "{params.html_p}" "{params.html_m}"

            # Create sample-level marker files in 00.QC directory
            local now_str
            now_str="$(date -u +"%Y-%m-%dT%H:%M:%SZ" 2>/dev/null || date)"
            printf "Sample:\t%s\nError_type:\t%s\nError_detail:\t%s\nRaw_R1:\t%s\nRaw_R2:\t%s\nTimestamp:\t%s\n" \
                "{wildcards.sample}" "$err_type" "$err_msg" "{input.r1_raw}" "{input.r2_raw}" "$now_str" > "$CORRUPT_MARKER"
            printf "%s: %s\n" "$err_type" "$err_msg" > "$CORRUPT_TXT"

            # Atomically record to central corrupted report file
            local report_file="{params.corrupt_report}"
            mkdir -p "$(dirname "$report_file")"
            (
                if command -v flock >/dev/null 2>&1; then
                    flock -x 200
                fi
                if [ ! -f "$report_file" ] || [ ! -s "$report_file" ]; then
                    printf "RUN_id\tStore_address\tRaw_R1\tRaw_R2\tError_type\tError_detail\tTimestamp\n" > "$report_file"
                fi
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                    "{wildcards.sample}" \
                    "{wildcards.store}" \
                    "{input.r1_raw}" \
                    "{input.r2_raw}" \
                    "$err_type" \
                    "$err_msg" \
                    "$now_str" >> "$report_file"
            ) 200>"$report_file.lock"

            touch "{params.path_out}/done"
            echo "=== Completed fastp QC (corrupted sample recorded) at $(date) ==="
            exit 0
        }}

        if [ -d "{wildcards.store}/05.Kraken2" ] && ls "{wildcards.store}/05.Kraken2"/*fq.gz >/dev/null 2>&1; then
            echo "Found existing Kraken2 reads, moving to 00.QC..."
            mv "{wildcards.store}/05.Kraken2"/*fq.gz "{params.path_out}/"
            [ -f "{wildcards.store}/05.Kraken2/{wildcards.sample}.fastp.html" ] && mv "{wildcards.store}/05.Kraken2/{wildcards.sample}.fastp.html" "{params.html_p}" || true
            [ -f "{wildcards.store}/05.Kraken2/{wildcards.sample}.fastp.json" ] && mv "{wildcards.store}/05.Kraken2/{wildcards.sample}.fastp.json" "{params.json_p}" || true
            [ -f "{params.path_out}/{wildcards.sample}.fastp.json" ] && mv "{params.path_out}/{wildcards.sample}.fastp.json" "{params.json_m}" || true
            [ -f "{params.path_out}/{wildcards.sample}.fastp.html" ] && mv "{params.path_out}/{wildcards.sample}.fastp.html" "{params.html_m}" || true
        else
            # 1. Check raw file existence and non-zero size
            if [ ! -f "{input.r1_raw}" ]; then
                handle_corrupted_sample "MissingRawR1" "Raw R1 file not found: {input.r1_raw}"
            elif [ ! -s "{input.r1_raw}" ]; then
                handle_corrupted_sample "EmptyRawR1" "Raw R1 file is empty (0 bytes): {input.r1_raw}"
            fi

            if [ ! -f "{input.r2_raw}" ]; then
                handle_corrupted_sample "MissingRawR2" "Raw R2 file not found: {input.r2_raw}"
            elif [ ! -s "{input.r2_raw}" ]; then
                handle_corrupted_sample "EmptyRawR2" "Raw R2 file is empty (0 bytes): {input.r2_raw}"
            fi

            # 2. Optional fast pre-check if enabled in config
            if [ "{params.precheck_gzip}" = "True" ] || [ "{params.precheck_gzip}" = "true" ]; then
                if ! gzip -t "{input.r1_raw}" >/dev/null 2>&1; then
                    GZ1_TEST=$(gzip -t "{input.r1_raw}" 2>&1 || true)
                    handle_corrupted_sample "CorruptedGzipR1" "gzip -t integrity check failed: $GZ1_TEST"
                fi
                if ! gzip -t "{input.r2_raw}" >/dev/null 2>&1; then
                    GZ2_TEST=$(gzip -t "{input.r2_raw}" 2>&1 || true)
                    handle_corrupted_sample "CorruptedGzipR2" "gzip -t integrity check failed: $GZ2_TEST"
                fi
            fi

            # 3. Run fastp merged (with error interception)
            echo "Running fastp merged..."
            set +e
            MERGED_ERR=$($FASTP_CMD -i "{input.r1_raw}" -I "{input.r2_raw}" \
                --out1 "{params.fq1_unmerged}" --out2 "{params.fq2_unmerged}" \
                -m --merged_out "{params.fq_merged}" \
                --adapter_sequence "{params.adapter1}" --adapter_sequence_r2 "{params.adapter2}" \
                --length_required "{params.length_required}" -w {threads} \
                -j "{params.json_m}" -h "{params.html_m}" \
                {params.extra} 2>&1)
            MERGED_RC=$?
            set -e

            if [ "$MERGED_RC" -ne 0 ]; then
                echo "fastp merged failed with return code $MERGED_RC. Checking gzip integrity..."
                GZ1_TEST=$(gzip -t "{input.r1_raw}" 2>&1 || true)
                GZ2_TEST=$(gzip -t "{input.r2_raw}" 2>&1 || true)
                if [ -n "$GZ1_TEST" ] && ! gzip -t "{input.r1_raw}" >/dev/null 2>&1; then
                    handle_corrupted_sample "CorruptedGzipR1" "gzip -t failed on {input.r1_raw}: $GZ1_TEST"
                elif [ -n "$GZ2_TEST" ] && ! gzip -t "{input.r2_raw}" >/dev/null 2>&1; then
                    handle_corrupted_sample "CorruptedGzipR2" "gzip -t failed on {input.r2_raw}: $GZ2_TEST"
                else
                    handle_corrupted_sample "FastpMergedError" "fastp merged failed (rc=$MERGED_RC): $MERGED_ERR"
                fi
            fi

            # 4. Run fastp paired-end clean (with error interception)
            echo "Running fastp paired-end clean..."
            set +e
            CLEAN_ERR=$($FASTP_CMD -i "{input.r1_raw}" -I "{input.r2_raw}" \
                --out1 "{output.r1_clean}" --out2 "{output.r2_clean}" \
                --adapter_sequence "{params.adapter1}" --adapter_sequence_r2 "{params.adapter2}" \
                --length_required "{params.length_required}" -w {threads} \
                -j "{params.json_p}" -h "{params.html_p}" \
                {params.extra} 2>&1)
            CLEAN_RC=$?
            set -e

            if [ "$CLEAN_RC" -ne 0 ]; then
                echo "fastp paired-end clean failed with return code $CLEAN_RC. Checking gzip integrity..."
                GZ1_TEST=$(gzip -t "{input.r1_raw}" 2>&1 || true)
                GZ2_TEST=$(gzip -t "{input.r2_raw}" 2>&1 || true)
                if [ -n "$GZ1_TEST" ] && ! gzip -t "{input.r1_raw}" >/dev/null 2>&1; then
                    handle_corrupted_sample "CorruptedGzipR1" "gzip -t failed on {input.r1_raw}: $GZ1_TEST"
                elif [ -n "$GZ2_TEST" ] && ! gzip -t "{input.r2_raw}" >/dev/null 2>&1; then
                    handle_corrupted_sample "CorruptedGzipR2" "gzip -t failed on {input.r2_raw}: $GZ2_TEST"
                else
                    handle_corrupted_sample "FastpCleanError" "fastp clean failed (rc=$CLEAN_RC): $CLEAN_ERR"
                fi
            fi
        fi

        touch "{params.path_out}/done"
        echo "=== Completed fastp QC at $(date) ==="
        """



rule read_profiling_sylph_sample:
    input:
        r1 = lambda w: get_sample_r1(w.sample),
        r2 = lambda w: get_sample_r2(w.sample)
    output:
        tsv = os.path.join(RP_RES, "{sample}/{sample}.sylph.tsv")
    log:
        os.path.join(RP_LOGS, "sylph/{sample}.log")
    threads: get_threads(RP_C.get("sample_threads", RP_C.get("threads", 8)))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = rp_env("sylph"),
        db = rp_full_path(RP_C.get("database", "")),
        extra = RP_C.get("params", "")
    wildcard_constraints:
        sample = r"[^/]+"
    shell:
        """
        set -euo pipefail

        mkdir -p "$(dirname "{output.tsv}")"
        mkdir -p "$(dirname "{log}")"

        exec > "{log}" 2>&1
        set -x

        ulimit -s unlimited || true
        ulimit -l unlimited || true

        CORRUPT_MARKER="$(dirname "{input.r1}")/{wildcards.sample}.corrupted"
        CORRUPT_TXT="$(dirname "{input.r1}")/CORRUPTED.txt"

        # Check if sample was marked as corrupted or has empty reads from fastp QC
        IS_SAMPLE_CORRUPTED=0
        if [ -f "$CORRUPT_MARKER" ] || [ -f "$CORRUPT_TXT" ]; then
            IS_SAMPLE_CORRUPTED=1
        elif [ ! -s "{input.r1}" ] || [ ! -s "{input.r2}" ]; then
            IS_SAMPLE_CORRUPTED=1
        elif [ $(gzip -dc "{input.r1}" 2>/dev/null | head -c 10 | wc -c) -eq 0 ]; then
            IS_SAMPLE_CORRUPTED=1
        fi

        if [ "$IS_SAMPLE_CORRUPTED" -eq 1 ]; then
            echo "================================================================================"
            echo "NOTICE: Sample {wildcards.sample} was corrupted/empty during QC."
            echo "Skipping sylph profile execution and writing empty profile placeholder."
            echo "================================================================================"
            printf "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_marker_coverage\tEff_cov\tTrue_cov\tAni\tContig_name\n" > "{output.tsv}"
            echo "=== Completed Sylph Sample Profiling placeholder at $(date) ==="
            exit 0
        fi

        if [ -x "{params.env}/bin/sylph" ]; then
            SYLPH_CMD="{params.env}/bin/sylph"
        elif command -v sylph >/dev/null 2>&1; then
            SYLPH_CMD="$(command -v sylph)"
        elif [ -d "{params.env}" ]; then
            SYLPH_CMD="conda run -p {params.env} sylph"
        else
            SYLPH_CMD="conda run -n {params.env} sylph"
        fi

        echo "=== Starting Sylph Sample Profiling: {wildcards.sample} ==="
        echo "Date: $(date)"
        echo "Host: $(hostname)"
        echo "Sylph Command: $SYLPH_CMD"
        echo "Database: {params.db}"
        echo "R1: {input.r1}"
        echo "R2: {input.r2}"
        echo "Threads: {threads}"

        if [ ! -s "{input.r1}" ]; then
            echo "ERROR: Input R1 does not exist or is empty (0 bytes): {input.r1}" >&2
            exit 1
        fi
        if [ ! -s "{input.r2}" ]; then
            echo "ERROR: Input R2 does not exist or is empty (0 bytes): {input.r2}" >&2
            exit 1
        fi
        if [ ! -f "{params.db}" ]; then
            echo "ERROR: Sylph database file does not exist: {params.db}" >&2
            exit 1
        fi

        $SYLPH_CMD profile \
            "{params.db}" \
            -1 "{input.r1}" \
            -2 "{input.r2}" \
            -t {threads} \
            {params.extra} > "{output.tsv}"

        echo "=== Completed Sylph Sample Profiling at $(date) ==="
        """


rule read_profiling_sylph_merge:
    input:
        tsvs = lambda wildcards: [os.path.join(RP_RES, s, f"{s}.sylph.tsv") for s in VALID_SAMPLE_KEYS]
    output:
        merged = os.path.join(RP_RES, "merged_sylph_profiling.tsv")
    log:
        os.path.join(RP_LOGS, "sylph_merge.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = rp_env("python")
    shell:
        """
        set -euo pipefail

        mkdir -p "$(dirname "{output.merged}")"
        mkdir -p "$(dirname "{log}")"

        exec > "{log}" 2>&1
        set -x

        first=1
        touch "{output.merged}"
        for f in {input.tsvs}; do
            [ ! -f "$f" ] && continue
            if [ "$first" -eq 1 ]; then
                cat "$f" > "{output.merged}"
                first=0
            else
                tail -n +2 "$f" >> "{output.merged}"
            fi
        done

        if [ ! -s "{output.merged}" ]; then
            printf "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_marker_coverage\tEff_cov\tTrue_cov\tAni\tContig_name\n" > "{output.merged}"
        fi
        """


rule read_profiling_sylph_taxprof_sample:
    input:
        tsv = os.path.join(RP_RES, "{sample}/{sample}.sylph.tsv")
    output:
        sylphmpa = os.path.join(RP_RES, f"{{sample}}/{RP_TAX_PREFIX}_{{sample}}.sylphmpa")
    log:
        os.path.join(RP_LOGS, f"sylph_taxprof/{{sample}}.{RP_TAX_PREFIX}.log")
    threads: get_threads(RP_C.get("taxprof_threads", 4))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = rp_env("sylph"),
        taxonomy = rp_resolve_path(RP_C.get("taxonomy", "/dssg/home/acct-schhui/gses0/liliuyang/osma/software_from_nature/sylph/taxonomy_file_folder/gtdb_r226_metadata.tsv.gz")),
        prefix = lambda w: os.path.join(RP_RES, w.sample, f"{RP_TAX_PREFIX}_{w.sample}."),
        extra = RP_C.get("taxprof_params", "")
    wildcard_constraints:
        sample = r"[^/]+"
    shell:
        """
        set -euo pipefail

        mkdir -p "$(dirname "{output.sylphmpa}")"
        mkdir -p "$(dirname "{log}")"

        exec > "{log}" 2>&1
        set -x

        ulimit -s unlimited || true
        ulimit -l unlimited || true

        # Source conda hook if available
        for conda_sh in \
            /dssg/home/acct-schhui/gses0/Soft/Conda/etc/profile.d/conda.sh \
            "$HOME/miniconda3/etc/profile.d/conda.sh" \
            "$HOME/anaconda3/etc/profile.d/conda.sh"; do
            if [ -f "$conda_sh" ]; then
                source "$conda_sh"
                break
            fi
        done

        if [ -d "{params.env}/bin" ]; then
            export PATH="{params.env}/bin:$PATH"
        fi

        if command -v conda >/dev/null 2>&1; then
            conda activate "{params.env}" 2>/dev/null || conda activate sylph 2>/dev/null || true
        fi

        if [ -x "{params.env}/bin/python" ] && [ -f "{params.env}/bin/sylph-tax" ]; then
            SYLPH_TAX_CMD="{params.env}/bin/python {params.env}/bin/sylph-tax"
        elif [ -x "{params.env}/bin/sylph-tax" ]; then
            SYLPH_TAX_CMD="{params.env}/bin/sylph-tax"
        elif command -v sylph-tax >/dev/null 2>&1; then
            SYLPH_TAX_CMD="$(command -v sylph-tax)"
        elif [ -d "{params.env}" ]; then
            SYLPH_TAX_CMD="conda run -p {params.env} sylph-tax"
        else
            SYLPH_TAX_CMD="conda run -n {params.env} sylph-tax"
        fi

        echo "=== Starting sylph-tax Taxonomic Profiling: {wildcards.sample} ==="
        echo "Date: $(date)"
        echo "Host: $(hostname)"
        echo "Sylph-tax Command: $SYLPH_TAX_CMD"
        echo "Input TSV: {input.tsv}"
        echo "Taxonomy: {params.taxonomy}"

        # If input profiling TSV has only header (corrupted sample or 0 hits), write empty sylphmpa profile
        NUM_LINES=$(wc -l < "{input.tsv}" 2>/dev/null || echo 0)
        if [ "$NUM_LINES" -le 1 ]; then
            echo "Notice: Input profiling TSV has no taxonomic hits (lines <= 1). Creating empty sylphmpa profile for {wildcards.sample}."
            printf "#SampleID\t{wildcards.sample}\tTaxonomies_used:{params.taxonomy}\nclade_name\trelative_abundance\tsequence_abundance\tANI (if strain-level)\tCoverage (if strain-level)\n" > "{output.sylphmpa}"
            echo "=== Completed sylph-tax placeholder at $(date) ==="
            exit 0
        fi

        if [ ! -s "{input.tsv}" ]; then
            echo "ERROR: Input profiling TSV is missing or empty (0 bytes): {input.tsv}" >&2
            exit 1
        fi
        if [ ! -f "{params.taxonomy}" ]; then
            echo "ERROR: Taxonomy metadata file not found: {params.taxonomy}" >&2
            exit 1
        fi

        # Clean any stale .sylphmpa files in this sample's folder to avoid conflict
        rm -f "$(dirname "{output.sylphmpa}")"/*.sylphmpa

        $SYLPH_TAX_CMD taxprof "{input.tsv}" \
            -t "{params.taxonomy}" \
            -o "{params.prefix}" \
            {params.extra}

        # Verify output file or normalize if sylph-tax appended sample/fq name
        if [ ! -f "{output.sylphmpa}" ]; then
            MATCH=$(find "$(dirname "{output.sylphmpa}")" -maxdepth 1 -name "*.sylphmpa" ! -name "$(basename "{output.sylphmpa}")" | head -n 1)
            if [ -n "$MATCH" ] && [ -f "$MATCH" ]; then
                mv -f "$MATCH" "{output.sylphmpa}"
            fi
        fi

        # If sylph.tsv had 0 hits (no taxa detected), sylph-tax produces no file.
        # Create a valid empty CAMI / sylphmpa profile so downstream merge succeeds.
        if [ ! -f "{output.sylphmpa}" ]; then
            echo "Notice: No taxonomy matches found in {input.tsv}. Creating empty sylphmpa profile."
            printf "#SampleID\t{wildcards.sample}\tTaxonomies_used:{params.taxonomy}\nclade_name\trelative_abundance\tsequence_abundance\tANI (if strain-level)\tCoverage (if strain-level)\n" > "{output.sylphmpa}"
        else
            # Ensure the first line has the clean sample name as SampleID for downstream merge table headers
            sed -i.bak -e "1s|^#SampleID[^\t]*\t[^\t]*|#SampleID\t{wildcards.sample}|" "{output.sylphmpa}" && rm -f "{output.sylphmpa}.bak"
        fi

        echo "=== Completed sylph-tax Taxonomic Profiling at $(date) ==="
        """


rule read_profiling_sylph_taxprof_merge:
    input:
        sylphmpas = lambda wildcards: [os.path.join(RP_RES, s, f"{RP_TAX_PREFIX}_{s}.sylphmpa") for s in VALID_SAMPLE_KEYS]
    output:
        rel_abun = os.path.join(RP_RES, f"merged_{RP_TAX_PREFIX}_relative_abundance.tsv"),
        seq_abun = os.path.join(RP_RES, f"merged_{RP_TAX_PREFIX}_sequence_abundance.tsv")
    log:
        os.path.join(RP_LOGS, f"sylph_taxprof_merge_{RP_TAX_PREFIX}.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = rp_env("sylph")
    shell:
        """
        set -euo pipefail

        mkdir -p "$(dirname "{output.rel_abun}")"
        mkdir -p "$(dirname "{log}")"

        exec > "{log}" 2>&1
        set -x

        for conda_sh in \
            /dssg/home/acct-schhui/gses0/Soft/Conda/etc/profile.d/conda.sh \
            "$HOME/miniconda3/etc/profile.d/conda.sh" \
            "$HOME/anaconda3/etc/profile.d/conda.sh"; do
            if [ -f "$conda_sh" ]; then
                source "$conda_sh"
                break
            fi
        done

        if [ -d "{params.env}/bin" ]; then
            export PATH="{params.env}/bin:$PATH"
        fi

        if command -v conda >/dev/null 2>&1; then
            conda activate "{params.env}" 2>/dev/null || conda activate sylph 2>/dev/null || true
        fi

        if [ -x "{params.env}/bin/python" ] && [ -f "{params.env}/bin/sylph-tax" ]; then
            SYLPH_TAX_CMD="{params.env}/bin/python {params.env}/bin/sylph-tax"
        elif [ -x "{params.env}/bin/sylph-tax" ]; then
            SYLPH_TAX_CMD="{params.env}/bin/sylph-tax"
        elif command -v sylph-tax >/dev/null 2>&1; then
            SYLPH_TAX_CMD="$(command -v sylph-tax)"
        elif [ -d "{params.env}" ]; then
            SYLPH_TAX_CMD="conda run -p {params.env} sylph-tax"
        else
            SYLPH_TAX_CMD="conda run -n {params.env} sylph-tax"
        fi

        echo "=== Starting sylph-tax merge across samples ==="
        echo "Date: $(date)"
        echo "Output Relative Abundance: {output.rel_abun}"
        echo "Output Sequence Abundance: {output.seq_abun}"

        # Try merging all profiles; if sylph-tax merge fails due to empty profiles, fallback to non-empty profiles
        if ! $SYLPH_TAX_CMD merge {input.sylphmpas} --column relative_abundance -o "{output.rel_abun}" 2>/dev/null; then
            VALID_MPAS=()
            for mpa in {input.sylphmpas}; do
                if [ -f "$mpa" ] && [ $(wc -l < "$mpa") -gt 2 ]; then
                    VALID_MPAS+=("$mpa")
                fi
            done
            if [ ${{#VALID_MPAS[@]}} -gt 0 ]; then
                $SYLPH_TAX_CMD merge "${{VALID_MPAS[@]}}" --column relative_abundance -o "{output.rel_abun}"
            else
                printf "clade_name\n" > "{output.rel_abun}"
            fi
        fi

        if ! $SYLPH_TAX_CMD merge {input.sylphmpas} --column sequence_abundance -o "{output.seq_abun}" 2>/dev/null; then
            VALID_MPAS=()
            for mpa in {input.sylphmpas}; do
                if [ -f "$mpa" ] && [ $(wc -l < "$mpa") -gt 2 ]; then
                    VALID_MPAS+=("$mpa")
                fi
            done
            if [ ${{#VALID_MPAS[@]}} -gt 0 ]; then
                $SYLPH_TAX_CMD merge "${{VALID_MPAS[@]}}" --column sequence_abundance -o "{output.seq_abun}"
            else
                printf "clade_name\n" > "{output.seq_abun}"
            fi
        fi

        echo "=== Completed sylph-tax merge at $(date) ==="
        """


rule read_profiling_sylph_taxprof_batch:
    input:
        tsv = os.path.join(RP_RES, "profiling.{dataset}.tsv")
    output:
        sylphmpa = os.path.join(RP_RES, f"{RP_TAX_PREFIX}_{{dataset}}.sylphmpa")
    log:
        os.path.join(RP_LOGS, f"sylph_taxprof_batch/{{dataset}}.{RP_TAX_PREFIX}.log")
    threads: get_threads(RP_C.get("taxprof_threads", 8))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = rp_env("sylph"),
        taxonomy = rp_resolve_path(RP_C.get("taxonomy", "/dssg/home/acct-schhui/gses0/liliuyang/osma/software_from_nature/sylph/taxonomy_file_folder/gtdb_r226_metadata.tsv.gz")),
        prefix = lambda w: os.path.join(RP_RES, f"{RP_TAX_PREFIX}_{w.dataset}."),
        extra = RP_C.get("taxprof_params", "")
    wildcard_constraints:
        dataset = r"[^/]+"
    shell:
        """
        set -euo pipefail

        mkdir -p "$(dirname "{output.sylphmpa}")"
        mkdir -p "$(dirname "{log}")"

        exec > "{log}" 2>&1
        set -x

        ulimit -s unlimited || true
        ulimit -l unlimited || true

        for conda_sh in \
            /dssg/home/acct-schhui/gses0/Soft/Conda/etc/profile.d/conda.sh \
            "$HOME/miniconda3/etc/profile.d/conda.sh" \
            "$HOME/anaconda3/etc/profile.d/conda.sh"; do
            if [ -f "$conda_sh" ]; then
                source "$conda_sh"
                break
            fi
        done

        if [ -d "{params.env}/bin" ]; then
            export PATH="{params.env}/bin:$PATH"
        fi

        if command -v conda >/dev/null 2>&1; then
            conda activate "{params.env}" 2>/dev/null || conda activate sylph 2>/dev/null || true
        fi

        if [ -x "{params.env}/bin/python" ] && [ -f "{params.env}/bin/sylph-tax" ]; then
            SYLPH_TAX_CMD="{params.env}/bin/python {params.env}/bin/sylph-tax"
        elif [ -x "{params.env}/bin/sylph-tax" ]; then
            SYLPH_TAX_CMD="{params.env}/bin/sylph-tax"
        elif command -v sylph-tax >/dev/null 2>&1; then
            SYLPH_TAX_CMD="$(command -v sylph-tax)"
        elif [ -d "{params.env}" ]; then
            SYLPH_TAX_CMD="conda run -p {params.env} sylph-tax"
        else
            SYLPH_TAX_CMD="conda run -n {params.env} sylph-tax"
        fi

        echo "=== Starting sylph-tax taxprof batch: {wildcards.dataset} ==="
        echo "Date: $(date)"
        echo "Dataset: {wildcards.dataset}"
        echo "Input TSV: {input.tsv}"
        echo "Taxonomy: {params.taxonomy}"

        rm -f "$(dirname "{output.sylphmpa}")"/*{wildcards.dataset}*.sylphmpa

        $SYLPH_TAX_CMD taxprof "{input.tsv}" \
            -t "{params.taxonomy}" \
            -o "{params.prefix}" \
            {params.extra}

        if [ ! -f "{output.sylphmpa}" ]; then
            MATCH=$(find "$(dirname "{output.sylphmpa}")" -maxdepth 1 -name "*.sylphmpa" ! -name "$(basename "{output.sylphmpa}")" | head -n 1)
            if [ -n "$MATCH" ] && [ -f "$MATCH" ]; then
                mv -f "$MATCH" "{output.sylphmpa}"
            fi
        fi

        if [ ! -f "{output.sylphmpa}" ]; then
            echo "Notice: No taxonomy matches found in {input.tsv}. Creating empty sylphmpa profile."
            printf "#SampleID\t{wildcards.dataset}\tTaxonomies_used:{params.taxonomy}\nclade_name\trelative_abundance\tsequence_abundance\tANI (if strain-level)\tCoverage (if strain-level)\n" > "{output.sylphmpa}"
        fi

        echo "=== Completed sylph-tax taxprof batch at $(date) ==="
        """


# Explicit post-run inventory; never part of rule all or an automatic cleanup.
if os.path.basename(CONF_DIR) == "config_read" and SAMPLE_TSV_PATH:
    rule read_profiling_archive_plan:
        input:
            tsv = resolved_tsv,
            results = READ_PROFILING_TARGETS
        output:
            plan = os.path.join(RP_LOGS, "archive_plan.tsv")
        params:
            script = os.path.join(WORKFLOW_ROOT, "projects/config_read/scripts/archive_plan.py"),
            results_dir = RP_RES
        threads: 1
        shell:
            "python3 {params.script:q} --tsv {input.tsv:q} --results-dir {params.results_dir:q} --output {output.plan:q}"
