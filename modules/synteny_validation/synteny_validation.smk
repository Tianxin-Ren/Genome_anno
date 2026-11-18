# ==============================================================================
# Module: Synteny & Operon Gene-Loss Validation Engine (config_synteny)
# ==============================================================================
import os

SYNTENY_C = config.get("synteny_validation", config)

SYN_ROOT = config.get("project_root", "")

def syn_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    return os.path.normpath(os.path.join(SYN_ROOT, rel_path))

SYN_RES = os.path.join(
    SYN_ROOT,
    SYNTENY_C.get("results_dir", "results/synteny_validation"),
)
SYN_LOGS = os.path.join(
    SYN_ROOT,
    SYNTENY_C.get("logs_dir", "logs/synteny_validation"),
)

def syn_env(name):
    env_name = config.get("environments", {}).get(name) or SYNTENY_C.get("conda_envs", {}).get(name, name)
    env_dir = config.get("conda_env_dir", "") or config.get("paths", {}).get("conda_envs", "")
    if not env_dir or os.path.isabs(env_name):
        return env_name
    candidate = os.path.join(env_dir, env_name)
    if os.path.exists(candidate):
        return candidate
    if os.path.exists(env_dir):
        try:
            for entry in os.listdir(env_dir):
                if entry == env_name or entry.startswith(f"{env_name}_"):
                    return os.path.join(env_dir, entry)
        except Exception:
            pass
    return candidate

# ------------------------------------------------------------------------------
# 待测样本加载与生物域预设解析 (Sample Loading & Dynamic Domain Resolution)
# ------------------------------------------------------------------------------
domain = str(config.get("domain") or SYNTENY_C.get("domain", "bacteria")).strip().lower()
domain_presets = SYNTENY_C.get("domain_presets", {})
active_domain_cfg = domain_presets.get(domain, {})

# Sample loading (supports projects/config_synteny/samples.yaml, config/samples.yaml, or direct config injection)
SYNTENY_SAMPLES = config.get("synteny_samples", {})
if not SYNTENY_SAMPLES:
    SYNTENY_SAMPLES = config.get("external_genomes", {})
if not SYNTENY_SAMPLES:
    _samples_cands = [
        config.get("data", {}).get("samples", ""),
        syn_full_path("projects/config_synteny/samples.yaml"),
        "projects/config_synteny/samples.yaml",
        "config/samples.yaml"
    ]
    for _sc in _samples_cands:
        if _sc and os.path.exists(_sc):
            try:
                import yaml
                with open(_sc, "r", encoding="utf-8") as _scf:
                    _sc_data = yaml.safe_load(_scf) or {}
                    if "synteny_samples" in _sc_data:
                        SYNTENY_SAMPLES = _sc_data.get("synteny_samples", {})
                        for _k in ["target_genomes_dir", "target_proteins_dir", "target_genes_dir",
                                   "target_classification", "target_genomes_list", "target_annotation",
                                   "target_annotation_dir", "query_operon_faa", "default_domain", "default_presets"]:
                            if _k in _sc_data and _k not in config:
                                config[_k] = _sc_data[_k]
                        break
            except Exception:
                pass

if SYNTENY_SAMPLES is None:
    SYNTENY_SAMPLES = {}

def get_sample_domain(sample=None):
    if not sample:
        return domain
    cl_s = sample.replace("GB_", "").replace("RS_", "")
    info = SYNTENY_SAMPLES.get(cl_s, {}) or SYNTENY_SAMPLES.get(sample, {})
    if isinstance(info, dict):
        s_dom = info.get("domain")
        if s_dom and str(s_dom).strip():
            return str(s_dom).strip().lower()
        tax = str(info.get("taxonomy", "")).lower()
        if "d__archaea" in tax:
            return "archaea"
        if "d__bacteria" in tax:
            return "bacteria"
    if config.get("default_domain"):
        return str(config["default_domain"]).strip().lower()
    return domain

def get_domain_cfg(key, default="", sample=None):
    s_dom = get_sample_domain(sample) if sample else domain
    d_cfg = domain_presets.get(s_dom, active_domain_cfg)
    if sample:
        cl_s = sample.replace("GB_", "").replace("RS_", "")
        info = SYNTENY_SAMPLES.get(cl_s, {}) or SYNTENY_SAMPLES.get(sample, {})
        if isinstance(info, dict) and key in info and info[key] is not None and info[key] != "":
            return info[key]
    if key in config and config[key] is not None and config[key] != "":
        return config[key]
    if key in d_cfg and d_cfg[key] is not None:
        return d_cfg[key]
    if key in SYNTENY_C and SYNTENY_C[key] is not None:
        return SYNTENY_C[key]
    return default

gtdb_meta = syn_full_path(get_domain_cfg("gtdb_metadata", ""))
gtdb_genomes = syn_full_path(get_domain_cfg("gtdb_genomes_dir", ""))
gtdb_prots = syn_full_path(get_domain_cfg("gtdb_proteins_dir", ""))
target_genomes = syn_full_path(get_domain_cfg("target_genomes_dir", ""))
target_prots = syn_full_path(get_domain_cfg("target_proteins_dir", ""))
target_genes = syn_full_path(get_domain_cfg("target_genes_dir", ""))
target_class_file = syn_full_path(get_domain_cfg("target_classification", ""))

# Automatic batch ingestion from target_genomes_list if provided and SYNTENY_SAMPLES is empty
target_list_cfg = get_domain_cfg("target_genomes_list", "")
if not SYNTENY_SAMPLES and target_list_cfg and os.path.exists(syn_full_path(target_list_cfg)):
    target_ids = []
    with open(syn_full_path(target_list_cfg), "r", encoding="utf-8") as tlf:
        for t_line in tlf:
            t_line = t_line.strip()
            if not t_line or t_line.startswith("#"):
                continue
            tok = t_line.split()[0].split("\t")[0].strip()
            if tok and tok.lower() not in ["genome_id", "accession", "id", "mag_id", "mag", "sample"]:
                target_ids.append(tok)

    meta_tax_map = {}

    # 1. Ingest target classification first (e.g. NS_MOB mag_classification)
    if target_class_file and os.path.exists(target_class_file):
        try:
            with open(target_class_file, "r", encoding="utf-8") as _tcf:
                for _tc_line in _tcf:
                    _tc_line = _tc_line.strip()
                    if not _tc_line or _tc_line.startswith("#"):
                        continue
                    _parts = _tc_line.split("\t")
                    if len(_parts) >= 2 and _parts[0].lower() not in ["mag_id", "genome_id", "accession", "id", "mag", "sample"]:
                        _gid = _parts[0].strip()
                        _gtax = _parts[1].strip()
                        _rec = {"taxonomy": _gtax, "completeness": 100.0, "contamination": 0.0}
                        meta_tax_map[_gid] = _rec
                        _cl_g = _gid.replace("GB_", "").replace("RS_", "")
                        meta_tax_map[_cl_g] = _rec
                        if "_" in _cl_g:
                            _p01 = "_".join(_cl_g.split("_")[:2])
                            meta_tax_map[_p01] = _rec
        except Exception:
            pass

    # 2. Check if any target_ids are still missing taxonomy
    missing_targets = set()
    for tid in target_ids:
        cl_id = tid.replace("GB_", "").replace("RS_", "")
        if tid in meta_tax_map or cl_id in meta_tax_map:
            continue
        if "_" in cl_id and "_".join(cl_id.split("_")[:2]) in meta_tax_map:
            continue
        missing_targets.add(tid)
        missing_targets.add(cl_id)
        if "_" in cl_id:
            missing_targets.add("_".join(cl_id.split("_")[:2]))

    # 3. Stream GTDB metadata ONLY for missing targets (zero-memory streaming, avoids OOM on large metadata)
    if missing_targets and gtdb_meta and os.path.exists(gtdb_meta):
        try:
            with open(gtdb_meta, "r", encoding="utf-8", errors="ignore") as _gmf:
                _header = _gmf.readline()
                if _header:
                    _cols = [c.strip() for c in _header.split("\t")]
                    _ac_idx = _cols.index("accession") if "accession" in _cols else 0
                    _tc_idx = _cols.index("gtdb_taxonomy") if "gtdb_taxonomy" in _cols else (_cols.index("taxonomy") if "taxonomy" in _cols else 1)
                    _comp_idx = _cols.index("checkm_completeness") if "checkm_completeness" in _cols else (_cols.index("completeness") if "completeness" in _cols else -1)
                    _contam_idx = _cols.index("checkm_contamination") if "checkm_contamination" in _cols else (_cols.index("contamination") if "contamination" in _cols else -1)

                    _rem = set(missing_targets)
                    for _line in _gmf:
                        if not _rem:
                            break
                        if not any(_k in _line for _k in _rem):
                            continue
                        _parts = _line.split("\t")
                        if len(_parts) > _ac_idx:
                            _acc = _parts[_ac_idx].strip()
                            _acc_cl = _acc.replace("GB_", "").replace("RS_", "")
                            _acc_p01 = "_".join(_acc_cl.split("_")[:2]) if "_" in _acc_cl else _acc_cl

                            if _acc in _rem or _acc_cl in _rem or _acc_p01 in _rem:
                                _gtax = _parts[_tc_idx].strip() if len(_parts) > _tc_idx else "d__Bacteria;p__Unknown;c__Unknown;o__Unknown;f__Unknown;g__Unknown;s__"
                                _comp = float(_parts[_comp_idx]) if (_comp_idx >= 0 and len(_parts) > _comp_idx and _parts[_comp_idx].strip()) else 100.0
                                _contam = float(_parts[_contam_idx]) if (_contam_idx >= 0 and len(_parts) > _contam_idx and _parts[_contam_idx].strip()) else 0.0
                                _rec = {"taxonomy": _gtax, "completeness": _comp, "contamination": _contam}
                                meta_tax_map[_acc] = _rec
                                meta_tax_map[_acc_cl] = _rec
                                meta_tax_map[_acc_p01] = _rec
                                _rem.discard(_acc)
                                _rem.discard(_acc_cl)
                                _rem.discard(_acc_p01)
        except Exception:
            pass

    # Fast index GTDB genome and protein files for instant O(1) matching with assembly suffixes
    import glob
    import re

    _in_memory_index_cache = {}

    def _build_gtdb_file_index(directory, valid_exts):
        idx = {}
        if not directory or not os.path.isdir(directory):
            return idx
        cache_key = (directory, tuple(valid_exts))
        if cache_key in _in_memory_index_cache:
            return _in_memory_index_cache[cache_key]
        try:
            visited_dirs = set()
            for root, dirs, files in os.walk(directory, followlinks=False):
                real_r = os.path.realpath(root)
                if real_r in visited_dirs:
                    dirs.clear()
                    continue
                visited_dirs.add(real_r)
                for fname in files:
                    matched_ext = None
                    for ext in valid_exts:
                        if fname.endswith(ext):
                            matched_ext = ext
                            break
                    if matched_ext is not None:
                        full_p = os.path.join(root, fname)
                        idx[fname] = full_p
                        base_name = fname[:-len(matched_ext)]
                        idx[base_name] = full_p
                        idx[base_name.replace("GB_", "").replace("RS_", "")] = full_p

                        m = re.search(r"(GC[AF])_(\d+\.\d+)", fname)
                        if m:
                            prefix, num = m.group(1), m.group(2)
                            acc_key = f"{prefix}_{num}"
                            idx[acc_key] = full_p
                            idx[f"GB_{acc_key}"] = full_p
                            idx[f"RS_{acc_key}"] = full_p
                            idx[num] = full_p

                        m2 = re.search(r"(GC[AF]_\d+\.\d+_[^.]+)", fname)
                        if m2:
                            full_acc = m2.group(1)
                            idx[full_acc] = full_p
                            idx[f"GB_{full_acc}"] = full_p
                            idx[f"RS_{full_acc}"] = full_p
            _in_memory_index_cache[cache_key] = idx
        except Exception:
            pass
        return idx

    target_genomes = syn_full_path(get_domain_cfg("target_genomes_dir", ""))
    target_prots = syn_full_path(get_domain_cfg("target_proteins_dir", ""))
    target_genes = syn_full_path(get_domain_cfg("target_genes_dir", ""))

    valid_genome_exts = (".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")
    genome_file_index = {}
    if target_genomes and os.path.isdir(target_genomes):
        genome_file_index.update(_build_gtdb_file_index(target_genomes, valid_genome_exts))
    if gtdb_genomes and os.path.isdir(gtdb_genomes):
        genome_file_index.update(_build_gtdb_file_index(gtdb_genomes, valid_genome_exts))

    valid_protein_exts = (".faa", ".faa.gz")
    search_prot_dirs = []
    if target_prots and os.path.isdir(target_prots):
        search_prot_dirs.append(target_prots)
        search_prot_dirs.append(os.path.join(target_prots, "all"))
    if gtdb_prots and os.path.isdir(gtdb_prots):
        search_prot_dirs.append(gtdb_prots)
        all_sub = os.path.join(gtdb_prots, "all")
        if os.path.isdir(all_sub):
            search_prot_dirs.append(all_sub)

    protein_file_index = {}
    for pd_dir in set(search_prot_dirs):
        if os.path.isdir(pd_dir):
            protein_file_index.update(_build_gtdb_file_index(pd_dir, valid_protein_exts))

    valid_gene_exts = (".gff", ".gff.gz", ".gff3", ".gff3.gz")
    gene_file_index = {}
    if target_genes and os.path.isdir(target_genes):
        gene_file_index.update(_build_gtdb_file_index(target_genes, valid_gene_exts))

    for tid in target_ids:
        cl_id = tid.replace("GB_", "").replace("RS_", "")
        acc_short = "_".join(cl_id.split("_")[:2]) if "_" in cl_id else cl_id

        m_inf = meta_tax_map.get(tid) or meta_tax_map.get(cl_id) or meta_tax_map.get(acc_short)
        if not m_inf:
            m_inf = {"taxonomy": "d__Bacteria;p__Unknown;c__Unknown;o__Unknown;f__Unknown;g__Unknown;s__", "completeness": 100.0, "contamination": 0.0}

        # 1. Genome FASTA resolution
        fa_cand = (
            genome_file_index.get(cl_id) or 
            genome_file_index.get(tid) or 
            genome_file_index.get(acc_short) or
            genome_file_index.get(cl_id.replace("GB_", "").replace("RS_", ""))
        )
        if not fa_cand and target_genomes and os.path.isdir(target_genomes):
            matches = (
                glob.glob(os.path.join(target_genomes, f"*{cl_id}*")) or
                glob.glob(os.path.join(target_genomes, f"**/*{cl_id}*"), recursive=True) or
                glob.glob(os.path.join(target_genomes, f"*{acc_short}*")) or
                glob.glob(os.path.join(target_genomes, f"**/*{acc_short}*"), recursive=True)
            )
            for mf in matches:
                if any(mf.endswith(ext) for ext in valid_genome_exts):
                    fa_cand = mf
                    break

        if not fa_cand and gtdb_genomes and os.path.isdir(gtdb_genomes):
            matches = (
                glob.glob(os.path.join(gtdb_genomes, f"*{cl_id}*")) or
                glob.glob(os.path.join(gtdb_genomes, f"*{acc_short}*"))
            )
            for mf in matches:
                if any(mf.endswith(ext) for ext in valid_genome_exts):
                    fa_cand = mf
                    break

        # Defense check: if genome FASTA is not physically present, skip to avoid MissingInputException
        if not fa_cand or not os.path.exists(fa_cand):
            print(f"[Synteny Warning] Genome FASTA for sample '{cl_id}' ({tid}) not found on disk in target_genomes_dir or gtdb_genomes_dir. Skipping sample from run.", file=sys.stderr)
            continue

        # 2. Proteome FAA resolution
        faa_cand = (
            protein_file_index.get(cl_id) or 
            protein_file_index.get(tid) or 
            protein_file_index.get(acc_short)
        )
        if not faa_cand and target_prots and os.path.isdir(target_prots):
            matches = (
                glob.glob(os.path.join(target_prots, f"*{cl_id}*.faa*")) or 
                glob.glob(os.path.join(target_prots, f"**/*{cl_id}*.faa*"), recursive=True) or
                glob.glob(os.path.join(target_prots, f"*{acc_short}*.faa*")) or 
                glob.glob(os.path.join(target_prots, f"**/*{acc_short}*.faa*"), recursive=True)
            )
            if matches:
                faa_cand = matches[0]
        if not faa_cand and gtdb_prots and os.path.isdir(gtdb_prots):
            matches = (
                glob.glob(os.path.join(gtdb_prots, f"*{cl_id}*.faa*")) or 
                glob.glob(os.path.join(gtdb_prots, f"**/*{cl_id}*.faa*"), recursive=True) or
                glob.glob(os.path.join(gtdb_prots, f"*{acc_short}*.faa*")) or 
                glob.glob(os.path.join(gtdb_prots, f"**/*{acc_short}*.faa*"), recursive=True)
            )
            if matches:
                faa_cand = matches[0]

        # 3. Gene GFF resolution
        gff_cand = (
            gene_file_index.get(cl_id) or 
            gene_file_index.get(tid) or 
            gene_file_index.get(acc_short)
        )
        if not gff_cand and target_genes and os.path.isdir(target_genes):
            matches = (
                glob.glob(os.path.join(target_genes, f"*{cl_id}*.gff*")) or 
                glob.glob(os.path.join(target_genes, f"**/*{cl_id}*.gff*"), recursive=True) or
                glob.glob(os.path.join(target_genes, f"*{acc_short}*.gff*")) or 
                glob.glob(os.path.join(target_genes, f"**/*{acc_short}*.gff*"), recursive=True)
            )
            if matches:
                gff_cand = matches[0]

        sample_entry = {
            "fasta": fa_cand or "",
            "faa": faa_cand or "",
            "gff": gff_cand or "",
            "taxonomy": m_inf["taxonomy"],
            "completeness": m_inf.get("completeness", 100.0),
            "contamination": m_inf.get("contamination", 0.0)
        }
        SYNTENY_SAMPLES[cl_id] = sample_entry

def get_syn_input(sample, key, default=""):
    cl_s = sample.replace("GB_", "").replace("RS_", "") if sample else ""
    info = SYNTENY_SAMPLES.get(cl_s, {}) or SYNTENY_SAMPLES.get(sample, {})
    if isinstance(info, dict):
        return info.get(key, default)
    return default

def get_target_fa(wildcards):
    s = wildcards.sample
    path = get_syn_input(s, "fasta") or get_syn_input(s, "path")
    if path and os.path.exists(syn_full_path(path)):
        return syn_full_path(path)
    cl_s = s.replace("GB_", "").replace("RS_", "")
    search_dirs = []
    if target_genomes:
        search_dirs.extend([target_genomes, os.path.join(target_genomes, "all")])
    if gtdb_genomes:
        search_dirs.extend([gtdb_genomes, os.path.join(gtdb_genomes, "all")])
    for d in set(search_dirs):
        if d and os.path.isdir(d):
            m = glob.glob(os.path.join(d, f"*{cl_s}*.fna*")) or glob.glob(os.path.join(d, f"*{cl_s}*.fa*")) or glob.glob(os.path.join(d, f"*{cl_s}*.fasta*"))
            if m:
                return m[0]
    return syn_full_path(path)

def get_target_faa(wildcards):
    s = wildcards.sample
    path = get_syn_input(s, "faa")
    if path and os.path.exists(syn_full_path(path)):
        return syn_full_path(path)
    cl_s = s.replace("GB_", "").replace("RS_", "")
    cand_dirs = []
    if target_prots:
        cand_dirs.extend([
            target_prots,
            os.path.join(target_prots, "all")
        ])
    if gtdb_prots:
        cand_dirs.extend([
            gtdb_prots,
            os.path.join(gtdb_prots, "all")
        ])
    for cd in set(cand_dirs):
        if cd and os.path.isdir(cd):
            m = glob.glob(os.path.join(cd, f"*{cl_s}*.faa*")) or glob.glob(os.path.join(cd, f"*{s}*.faa*"))
            if m:
                return m[0]
    fa = get_syn_input(s, "fasta") or get_syn_input(s, "path")
    if fa:
        sibling = fa.replace(".fna", ".faa").replace(".fasta", ".faa").replace(".fa", ".faa")
        if sibling and os.path.exists(syn_full_path(sibling)):
            return syn_full_path(sibling)
    return os.path.join(SYN_RES, f"{s}/predicted/{s}.faa")

# ------------------------------------------------------------------------------
# 操纵子预设与多轨模块动态解析引擎 (Operon Preset & Module Resolution Engine)
# ------------------------------------------------------------------------------
operon_presets = SYNTENY_C.get("operon_presets", {})

MODULE_DEFS = {}
# First register all module definitions across all known presets so any sample can look up module specs
for p_name, p_data in operon_presets.items():
    if isinstance(p_data, dict) and "modules" in p_data:
        MODULE_DEFS.update(p_data["modules"])
if "operon_modules" in SYNTENY_C and isinstance(SYNTENY_C["operon_modules"], dict):
    MODULE_DEFS.update(SYNTENY_C["operon_modules"])
if "atp_synthase_modules" in SYNTENY_C and isinstance(SYNTENY_C["atp_synthase_modules"], dict):
    MODULE_DEFS.update(SYNTENY_C["atp_synthase_modules"])

def parse_preset_list(val):
    if isinstance(val, list):
        return [str(p).strip() for p in val if str(p).strip()]
    elif isinstance(val, str):
        if "," in val:
            return [p.strip() for p in val.split(",") if p.strip()]
        elif val.strip():
            return [val.strip()]
    return []

# Support both active_presets (list/str) and active_preset (scalar/list)
raw_active = (
    config.get("active_presets")
    or config.get("active_preset")
    or config.get("default_presets")
    or config.get("default_preset")
    or (active_domain_cfg.get("default_preset") if "domain" in config else None)
    or SYNTENY_C.get("active_presets")
    or SYNTENY_C.get("active_preset")
    or active_domain_cfg.get("default_preset", "sox_system" if domain == "bacteria" else "M00159")
)

active_preset_names = parse_preset_list(raw_active)

def get_sample_presets(sample):
    if not sample:
        return active_preset_names
    cl_s = sample.replace("GB_", "").replace("RS_", "")
    info = SYNTENY_SAMPLES.get(cl_s, {}) or SYNTENY_SAMPLES.get(sample, {})
    if isinstance(info, dict):
        raw_p = info.get("presets") or info.get("preset") or info.get("active_presets") or info.get("active_preset")
        if raw_p:
            parsed = parse_preset_list(raw_p)
            if parsed:
                return parsed
    # Fallback to default_presets from samples.yaml or config
    if config.get("default_presets"):
        parsed = parse_preset_list(config["default_presets"])
        if parsed:
            return parsed
    # Fallback to global active_presets
    if active_preset_names:
        return active_preset_names
    # Fallback to sample's domain default
    s_dom = get_sample_domain(sample)
    d_default = domain_presets.get(s_dom, {}).get("default_preset", "sox_system" if s_dom == "bacteria" else "M00159")
    return parse_preset_list(d_default)

def get_sample_modules(sample):
    # Explicit CLI override from config takes highest precedence if provided
    user_enabled = config.get("enabled_modules") or SYNTENY_C.get("enabled_modules")
    if user_enabled:
        return parse_preset_list(user_enabled)

    p_names = get_sample_presets(sample)
    mods = []
    for p in p_names:
        if p in operon_presets:
            p_data = operon_presets[p]
            p_mods = p_data.get("enabled_modules", list(p_data.get("modules", {}).keys()))
            for m in p_mods:
                if m not in mods:
                    mods.append(m)
        elif p in MODULE_DEFS:
            if p not in mods:
                mods.append(p)
        elif p in ["atp", "atp_synthase"]:
            for m in ["M00159", "M00157", "M00158", "M00160"]:
                if m not in mods:
                    mods.append(m)
        elif p in ["sox", "sox_system"]:
            if "M00595" not in mods:
                mods.append("M00595")
        elif re.match(r"^M\d{5}$", p):
            if p not in mods:
                mods.append(p)
    return mods if mods else ["M00159"]

# Aggregate all enabled modules across all samples for global fallbacks
ENABLED_MODULES = []
for _s in sorted(SYNTENY_SAMPLES.keys()):
    for _m in get_sample_modules(_s):
        if _m not in ENABLED_MODULES:
            ENABLED_MODULES.append(_m)

if not ENABLED_MODULES:
    for p in active_preset_names:
        if p in operon_presets:
            for m in operon_presets[p].get("enabled_modules", []):
                if m not in ENABLED_MODULES:
                    ENABLED_MODULES.append(m)
        elif p in MODULE_DEFS and p not in ENABLED_MODULES:
            ENABLED_MODULES.append(p)

if not ENABLED_MODULES:
    ENABLED_MODULES = ["M00159"]

def get_module_kos(wildcards):
    mod = getattr(wildcards, "module", ENABLED_MODULES[0] if ENABLED_MODULES else "M00595")
    if mod in MODULE_DEFS and "target_kos" in MODULE_DEFS[mod]:
        kos = MODULE_DEFS[mod]["target_kos"]
        if isinstance(kos, list):
            return ",".join(kos)
        return str(kos)
    kos = SYNTENY_C.get("target_kos", [])
    if isinstance(kos, list):
        return ",".join(kos)
    return str(kos)

def get_module_alt_kos(wildcards):
    mod = getattr(wildcards, "module", ENABLED_MODULES[0] if ENABLED_MODULES else "M00595")
    if mod in MODULE_DEFS:
        if "alternative_kos" in MODULE_DEFS[mod]:
            kos = MODULE_DEFS[mod]["alternative_kos"]
            return ",".join(kos) if isinstance(kos, list) else str(kos)
        alt_to = MODULE_DEFS[mod].get("alternative_to", "")
        if alt_to and alt_to in MODULE_DEFS and "target_kos" in MODULE_DEFS[alt_to]:
            kos = MODULE_DEFS[alt_to]["target_kos"]
            return ",".join(kos) if isinstance(kos, list) else str(kos)
    return ""

def get_module_name(wildcards):
    mod = getattr(wildcards, "module", ENABLED_MODULES[0] if ENABLED_MODULES else "M00595")
    if mod in MODULE_DEFS and "name" in MODULE_DEFS[mod]:
        return MODULE_DEFS[mod]["name"]
    return mod

def get_module_expected_len(wildcards):
    mod = getattr(wildcards, "module", ENABLED_MODULES[0] if ENABLED_MODULES else "M00595")
    if mod in MODULE_DEFS and "expected_locus_length" in MODULE_DEFS[mod]:
        return int(MODULE_DEFS[mod]["expected_locus_length"])
    return int(SYNTENY_C.get("expected_locus_length", 10000))

def get_query_operon_input(wildcards):
    mod = getattr(wildcards, "module", "")
    q = ""
    if mod and mod in MODULE_DEFS:
        q = MODULE_DEFS[mod].get("query_operon_faa", "")
    if not q:
        q = get_domain_cfg("query_operon_faa", "")
    if q and os.path.exists(syn_full_path(q)):
        return syn_full_path(q)
    return []

def get_target_gff(wildcards):
    s = wildcards.sample
    path = get_syn_input(s, "gff") or get_syn_input(s, "gff3")
    if path and os.path.exists(syn_full_path(path)):
        return syn_full_path(path)
    cl_s = s.replace("GB_", "").replace("RS_", "")
    target_genes = syn_full_path(SYNTENY_C.get("target_genes_dir", ""))
    cand_roots = []
    if target_genes:
        cand_roots.append(target_genes)
    if gtdb_prots:
        cand_roots.extend([gtdb_prots, os.path.dirname(gtdb_prots.rstrip("/")), os.path.dirname(os.path.dirname(gtdb_prots.rstrip("/")))])
    if gtdb_genomes:
        cand_roots.extend([gtdb_genomes, os.path.dirname(gtdb_genomes.rstrip("/"))])
    for cr in set(cand_roots):
        if cr and os.path.isdir(cr):
            for sub in ["genes", "gff", "metadata"]:
                sub_dir = os.path.join(cr, sub)
                if os.path.isdir(sub_dir):
                    m = glob.glob(os.path.join(sub_dir, f"*{cl_s}*.gff*")) or glob.glob(os.path.join(sub_dir, f"*{s}*.gff*"))
                    if m:
                        return m[0]
    return ""

def get_target_bam(wildcards):
    s = wildcards.sample
    path = get_syn_input(s, "reads_bam") or get_syn_input(s, "bam")
    return syn_full_path(path) if path else ""

def get_target_annotation(wildcards):
    s = wildcards.sample
    # 1. From sample-specific config: annotation or kofam_tsv
    anno = get_syn_input(s, "annotation") or get_syn_input(s, "kofam_tsv")
    if anno and os.path.exists(syn_full_path(anno)):
        return syn_full_path(anno)

    # 2. From global template in synteny_validation.yaml (e.g. "functional/{sample}.kofam.tsv")
    tmpl = SYNTENY_C.get("target_annotation", "")
    if tmpl:
        formatted = tmpl.format(sample=s)
        if os.path.exists(syn_full_path(formatted)):
            return syn_full_path(formatted)
        if "_" in s:
            suf = s.split("_")[-1]
            formatted_suf = tmpl.format(sample=suf)
            if os.path.exists(syn_full_path(formatted_suf)):
                return syn_full_path(formatted_suf)

    # 3. From global directory: functional/{sample}.kofam.tsv or functional/{suf}.kofam.tsv
    anno_dir = SYNTENY_C.get("target_annotation_dir", "functional")
    for cand_name in [s, s.split("_")[-1] if "_" in s else s]:
        for ext in [".kofam.tsv", ".tsv"]:
            cand = os.path.join(anno_dir, f"{cand_name}{ext}")
            if os.path.exists(syn_full_path(cand)):
                return syn_full_path(cand)

    return ""

def get_target_annotation_input(wildcards):
    path = get_target_annotation(wildcards)
    return [path] if path and os.path.exists(path) else []

ENABLE_CLINKER = bool(SYNTENY_C.get("run_clinker", SYNTENY_C.get("enable_clinker", False)))
SYNTENY_TARGETS = []

for cl_id in sorted(SYNTENY_SAMPLES.keys()):
    for mod in get_sample_modules(cl_id):
        SYNTENY_TARGETS.append(os.path.join(SYN_RES, f"{cl_id}/{mod}/stage4/{cl_id}_locus_validation_summary.tsv"))
        SYNTENY_TARGETS.append(os.path.join(SYN_RES, f"{cl_id}/{mod}/stage1/{cl_id}_100kb_locus_proteins.tsv"))
        if ENABLE_CLINKER:
            SYNTENY_TARGETS.append(os.path.join(SYN_RES, f"{cl_id}/{mod}/visualization/{cl_id}_clinker.html"))

if SYNTENY_TARGETS:
    SYNTENY_TARGETS.append(os.path.join(SYN_RES, "locus_synteny_all_modules_summary.tsv"))
    SYNTENY_TARGETS.append(os.path.join(SYN_RES, "atp_synthase_synteny_matrix.tsv"))
    SYNTENY_TARGETS.append(os.path.join(SYN_RES, "locus_synteny_validation_summary.tsv"))

# ------------------------------------------------------------------------------
# Rules
# ------------------------------------------------------------------------------

rule synteny_predict_target_faa:
    input:
        fa = get_target_fa
    output:
        faa = os.path.join(SYN_RES, "{sample}/predicted/{sample}.faa")
    log:
        os.path.join(SYN_LOGS, "prodigal/{sample}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = syn_env("diamond")
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p $(dirname "{output.faa}")
        cmd_prefix=""
        if [ -n "{params.env}" ]; then
            if [ -d "{params.env}" ] || [[ "{params.env}" == *"/"* ]]; then
                cmd_prefix="conda run --no-capture-output -p {params.env} "
            else
                cmd_prefix="conda run --no-capture-output -n {params.env} "
            fi
        fi
        $cmd_prefix prodigal -i "{input.fa}" -a "{output.faa}" -o /dev/null -p single -q
        }} 2>&1 | tee {log}
        """

rule synteny_stage0_select_refs:
    input:
        target_fa = get_target_fa,
        target_faa = get_target_faa,
        query_operon = get_query_operon_input
    output:
        refs_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage0/selected_validated_references.tsv"),
        anchors_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage0/consensus_anchors.tsv"),
        seed_query_faa = os.path.join(SYN_RES, "{sample}/{module}/stage0/{sample}_seed_query_operon.faa"),
        stage0_dir = directory(os.path.join(SYN_RES, "{sample}/{module}/stage0/ref_data"))
    log:
        os.path.join(SYN_LOGS, "stage0/{sample}_{module}.log")
    threads: get_threads(8)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        domain = lambda wildcards: get_sample_domain(wildcards.sample),
        gtdb_metadata = lambda wildcards: syn_full_path(get_domain_cfg("gtdb_metadata", "", wildcards.sample)),
        gtdb_genomes_dir = lambda wildcards: syn_full_path(get_domain_cfg("gtdb_genomes_dir", "", wildcards.sample)),
        taxonomy = lambda wildcards: get_syn_input(wildcards.sample, "taxonomy", "d__Bacteria;" if get_sample_domain(wildcards.sample) == "bacteria" else "d__Archaea;"),
        completeness = lambda wildcards: float(get_syn_input(wildcards.sample, "completeness", 100.0)),
        contamination = lambda wildcards: float(get_syn_input(wildcards.sample, "contamination", 0.0)),
        min_af = lambda wildcards: float(get_domain_cfg("min_af", SYNTENY_C.get("min_af", 0.50), wildcards.sample)),
        adaptive_af = lambda wildcards: "--adaptive_af" if get_domain_cfg("adaptive_af_for_mags", SYNTENY_C.get("adaptive_af_for_mags", True), wildcards.sample) else "",
        min_ani = lambda wildcards: float(get_domain_cfg("min_ani", SYNTENY_C.get("min_ani", 80.0), wildcards.sample)),
        max_pos = lambda wildcards: int(get_domain_cfg("max_positive_refs", 10 if get_sample_domain(wildcards.sample) == "bacteria" else 5, wildcards.sample)),
        max_neg = lambda wildcards: int(get_domain_cfg("max_negative_refs", 1 if get_sample_domain(wildcards.sample) == "bacteria" else 2, wildcards.sample)),
        skani_env = lambda wildcards: syn_env("skani"),
        diamond_env = lambda wildcards: syn_env("diamond"),
        py_env = lambda wildcards: syn_env("biopython"),
        module_id = lambda wildcards: wildcards.module,
        module_name = lambda wildcards: get_module_name(wildcards),
        target_kos = lambda wildcards: get_module_kos(wildcards),
        gtdb_ko_table = lambda wildcards: syn_full_path(get_domain_cfg("gtdb_ko_table", "", wildcards.sample)),
        ko_index_cache = lambda wildcards: os.path.join(SYN_RES, f"cache/{wildcards.module}_gtdb_index.tsv"),
        file_index_cache = lambda wildcards: os.path.join(SYN_RES, "cache/gtdb_file_index.json"),
        gtdb_proteins_dir = lambda wildcards: syn_full_path(get_domain_cfg("gtdb_proteins_dir", "", wildcards.sample)),
        query_operon_faa = lambda wildcards: syn_full_path(MODULE_DEFS.get(wildcards.module, {}).get("query_operon_faa", get_domain_cfg("query_operon_faa", "", wildcards.sample))),
        only_representative = lambda wildcards: "--only_representative" if get_domain_cfg("only_representative", True if get_sample_domain(wildcards.sample) == "bacteria" else False, wildcards.sample) else "",
        auto_fetch_flag = lambda wildcards: "--auto_fetch" if get_domain_cfg("auto_fetch", True if get_sample_domain(wildcards.sample) == "bacteria" else False, wildcards.sample) else "--no_auto_fetch"
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p "{output.stage0_dir}"
        
        target_kos_arg=""
        if [ -n "{params.target_kos}" ]; then
            target_kos_arg="--target_kos {params.target_kos}"
        fi
        gtdb_ko_arg=""
        if [ -n "{params.gtdb_ko_table}" ]; then
            gtdb_ko_arg="--gtdb_ko_table {params.gtdb_ko_table}"
        fi
        cache_arg=""
        if [ -n "{params.ko_index_cache}" ]; then
            cache_arg="--ko_index_cache {params.ko_index_cache}"
        fi
        file_cache_arg=""
        if [ -n "{params.file_index_cache}" ]; then
            file_cache_arg="--file_index_cache {params.file_index_cache}"
        fi
        gtdb_proteins_arg=""
        if [ -n "{params.gtdb_proteins_dir}" ]; then
            gtdb_proteins_arg="--gtdb_proteins_dir {params.gtdb_proteins_dir}"
        fi
        query_arg=""
        if [ -n "{params.query_operon_faa}" ]; then
            query_arg="--query_operon_faa {params.query_operon_faa}"
        fi

        py_cmd=""
        if [ -n "{params.py_env}" ]; then
            if [ -d "{params.py_env}" ] || [[ "{params.py_env}" == *"/"* ]]; then
                py_cmd="conda run --no-capture-output -p {params.py_env} "
            else
                py_cmd="conda run --no-capture-output -n {params.py_env} "
            fi
        fi

        $py_cmd python3 scripts/select_and_validate_refs.py \
            --domain "{params.domain}" \
            {params.auto_fetch_flag} \
            --target_name "{wildcards.sample}" \
            --target_fa "{input.target_fa}" \
            --target_faa "{input.target_faa}" \
            --target_taxonomy "{params.taxonomy}" \
            --target_completeness {params.completeness} \
            --target_contamination {params.contamination} \
            --module_id "{params.module_id}" \
            --module_name "{params.module_name}" \
            --gtdb_metadata "{params.gtdb_metadata}" \
            --gtdb_genomes_dir "{params.gtdb_genomes_dir}" \
            $query_arg \
            --min_af {params.min_af} \
            {params.adaptive_af} \
            --min_ani {params.min_ani} \
            --max_positive_refs {params.max_pos} \
            --max_negative_refs {params.max_neg} \
            --skani_env "{params.skani_env}" \
            --diamond_env "{params.diamond_env}" \
            --threads {threads} \
            $target_kos_arg \
            $gtdb_ko_arg \
            $cache_arg \
            $file_cache_arg \
            $gtdb_proteins_arg \
            {params.only_representative} \
            --out_refs_tsv "{output.refs_tsv}" \
            --out_anchors_tsv "{output.anchors_tsv}" \
            --out_seed_operon_faa "{output.seed_query_faa}" \
            --out_dir "{output.stage0_dir}"
        }} 2>&1 | tee {log}
        """

rule synteny_stage1_extract_anchor_loci:
    input:
        target_fa = get_target_fa,
        target_faa = get_target_faa,
        refs_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage0/selected_validated_references.tsv"),
        anchors_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage0/consensus_anchors.tsv"),
        seed_operon_faa = os.path.join(SYN_RES, "{sample}/{module}/stage0/{sample}_seed_query_operon.faa"),
        target_annotation = get_target_annotation_input
    output:
        topology_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage1/anchor_topology.tsv"),
        target_locus_fa = os.path.join(SYN_RES, "{sample}/{module}/stage1/{sample}_locus.fa"),
        target_locus_gbk = os.path.join(SYN_RES, "{sample}/{module}/stage1/{sample}_locus.gbk"),
        extended_proteins_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage1/{sample}_100kb_locus_proteins.tsv"),
        refs_locus_dir = directory(os.path.join(SYN_RES, "{sample}/{module}/stage1/refs_loci"))
    log:
        os.path.join(SYN_LOGS, "stage1/{sample}_{module}.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        target_gff = lambda wildcards: get_target_gff(wildcards),
        target_annotation = lambda wildcards: get_target_annotation(wildcards),
        module_id = lambda wildcards: wildcards.module,
        module_name = lambda wildcards: get_module_name(wildcards),
        expected_len = lambda wildcards: get_module_expected_len(wildcards),
        min_edge = lambda wildcards: int(SYNTENY_C.get("min_edge_bp", 10000)),
        flank_bp = lambda wildcards: int(SYNTENY_C.get("flank_bp", 20000)),
        n_flank = lambda wildcards: int(SYNTENY_C.get("n_flank_genes", SYNTENY_C.get("flank_genes", 25))),
        max_vis_refs = lambda wildcards: int(SYNTENY_C.get("max_visualized_refs", 10)),
        max_candidate_loci = lambda wildcards: int(SYNTENY_C.get("operon_fate_classifier", {}).get("max_candidate_loci", 5)),
        consensus_contig_threshold = lambda wildcards: float(SYNTENY_C.get("consensus_contig_threshold", 0.75)),
        significance_mode = lambda wildcards: SYNTENY_C.get("significance_mode", "empirical_max"),
        significance_alpha = lambda wildcards: float(SYNTENY_C.get("significance_alpha", 0.01)),
        n_permutations = lambda wildcards: int(SYNTENY_C.get("n_permutations", 1000)),
        extended_window_bp = lambda wildcards: int(SYNTENY_C.get("extended_window_bp", 100000)),
        target_kos = lambda wildcards: get_module_kos(wildcards),
        gtdb_ko_table = lambda wildcards: syn_full_path(get_domain_cfg("gtdb_ko_table", "", wildcards.sample)),
        ko_index_cache = lambda wildcards: os.path.join(SYN_RES, f"cache/{wildcards.module}_gtdb_index.tsv"),
        diamond_env = lambda wildcards: syn_env("diamond"),
        py_env = lambda wildcards: syn_env("biopython")
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p "{output.refs_locus_dir}"

        anno_arg=""
        if [ -n "{params.target_annotation}" ] && [ -f "{params.target_annotation}" ]; then
            anno_arg="--target_annotation_tsv {params.target_annotation}"
        fi

        gtdb_ko_arg=""
        if [ -n "{params.gtdb_ko_table}" ] && [ -f "{params.gtdb_ko_table}" ]; then
            gtdb_ko_arg="--gtdb_ko_table {params.gtdb_ko_table}"
        fi

        cache_arg=""
        if [ -n "{params.ko_index_cache}" ]; then
            cache_arg="--ko_index_cache {params.ko_index_cache}"
        fi

        py_cmd=""
        if [ -n "{params.py_env}" ]; then
            if [ -d "{params.py_env}" ] || [[ "{params.py_env}" == *"/"* ]]; then
                py_cmd="conda run --no-capture-output -p {params.py_env} "
            else
                py_cmd="conda run --no-capture-output -n {params.py_env} "
            fi
        fi

        $py_cmd python3 scripts/extract_anchor_loci.py \
            --target_name "{wildcards.sample}" \
            --target_fa "{input.target_fa}" \
            --target_faa "{input.target_faa}" \
            --target_gff "{params.target_gff}" \
            --module_id "{params.module_id}" \
            --module_name "{params.module_name}" \
            --target_kos "{params.target_kos}" \
            $anno_arg \
            $gtdb_ko_arg \
            $cache_arg \
            --threads {threads} \
            --anchors_tsv "{input.anchors_tsv}" \
            --refs_tsv "{input.refs_tsv}" \
            --seed_operon_faa "{input.seed_operon_faa}" \
            --expected_locus_length {params.expected_len} \
            --min_edge_bp {params.min_edge} \
            --flank_bp {params.flank_bp} \
            --n_flank_genes {params.n_flank} \
            --max_visualized_refs {params.max_vis_refs} \
            --max_candidate_loci {params.max_candidate_loci} \
            --consensus_contig_threshold {params.consensus_contig_threshold} \
            --significance_mode "{params.significance_mode}" \
            --significance_alpha {params.significance_alpha} \
            --n_permutations {params.n_permutations} \
            --extended_window_bp {params.extended_window_bp} \
            --diamond_env "{params.diamond_env}" \
            --out_topology_tsv "{output.topology_tsv}" \
            --out_target_locus_fa "{output.target_locus_fa}" \
            --out_target_locus_gbk "{output.target_locus_gbk}" \
            --out_refs_locus_dir "{output.refs_locus_dir}" \
            --out_extended_proteins_tsv "{output.extended_proteins_tsv}"
        }} 2>&1 | tee {log}
        """

rule synteny_stage2_rescue_genes:
    input:
        target_fa = get_target_fa,
        target_locus_fa = os.path.join(SYN_RES, "{sample}/{module}/stage1/{sample}_locus.fa"),
        topology_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage1/anchor_topology.tsv"),
        target_faa = get_target_faa,
        query_operon = os.path.join(SYN_RES, "{sample}/{module}/stage0/{sample}_seed_query_operon.faa")
    output:
        gene_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage2/gene_level_rescue_summary.tsv"),
        cluster_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage2/cluster_state.tsv")
    log:
        os.path.join(SYN_LOGS, "stage2/{sample}_{module}.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        target_gff = lambda wildcards: get_target_gff(wildcards),
        target_annotation = lambda wildcards: get_target_annotation(wildcards),
        module_id = lambda wildcards: wildcards.module,
        module_name = lambda wildcards: get_module_name(wildcards),
        blastp_id = lambda wildcards: float(SYNTENY_C.get("blastp_identity", 25.0)),
        blastp_cov = lambda wildcards: float(SYNTENY_C.get("blastp_coverage", 30.0)),
        blastx_evalue = lambda wildcards: float(SYNTENY_C.get("blastx_evalue", 1e-3)),
        max_operon_gap = lambda wildcards: int(SYNTENY_C.get("max_operon_intergenic_bp", 300)),
        alt_kos = lambda wildcards: get_module_alt_kos(wildcards),
        diamond_env = lambda wildcards: syn_env("diamond"),
        py_env = lambda wildcards: syn_env("biopython")
    shell:
        r"""
        {{
        set -euo pipefail

        anno_arg=""
        if [ -n "{params.target_annotation}" ] && [ -f "{params.target_annotation}" ]; then
            anno_arg="--target_annotation_tsv {params.target_annotation}"
        fi

        gff_arg=""
        if [ -n "{params.target_gff}" ] && [ -f "{params.target_gff}" ]; then
            gff_arg="--target_gff {params.target_gff}"
        fi

        alt_kos_arg=""
        if [ -n "{params.alt_kos}" ]; then
            alt_kos_arg="--alternative_kos '{params.alt_kos}'"
        fi

        py_cmd=""
        if [ -n "{params.py_env}" ]; then
            if [ -d "{params.py_env}" ] || [[ "{params.py_env}" == *"/"* ]]; then
                py_cmd="conda run --no-capture-output -p {params.py_env} "
            else
                py_cmd="conda run --no-capture-output -n {params.py_env} "
            fi
        fi

        $py_cmd python3 scripts/rescue_and_classify_gene_state.py \
            --target_name "{wildcards.sample}" \
            --target_locus_fa "{input.target_locus_fa}" \
            --target_fa "{input.target_fa}" \
            --topology_tsv "{input.topology_tsv}" \
            --target_faa "{input.target_faa}" \
            --query_operon_faa "{input.query_operon}" \
            --module_id "{params.module_id}" \
            --module_name "{params.module_name}" \
            $gff_arg \
            $anno_arg \
            $alt_kos_arg \
            --diamond_env "{params.diamond_env}" \
            --blastp_identity {params.blastp_id} \
            --blastp_coverage {params.blastp_cov} \
            --blastx_evalue {params.blastx_evalue} \
            --max_operon_intergenic_bp {params.max_operon_gap} \
            --out_gene_tsv "{output.gene_tsv}" \
            --out_cluster_tsv "{output.cluster_tsv}"
        }} 2>&1 | tee {log}
        """

rule synteny_stage3_validate_absence:
    input:
        target_locus_fa = os.path.join(SYN_RES, "{sample}/{module}/stage1/{sample}_locus.fa"),
        topology_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage1/anchor_topology.tsv"),
        cluster_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage2/cluster_state.tsv"),
        refs_locus_dir = os.path.join(SYN_RES, "{sample}/{module}/stage1/refs_loci")
    output:
        paf = os.path.join(SYN_RES, "{sample}/{module}/stage3/{sample}_vs_ref.paf"),
        locus_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage3/locus_absence_validation.tsv")
    log:
        os.path.join(SYN_LOGS, "stage3/{sample}_{module}.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        reads_bam = lambda wildcards: get_target_bam(wildcards),
        minimap2_env = lambda wildcards: syn_env("mapping"),
        py_env = lambda wildcards: syn_env("biopython"),
        refs_locus_dir = lambda wildcards: os.path.join(SYN_RES, f"{wildcards.sample}/{wildcards.module}/stage1/refs_loci"),
        module_id = lambda wildcards: wildcards.module,
        module_name = lambda wildcards: get_module_name(wildcards)
    shell:
        r"""
        {{
        set -euo pipefail
        py_cmd=""
        if [ -n "{params.py_env}" ]; then
            if [ -d "{params.py_env}" ] || [[ "{params.py_env}" == *"/"* ]]; then
                py_cmd="conda run --no-capture-output -p {params.py_env} "
            else
                py_cmd="conda run --no-capture-output -n {params.py_env} "
            fi
        fi

        $py_cmd python3 scripts/validate_syntenic_absence.py \
            --target_name "{wildcards.sample}" \
            --target_locus_fa "{input.target_locus_fa}" \
            --refs_locus_dir "{params.refs_locus_dir}" \
            --topology_tsv "{input.topology_tsv}" \
            --cluster_tsv "{input.cluster_tsv}" \
            --module_id "{params.module_id}" \
            --module_name "{params.module_name}" \
            --reads_bam "{params.reads_bam}" \
            --minimap2_env "{params.minimap2_env}" \
            --out_paf "{output.paf}" \
            --out_locus_tsv "{output.locus_tsv}"
        }} 2>&1 | tee {log}
        """

rule synteny_stage4_infer_loss:
    input:
        refs_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage0/selected_validated_references.tsv"),
        topology_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage1/anchor_topology.tsv"),
        cluster_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage2/cluster_state.tsv"),
        locus_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage3/locus_absence_validation.tsv")
    output:
        summary_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage4/{sample}_locus_validation_summary.tsv"),
        tree_nwk = os.path.join(SYN_RES, "{sample}/{module}/stage4/{sample}_local_tree.nwk")
    log:
        os.path.join(SYN_LOGS, "stage4/{sample}_{module}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        taxonomy = lambda wildcards: get_syn_input(wildcards.sample, "taxonomy", "d__Archaea;"),
        completeness = lambda wildcards: float(get_syn_input(wildcards.sample, "completeness", 100.0)),
        contamination = lambda wildcards: float(get_syn_input(wildcards.sample, "contamination", 0.0)),
        tree_env = lambda wildcards: syn_env("tree"),
        species_tree = lambda wildcards: syn_full_path(get_syn_input(wildcards.sample, "species_tree_nwk", "")),
        species_outgroup = lambda wildcards: get_syn_input(wildcards.sample, "species_tree_outgroup", ""),
        py_env = lambda wildcards: syn_env("biopython"),
        module_id = lambda wildcards: wildcards.module,
        module_name = lambda wildcards: get_module_name(wildcards),
        localizable_locator_min = lambda wildcards: float(SYNTENY_C.get("operon_fate_classifier", {}).get("localizable_locator_min", 0.4)),
        bilateral_context_min = lambda wildcards: float(SYNTENY_C.get("operon_fate_classifier", {}).get("bilateral_context_min", 0.5)),
        target_retained_min = lambda wildcards: float(SYNTENY_C.get("operon_fate_classifier", {}).get("target_retained_min", 0.6)),
        known_replacement_min = lambda wildcards: float(SYNTENY_C.get("operon_fate_classifier", {}).get("known_replacement_min", 0.6)),
        replacement_absence_max = lambda wildcards: float(SYNTENY_C.get("operon_fate_classifier", {}).get("replacement_absence_max", 0.2)),
        novel_insertion_min = lambda wildcards: float(SYNTENY_C.get("operon_fate_classifier", {}).get("novel_insertion_min", 0.6)),
        loss_confidence_min = lambda wildcards: float(SYNTENY_C.get("operon_fate_classifier", {}).get("loss_confidence_min", 0.6)),
        clean_closure_min = lambda wildcards: float(SYNTENY_C.get("operon_fate_classifier", {}).get("clean_closure_min", 0.6))
    shell:
        r"""
        {{
        set -euo pipefail

        py_cmd=""
        if [ -n "{params.py_env}" ]; then
            if [ -d "{params.py_env}" ] || [[ "{params.py_env}" == *"/"* ]]; then
                py_cmd="conda run --no-capture-output -p {params.py_env} "
            else
                py_cmd="conda run --no-capture-output -n {params.py_env} "
            fi
        fi

        $py_cmd python3 scripts/infer_evolutionary_loss.py \
            --target_name "{wildcards.sample}" \
            --target_taxonomy "{params.taxonomy}" \
            --target_completeness {params.completeness} \
            --target_contamination {params.contamination} \
            --module_id "{params.module_id}" \
            --module_name "{params.module_name}" \
            --refs_tsv "{input.refs_tsv}" \
            --topology_tsv "{input.topology_tsv}" \
            --cluster_tsv "{input.cluster_tsv}" \
            --locus_tsv "{input.locus_tsv}" \
            --localizable_locator_min {params.localizable_locator_min} \
            --bilateral_context_min {params.bilateral_context_min} \
            --target_retained_min {params.target_retained_min} \
            --known_replacement_min {params.known_replacement_min} \
            --replacement_absence_max {params.replacement_absence_max} \
            --novel_insertion_min {params.novel_insertion_min} \
            --loss_confidence_min {params.loss_confidence_min} \
            --clean_closure_min {params.clean_closure_min} \
            --tree_tool_env "{params.tree_env}" \
            --species_tree_nwk "{params.species_tree}" \
            --species_tree_outgroup "{params.species_outgroup}" \
            --out_summary_tsv "{output.summary_tsv}" \
            --out_tree_nwk "{output.tree_nwk}"
        }} 2>&1 | tee {log}
        """

rule synteny_clinker_visualization:
    input:
        target_locus_gbk = os.path.join(SYN_RES, "{sample}/{module}/stage1/{sample}_locus.gbk"),
        topology_tsv = os.path.join(SYN_RES, "{sample}/{module}/stage1/anchor_topology.tsv"),
        refs_locus_dir = os.path.join(SYN_RES, "{sample}/{module}/stage1/refs_loci")
    output:
        html = os.path.join(SYN_RES, "{sample}/{module}/visualization/{sample}_clinker.html"),
        focused_loci = directory(os.path.join(SYN_RES, "{sample}/{module}/visualization/focused_loci"))
    log:
        os.path.join(SYN_LOGS, "clinker/{sample}_{module}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = syn_env("clinker"),
        biopython_env = syn_env("biopython"),
        refs_dir = lambda wildcards: os.path.join(SYN_RES, f"{wildcards.sample}/{wildcards.module}/stage1/refs_loci"),
        pairwise_mode = lambda wildcards: str(SYNTENY_C.get("clinker_pairwise_mode", "all")).lower(),
        context_bp = lambda wildcards: int(SYNTENY_C.get("clinker_context_bp", 10000)),
        flank_genes = lambda wildcards: int(SYNTENY_C.get("clinker_flank_genes", 25)),
        target_kos = lambda wildcards: get_module_kos(wildcards)
    shell:
        r"""
        {{
        set -euo pipefail
        vis_dir=$(dirname "{output.html}")
        mkdir -p "$vis_dir"
        focus_dir="{output.focused_loci}"
        mkdir -p "$focus_dir"
        gbk_count=$(ls -1 "{params.refs_dir}"/*.gbk 2>/dev/null | wc -l || echo 0)
        if [ ! -s "{input.target_locus_gbk}" ] || [ "$gbk_count" -eq 0 ]; then
            echo "[WARNING] Synteny Clinker early exit triggered for sample '{wildcards.sample}', module '{wildcards.module}':" >&2
            echo "          - target_locus_gbk: '{input.target_locus_gbk}' (non-empty: $([ -s '{input.target_locus_gbk}' ] && echo 'YES' || echo 'NO/EMPTY'))" >&2
            echo "          - refs_dir: '{params.refs_dir}' (*.gbk count = $gbk_count)" >&2
            echo "          Writing placeholder HTML and exiting." >&2
            echo '<html><body><p>No validated reference locus is available for this sample.</p></body></html>' > "{output.html}"
            exit 0
        fi

        py_cmd="python3"
        if [ -d "{params.biopython_env}" ] || [[ "{params.biopython_env}" == *"/"* ]]; then
            py_cmd="conda run --no-capture-output -p {params.biopython_env} python3"
        elif [ -n "{params.biopython_env}" ]; then
            py_cmd="conda run --no-capture-output -n {params.biopython_env} python3"
        fi
        $py_cmd scripts/prepare_clinker_inputs.py \
            --target-gbk "{input.target_locus_gbk}" \
            --topology-tsv "{input.topology_tsv}" \
            --refs-dir "{params.refs_dir}" \
            --out-dir "$focus_dir" \
            --module-id "{wildcards.module}" \
            --target-kos "{params.target_kos}" \
            --context-bp {params.context_bp} \
            --flank-genes {params.flank_genes}
        target_focus="$focus_dir/target_locus.focused.gbk"
        if [ ! -s "$target_focus" ]; then
            echo "ERROR: target focused GenBank was not created: $target_focus" >&2
            ls -la "$focus_dir" >&2 || true
            exit 1
        fi

        clinker_cmd="clinker"
        if [ -d "{params.env}" ] || [[ "{params.env}" == *"/"* ]]; then
            clinker_cmd="conda run --no-capture-output -p {params.env} clinker"
        elif [ -n "{params.env}" ]; then
            clinker_cmd="conda run --no-capture-output -n {params.env} clinker"
        fi

        # 1. Generate pairwise Clinker plots based on pairwise_mode ('none', 'nearest_only', 'all')
        if [ "{params.pairwise_mode}" = "all" ]; then
            order_file="$focus_dir/ordered_clinker_files.txt"
            ref_list=""
            if [ -f "$order_file" ] && [ -s "$order_file" ]; then
                ref_list=$(cat "$order_file")
            else
                ref_list=$(ls -1 "$focus_dir"/*_locus.focused.gbk 2>/dev/null | grep -v '_target_locus.focused.gbk' || true)
            fi
            for ref_gbk in $ref_list; do
                if [ -f "$ref_gbk" ]; then
                    ref_base=$(basename "$ref_gbk" _locus.focused.gbk)
                    pair_html="$vis_dir/{wildcards.sample}_vs_${{ref_base}}_clinker.html"
                    ref_target_gbk="$focus_dir/{wildcards.sample}_vs_${{ref_base}}_target_locus.focused.gbk"
                    if [ ! -f "$ref_target_gbk" ]; then
                        ref_target_gbk="$target_focus"
                    fi
                    echo "Generating pairwise Clinker plot: $pair_html"
                    $clinker_cmd "$ref_gbk" "$ref_target_gbk" -p "$pair_html" -i 0.25 || true
                    if [ -f "$pair_html" ]; then
                        python3 scripts/refine_clinker_labels.py "$pair_html" || true
                    fi
                fi
            done
        elif [ "{params.pairwise_mode}" = "nearest_only" ]; then
            nearest_gbk=""
            order_file="$focus_dir/ordered_clinker_files.txt"
            if [ -f "$order_file" ] && [ -s "$order_file" ]; then
                nearest_gbk=$(tail -n 1 "$order_file" | tr -d '\r\n')
            else
                nearest_gbk=$(ls -1 "$focus_dir"/*_locus.focused.gbk 2>/dev/null | grep -v '_target_locus.focused.gbk' | tail -n 1 || true)
            fi
            if [ -n "$nearest_gbk" ] && [ -f "$nearest_gbk" ]; then
                ref_base=$(basename "$nearest_gbk" _locus.focused.gbk)
                pair_html="$vis_dir/{wildcards.sample}_vs_${{ref_base}}_clinker.html"
                ref_target_gbk="$focus_dir/{wildcards.sample}_vs_${{ref_base}}_target_locus.focused.gbk"
                if [ ! -f "$ref_target_gbk" ]; then
                    ref_target_gbk="$target_focus"
                fi
                echo "Generating nearest-ref pairwise Clinker plot: $pair_html"
                $clinker_cmd "$nearest_gbk" "$ref_target_gbk" -p "$pair_html" -i 0.25 || true
                if [ -f "$pair_html" ]; then
                    python3 scripts/refine_clinker_labels.py "$pair_html" || true
                fi
            fi

            # In addition to nearest, automatically generate pairwise plot for any reference matching an alternative contig
            for diff_target_gbk in "$focus_dir"/{wildcards.sample}_vs_*_target_locus.focused.gbk; do
                if [ -f "$diff_target_gbk" ]; then
                    diff_base=$(basename "$diff_target_gbk" _target_locus.focused.gbk | sed "s/^{wildcards.sample}_vs_//")
                    diff_ref_gbk="$focus_dir/${{diff_base}}_locus.focused.gbk"
                    diff_pair_html="$vis_dir/{wildcards.sample}_vs_${{diff_base}}_clinker.html"
                    if [ -f "$diff_ref_gbk" ] && [ ! -f "$diff_pair_html" ]; then
                        echo "Generating cross-contig ref pairwise Clinker plot: $diff_pair_html"
                        $clinker_cmd "$diff_ref_gbk" "$diff_target_gbk" -p "$diff_pair_html" -i 0.25 || true
                        if [ -f "$diff_pair_html" ]; then
                            python3 scripts/refine_clinker_labels.py "$diff_pair_html" || true
                        fi
                    fi
                fi
            done
        fi

        # 2. Generate multi-track Clinker plot containing all references + target
        order_file="$focus_dir/ordered_clinker_files.txt"
        if [ -f "$order_file" ] && [ -s "$order_file" ]; then
            ref_args=$(cat "$order_file" | tr '\n' ' ')
            $clinker_cmd $ref_args "$target_focus" -p "{output.html}" -i 0.25
        elif [ "$gbk_count" -gt 0 ]; then
            focused_refs=$(ls -1 "$focus_dir"/*_locus.focused.gbk 2>/dev/null | grep -v '_target_locus.focused.gbk' | tr '\n' ' ' || true)
            $clinker_cmd $focused_refs "$target_focus" -p "{output.html}" -i 0.25
        else
            $clinker_cmd "$target_focus" -p "{output.html}" -i 0.25
        fi

        # Post-process multi-track Clinker HTML to refine gene labels and eliminate qualifier redundancy
        if [ -f "{output.html}" ]; then
            python3 scripts/refine_clinker_labels.py "{output.html}" || true
        fi
        }} 2>&1 | tee {log}
        """

rule synteny_compile_atp_matrix:
    input:
        summaries = [os.path.join(SYN_RES, f"{s}/{mod}/stage4/{s}_locus_validation_summary.tsv") for s in sorted(SYNTENY_SAMPLES.keys()) for mod in get_sample_modules(s)]
    output:
        stacked_tsv = os.path.join(SYN_RES, "locus_synteny_all_modules_summary.tsv"),
        matrix_tsv = os.path.join(SYN_RES, "atp_synthase_synteny_matrix.tsv"),
        legacy_summary = os.path.join(SYN_RES, "locus_synteny_validation_summary.tsv")
    log:
        os.path.join(SYN_LOGS, "compile_atp_matrix.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        py_env = lambda wildcards: syn_env("biopython")
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p $(dirname "{output.stacked_tsv}")
        if [ -z "{input.summaries}" ]; then
            touch "{output.stacked_tsv}" "{output.matrix_tsv}" "{output.legacy_summary}"
            exit 0
        fi
        
        py_cmd=""
        if [ -n "{params.py_env}" ]; then
            if [ -d "{params.py_env}" ] || [[ "{params.py_env}" == *"/"* ]]; then
                py_cmd="conda run --no-capture-output -p {params.py_env} "
            else
                py_cmd="conda run --no-capture-output -n {params.py_env} "
            fi
        fi

        $py_cmd python3 scripts/compile_atp_synteny_matrix.py \
            --input_tsvs {input.summaries} \
            --out_stacked_tsv "{output.stacked_tsv}" \
            --out_matrix_tsv "{output.matrix_tsv}"

        cp -f "{output.stacked_tsv}" "{output.legacy_summary}" || true
        }} 2>&1 | tee {log}
        """
