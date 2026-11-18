"""Runtime configuration helpers shared by the top-level Snakefile.

The workflow has three independent configuration layers:

1. project: samples, enabled modules, and scientific parameters;
2. platform: paths to pre-installed Conda environments and databases;
3. profile: local/cluster scheduling parameters.

This module translates the new platform schema to the legacy keys currently
consumed by the individual workflow modules.  Keeping that translation here
lets modules migrate incrementally without breaking existing project configs.
"""

from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any, Mapping, MutableMapping


MODULE_NAMES = (
    "preprocess",
    "assembly",
    "binning",
    "annotation",
    "postprocess",
    "unified_annotation",
    "advanced_analysis",
    "hmm_build",
    "plasmid_analysis",
    "abundance",
    "data_prep",
    "phylogeny",
    "visualize",
    "protein_phylogeny",
    "mag_traits",
    "synteny_validation",
    "pathway_reconstruction",
    "mag_annotation",
    "cog_annotation",
    "read_profiling",
    "hydrogenase_annotation",
)

UNIFIED_IMPLEMENTED_STEPS = {
    "genome_qc",
    "taxonomy",
    "cds_prediction",
    "rna_prediction",
    "protein_similarity",
    "orthology",
    "ko_profiles",
    "protein_domains",
    "localization_topology",
    "signal_peptide_prediction",
    "biob_hmm",
    "candidate_integration",
}

UNIFIED_PROTEIN_DEPENDENT_STEPS = {
    "protein_similarity",
    "orthology",
    "ko_profiles",
    "protein_domains",
    "localization_topology",
    "signal_peptide_prediction",
    "biob_hmm",
    "candidate_integration",
}


def deep_merge(*layers: Mapping[str, Any] | None) -> dict[str, Any]:
    """Return a recursive merge; mappings on the right have precedence."""
    result: dict[str, Any] = {}
    for layer in layers:
        if not layer:
            continue
        _merge_into(result, layer)
    return result


def _merge_into(target: MutableMapping[str, Any], source: Mapping[str, Any]) -> None:
    for key, value in source.items():
        if isinstance(value, Mapping) and isinstance(target.get(key), Mapping):
            nested = copy.deepcopy(dict(target[key]))
            _merge_into(nested, value)
            target[key] = nested
        else:
            target[key] = copy.deepcopy(value)


def expand_path(value: Any) -> Any:
    """Expand environment variables and ``~`` in a scalar path."""
    if not isinstance(value, str) or not value:
        return value
    return os.path.expanduser(os.path.expandvars(value))


def project_path(config: Mapping[str, Any], value: str) -> str:
    """Resolve a path relative to the configured project root."""
    value = str(expand_path(value))
    if not value or os.path.isabs(value):
        return value
    root = str(config.get("project_root", ".") or ".")
    return os.path.normpath(os.path.join(root, value))


def _relative_to_project(project_root: str, value: str, fallback: str) -> str:
    value = str(expand_path(value or fallback))
    if os.path.isabs(value):
        try:
            if os.path.commonpath([project_root, value]) == project_root:
                return os.path.relpath(value, project_root)
            return value
        except ValueError:
            return value
    return os.path.normpath(value)


def _rebase_legacy_result(path: Any, results_root: str, module: str) -> str:
    """Map legacy ``results/...`` values onto a custom results root."""
    if not path:
        return os.path.join(results_root, module)
    path = str(expand_path(path))
    if os.path.isabs(path):
        return path
    normalized = os.path.normpath(path)
    parts = Path(normalized).parts
    if parts and parts[0] == "results":
        return os.path.normpath(os.path.join(results_root, *parts[1:]))
    return normalized


def _set_database_aliases(config: MutableMapping[str, Any]) -> None:
    """Expose platform database aliases under legacy module keys.

    Platform values deliberately override module-level paths: database
    locations describe a machine/site, not a scientific task.
    """
    databases = config.get("databases", {}) or {}
    if not isinstance(databases, Mapping):
        raise TypeError("'databases' must be a mapping")
    databases = {key: expand_path(value) for key, value in databases.items()}
    config["databases"] = databases

    def set_nested(target: MutableMapping[str, Any], keys: tuple[str, ...], value: Any) -> None:
        current = target
        for key in keys[:-1]:
            current = current.setdefault(key, {})
        current[keys[-1]] = value

    annotation_aliases = {
        "bakta": "bakta_db",
        "gtdbtk": "gtdbtk_db",
        "checkm2": "checkm2_db",
        "eggnog": "eggnog_db",
        "kofam": "kofam_db",
        "dbcan": "dbcan_db",
        "unified_hmm": "unified_hmm_db",
        "pfam_hmm": "pfam_hmm_db",
        "ko_list": "ko_list",
    }
    annotation = config.setdefault("annotation", {})
    for alias, legacy_key in annotation_aliases.items():
        if databases.get(alias):
            annotation[legacy_key] = expand_path(databases[alias])

    if databases.get("checkm"):
        config.setdefault("binning", {})["checkm_db"] = expand_path(databases["checkm"])
        config.setdefault("mag_traits", {}).setdefault("grodon2", {})["ribo_db"] = expand_path(
            databases["checkm"]
        )
    if databases.get("sortmerna"):
        config.setdefault("abundance", {})["sortmerna_db"] = expand_path(databases["sortmerna"])
    if databases.get("genomad"):
        config.setdefault("plasmid_analysis", {}).setdefault("genomad", {})["db"] = expand_path(
            databases["genomad"]
        )
    if databases.get("dbcan"):
        config.setdefault("mag_traits", {}).setdefault("dbcan", {})["db"] = expand_path(
            databases["dbcan"]
        )
    if databases.get("camper"):
        config.setdefault("mag_traits", {}).setdefault("camper", {})["hmm_db"] = expand_path(
            databases["camper"]
        )

    # COG2024 / CDD Cog_LE annotation aliases
    cog_aliases = {
        "cddid": "cddid",
        "cog_def": "cog_def",
        "cog_fun": "cog_fun",
        "cogle_db": "cogle_db",
        "bitscore_specific": "bitscore_specific",
    }
    cog_anno = config.setdefault("cog_annotation", {}).setdefault("databases", {})
    for alias, target_key in cog_aliases.items():
        if databases.get(alias):
            cog_anno[target_key] = expand_path(databases[alias])
    if databases.get("cdd"):
        cdd_dir = str(expand_path(databases["cdd"]))
        if not cog_anno.get("cddid"):
            cand = os.path.join(cdd_dir, "cddid.tbl")
            if os.path.exists(cand):
                cog_anno["cddid"] = cand
        if not cog_anno.get("cogle_db"):
            for cand_name in ("Cog_LE", "Cog"):
                cand = os.path.join(cdd_dir, cand_name)
                if os.path.exists(cand):
                    cog_anno["cogle_db"] = cand
                    break
        if not cog_anno.get("bitscore_specific"):
            cand = os.path.join(cdd_dir, "bitscore_specific.txt")
            if os.path.exists(cand):
                cog_anno["bitscore_specific"] = cand

    # Unified annotation has a deeper configuration contract. Keep the public
    # platform aliases flat while mapping them to the nested rule parameters.
    unified = config.setdefault("unified_annotation", {})
    if databases.get("root"):
        config["database_root"] = databases["root"]
        unified["database_root"] = databases["root"]
    unified_aliases = {
        "checkm2": ("quality_and_taxonomy", "completeness_contamination", "database"),
        "gtdbtk": ("quality_and_taxonomy", "taxonomy", "database"),
        "rfam": ("rna_prediction", "rfam", "database"),
        "swissprot": ("protein_function", "similarity", "databases", "swissprot"),
        "nr": ("protein_function", "similarity", "databases", "nr"),
        "uniref90": ("protein_function", "similarity", "databases", "uniref90_optional"),
        "refseq_protein": ("protein_function", "similarity", "databases", "refseq_protein_optional"),
        "eggnog": ("protein_function", "orthology", "database"),
        "kofam_profiles": ("protein_function", "ko_profiles", "profiles"),
        "ko_list": ("protein_function", "ko_profiles", "ko_list"),
        "pfam": ("protein_function", "domains", "existing_member_databases", "pfam"),
        "cdd": ("protein_function", "domains", "existing_member_databases", "cdd"),
        "biob_hmm": ("protein_function", "biob_hmm", "hmm_db"),
        "biob_reference_proteins": ("protein_function", "biob_hmm", "reference_proteins"),
        "biob_bacteria_refs": ("deep_disambiguation", "references_by_domain", "bacteria"),
        "biob_archaea_refs": ("deep_disambiguation", "references_by_domain", "archaea"),
        "bakta": ("reference_checks", "bakta", "database"),
    }
    for alias, nested_keys in unified_aliases.items():
        if databases.get(alias):
            set_nested(unified, nested_keys, databases[alias])


def _set_unified_environment_aliases(config: MutableMapping[str, Any]) -> None:
    environments = config.get("environments", {}) or {}
    unified = config.setdefault("unified_annotation", {})

    def set_if_present(alias: str, keys: tuple[str, ...]) -> None:
        if not environments.get(alias):
            return
        current = unified
        for key in keys[:-1]:
            current = current.setdefault(key, {})
        current[keys[-1]] = environments[alias]

    aliases = {
        "quast": ("quality_and_taxonomy", "assembly_stats", "env"),
        "checkm2": ("quality_and_taxonomy", "completeness_contamination", "env"),
        "gtdbtk": ("quality_and_taxonomy", "taxonomy", "env"),
        "pyrodigal": ("cds_prediction", "env"),
        "trnascan": ("rna_prediction", "trna", "env"),
        "infernal": ("rna_prediction", "rfam", "env"),
        "aragorn": ("rna_prediction", "tmrna", "env"),
        "mmseqs2": ("protein_function", "similarity", "env"),
        "eggnog": ("protein_function", "orthology", "env"),
        "kofamscan": ("protein_function", "ko_profiles", "env"),
        "interproscan": ("protein_function", "domains", "env"),
        "signalp6": ("protein_function", "signal_peptide", "env"),
        "deeptmhmm": ("protein_function", "topology", "env"),
        "hmmer": ("protein_function", "biob_hmm", "env"),
        "mafft": ("deep_disambiguation", "alignment", "env"),
        "iqtree": ("deep_disambiguation", "phylogeny", "env"),
        "foldseek": ("deep_disambiguation", "structure_search", "env"),
        "bakta": ("reference_checks", "bakta", "env"),
        "prokka": ("reference_checks", "prokka", "env"),
    }
    for alias, nested_keys in aliases.items():
        set_if_present(alias, nested_keys)


def validate_config_contract(config: Mapping[str, Any]) -> None:
    """Reject active settings that the current workflow cannot honor."""
    module_info = (config.get("modules", {}) or {}).get("unified_annotation", {}) or {}
    if not isinstance(module_info, Mapping) or not module_info.get("enable", False):
        return

    steps = (config.get("unified_annotation", {}) or {}).get("steps", {}) or {}
    if not isinstance(steps, Mapping):
        raise TypeError("'unified_annotation.steps' must be a mapping")

    unsupported = sorted(
        name for name, enabled in steps.items()
        if enabled and name not in UNIFIED_IMPLEMENTED_STEPS
    )
    if unsupported:
        raise ValueError(
            "Unified annotation steps are enabled but not implemented: "
            + ", ".join(unsupported)
        )

    if not steps.get("cds_prediction", True) and any(
        steps.get(name, False) for name in UNIFIED_PROTEIN_DEPENDENT_STEPS
    ):
        raise ValueError(
            "unified_annotation.steps.cds_prediction must be true when protein evidence "
            "or candidate_integration is enabled"
        )


def apply_runtime_config(config: MutableMapping[str, Any]) -> MutableMapping[str, Any]:
    """Normalize the layered schema and populate legacy compatibility keys."""
    paths = config.setdefault("paths", {})
    if not isinstance(paths, MutableMapping):
        raise TypeError("'paths' must be a mapping")

    project_root = expand_path(paths.get("project_root") or config.get("project_root") or ".")
    project_root = os.path.abspath(str(project_root))
    config["project_root"] = project_root
    paths["project_root"] = project_root

    conda_envs = expand_path(paths.get("conda_envs") or config.get("conda_env_dir") or "")
    config["conda_env_dir"] = str(conda_envs)
    paths["conda_envs"] = str(conda_envs)

    results_root = _relative_to_project(
        project_root, paths.get("results") or config.get("results_dir", "results"), "results"
    )
    logs_root = _relative_to_project(
        project_root, paths.get("logs") or config.get("logs_dir", "logs"), "logs"
    )
    work_root = _relative_to_project(
        project_root, paths.get("work") or config.get("work_dir", "work"), "work"
    )
    paths.update({"results": results_root, "logs": logs_root, "work": work_root})
    config.update({"results_dir": results_root, "logs_dir": logs_root, "work_dir": work_root})

    global_envs = config.get("environments", {}) or {}
    if not isinstance(global_envs, Mapping):
        raise TypeError("'environments' must be a mapping")

    for module in MODULE_NAMES:
        module_config = config.setdefault(module, {})
        if not isinstance(module_config, MutableMapping):
            continue
        module_config["results_dir"] = _rebase_legacy_result(
            module_config.get("results_dir"), results_root, module
        )
        module_config.setdefault("logs_dir", os.path.join(logs_root, module))
        module_config.setdefault("tmp_dir", os.path.join(work_root, module))
        module_envs = module_config.get("conda_envs", {}) or {}
        # Global aliases are platform-specific and therefore take precedence.
        module_config["conda_envs"] = deep_merge(module_envs, global_envs)
        if "gtdbtk" in module_config and isinstance(module_config["gtdbtk"], MutableMapping):
            gtdb_sub = module_config["gtdbtk"]
            if global_envs.get("gtdbtk"):
                gtdb_envs = gtdb_sub.get("conda_envs", {}) or {}
                if isinstance(gtdb_envs, MutableMapping):
                    gtdb_envs["gtdbtk"] = global_envs["gtdbtk"]
                    gtdb_sub["conda_envs"] = gtdb_envs

    annotation = config.setdefault("annotation", {})
    annotation["external_results_dir"] = _rebase_legacy_result(
        annotation.get("external_results_dir", "results/external_annotation"),
        results_root,
        "external_annotation",
    )

    _set_database_aliases(config)
    _set_unified_environment_aliases(config)
    validate_config_contract(config)
    return config


def unresolved_variables(config: Mapping[str, Any]) -> list[str]:
    """Return config locations that still contain an unexpanded ``${VAR}``."""
    unresolved: list[str] = []

    def walk(value: Any, prefix: str) -> None:
        if isinstance(value, Mapping):
            for key, child in value.items():
                walk(child, f"{prefix}.{key}" if prefix else str(key))
        elif isinstance(value, list):
            for index, child in enumerate(value):
                walk(child, f"{prefix}[{index}]")
        elif isinstance(value, str) and "${" in value:
            unresolved.append(prefix)

    walk(config, "")
    return unresolved
