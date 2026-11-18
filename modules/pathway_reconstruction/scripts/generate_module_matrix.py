#!/usr/bin/env python3
"""
generate_module_matrix.py
Aggregates module evaluations into Level 3 matrices and Level 4 energy/metabolic capabilities profile.
Outputs:
  - genome_module_matrix.tsv.gz: Genome x Module completeness matrix
  - genome_metabolism_capabilities.json: High-level energy & metabolic interpretation for BVS
"""

import argparse
import gzip
import json
import os
import sys

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def main():
    parser = argparse.ArgumentParser(description="Generate module matrices and Level 4 metabolic capability profile.")
    parser.add_argument("-c", "--completeness_tsv", required=True, help="Path to genome_module_completeness.tsv.gz")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    matrix_path = os.path.join(args.out_dir, "genome_module_matrix.tsv")
    capabilities_json = os.path.join(args.out_dir, "genome_metabolism_capabilities.json")

    # Read completeness records
    # (genome_id, module_id) -> (completeness, status)
    data = {}
    all_genomes = set()
    all_modules = set()

    if os.path.isfile(args.completeness_tsv):
        with open_file(args.completeness_tsv, "rt") as f:
            first_line = f.readline()
            if first_line:
                header = first_line.strip().split("\t")
                if "genome_id" in header and "module_id" in header and "completeness" in header and "status" in header:
                    g_idx = header.index("genome_id")
                    m_idx = header.index("module_id")
                    c_idx = header.index("completeness")
                    s_idx = header.index("status")

                    for line in f:
                        parts = line.strip().split("\t")
                        if len(parts) > max(g_idx, m_idx, c_idx, s_idx):
                            gid = parts[g_idx]
                            mid = parts[m_idx]
                            comp = float(parts[c_idx])
                            stat = parts[s_idx]

                            all_genomes.add(gid)
                            all_modules.add(mid)
                            data[(gid, mid)] = (comp, stat)

    sorted_genomes = sorted(all_genomes)
    sorted_modules = sorted(all_modules)

    # 1. Output Genome x Module Completeness Matrix
    print(f"Writing {len(sorted_genomes)} x {len(sorted_modules)} completeness matrix to {matrix_path}...")
    with open(matrix_path, "w", encoding="utf-8") as out:

        out.write("genome_id\t" + "\t".join(sorted_modules) + "\n")
        for gid in sorted_genomes:
            row = [str(round(data.get((gid, mid), (0.0, "absent"))[0], 4)) for mid in sorted_modules]
            out.write(gid + "\t" + "\t".join(row) + "\n")

    # 2. Level 4 Metabolic Capability Profile
    print(f"Constructing Level 4 Bioenergetics & Metabolic Capability profile -> {capabilities_json}...")
    genome_capabilities = []

    for gid in sorted_genomes:
        # Helper to get status
        def get_stat(mid):
            return data.get((gid, mid), (0.0, "absent"))[1]
        def get_comp(mid):
            return data.get((gid, mid), (0.0, "absent"))[0]

        cap = {
            "genome_id": gid,
            "energy_modules": {
                "ATP_synthesis": {
                    "V_A_type_M00159":                        {"status": get_stat("M00159"), "completeness": get_comp("M00159")},
                    "F_type_M00157":                          {"status": get_stat("M00157"), "completeness": get_comp("M00157")},
                    "NADH_quinone_oxidoreductase_M00144":     {"status": get_stat("M00144"), "completeness": get_comp("M00144")},
                    "succinate_dehydrogenase_M00149":         {"status": get_stat("M00149"), "completeness": get_comp("M00149")},
                    "cytochrome_bc1_complex_M00151":          {"status": get_stat("M00151"), "completeness": get_comp("M00151")},
                    "cytochrome_c_oxidase_prok_M00155":       {"status": get_stat("M00155"), "completeness": get_comp("M00155")},
                    "cytochrome_bd_ubiquinol_oxidase_M00153": {"status": get_stat("M00153"), "completeness": get_comp("M00153")}
                },
                "methanogenesis": {
                    "CO2_to_methane_M00567":                      {"status": get_stat("M00567"), "completeness": get_comp("M00567")},
                    "methanol_to_methane_M00356":                 {"status": get_stat("M00356"), "completeness": get_comp("M00356")},
                    "methylamine_to_methane_M00563":              {"status": get_stat("M00563"), "completeness": get_comp("M00563")},
                    "acetate_to_methane_M00357":                  {"status": get_stat("M00357"), "completeness": get_comp("M00357")},
                    "methane_oxidation_methanotroph_M00174":      {"status": get_stat("M00174"), "completeness": get_comp("M00174")},
                    "coenzyme_M_biosynthesis_M00358":             {"status": get_stat("M00358"), "completeness": get_comp("M00358")},
                    "coenzyme_F420_biosynthesis_archaea_M00378":   {"status": get_stat("M00378"), "completeness": get_comp("M00378")},
                    "coenzyme_F430_biosynthesis_M00836":          {"status": get_stat("M00836"), "completeness": get_comp("M00836")},
                    "methanofuran_biosynthesis_M00935":           {"status": get_stat("M00935"), "completeness": get_comp("M00935")}
                },
                "nitrogen_metabolism": {
                    "nitrogen_fixation_M00175":                    {"status": get_stat("M00175"), "completeness": get_comp("M00175")},
                    "nitrification_ammonia_to_nitrite_M00528":     {"status": get_stat("M00528"), "completeness": get_comp("M00528")},
                    "complete_nitrification_comammox_M00804":      {"status": get_stat("M00804"), "completeness": get_comp("M00804")},
                    "dissimilatory_nitrate_reduction_DNRA_M00530": {"status": get_stat("M00530"), "completeness": get_comp("M00530")},
                    "assimilatory_nitrate_reduction_M00531":       {"status": get_stat("M00531"), "completeness": get_comp("M00531")},
                    "denitrification_M00529":                      {"status": get_stat("M00529"), "completeness": get_comp("M00529")},
                    "anammox_M00973":                              {"status": get_stat("M00973"), "completeness": get_comp("M00973")}
                },
                "sulfur_metabolism": {
                    "dissimilatory_sulfate_reduction_M00596": {"status": get_stat("M00596"), "completeness": get_comp("M00596")},
                    "assimilatory_sulfate_reduction_M00176":  {"status": get_stat("M00176"), "completeness": get_comp("M00176")},
                    "thiosulfate_oxidation_sox_M00595":       {"status": get_stat("M00595"), "completeness": get_comp("M00595")},
                    "sulfide_oxidation_M00985":               {"status": get_stat("M00985"), "completeness": get_comp("M00985")},
                    "sulfur_reduction_M00986":                {"status": get_stat("M00986"), "completeness": get_comp("M00986")},
                    "thiosulfate_oxidation_S4I_M00984":       {"status": get_stat("M00984"), "completeness": get_comp("M00984")}
                },
                "carbon_fixation": {
                    "wood_ljungdahl_reductive_acetyl_CoA_M00377":       {"status": get_stat("M00377"), "completeness": get_comp("M00377")},
                    "hydroxypropionate_hydroxybutylate_cycle_M00375":   {"status": get_stat("M00375"), "completeness": get_comp("M00375")},
                    "dicarboxylate_hydroxybutyrate_cycle_M00374":       {"status": get_stat("M00374"), "completeness": get_comp("M00374")},
                    "reductive_citrate_cycle_rTCA_M00173":              {"status": get_stat("M00173"), "completeness": get_comp("M00173")},
                    "hydroxypropionate_bicycle_M00376":                 {"status": get_stat("M00376"), "completeness": get_comp("M00376")},
                    "calvin_cycle_reductive_pentose_phosphate_M00165":  {"status": get_stat("M00165"), "completeness": get_comp("M00165")},
                    "phosphate_acetyltransferase_acetate_kinase_M00579": {"status": get_stat("M00579"), "completeness": get_comp("M00579")}
                },
                "fermentation_and_central_carbon": {
                    "pyruvate_oxidation_M00307":                      {"status": get_stat("M00307"), "completeness": get_comp("M00307")},
                    "pentose_phosphate_pathway_archaea_M00580":       {"status": get_stat("M00580"), "completeness": get_comp("M00580")},
                    "semi_phosphorylative_ED_pathway_M00308":         {"status": get_stat("M00308"), "completeness": get_comp("M00308")},
                    "non_phosphorylative_ED_pathway_M00309":          {"status": get_stat("M00309"), "completeness": get_comp("M00309")},
                    "pentose_bisphosphate_pathway_archaea_M00968":   {"status": get_stat("M00968"), "completeness": get_comp("M00968")}
                }
            }
        }
        genome_capabilities.append(cap)

    with open(capabilities_json, "w", encoding="utf-8") as jf:
        json.dump(genome_capabilities, jf, indent=2)

    print("Matrix and metabolic capabilities profile successfully generated.")

if __name__ == "__main__":
    main()
