#!/usr/bin/env python3
import os
import re
import argparse
import pandas as pd

BROAD_SUBTYPE_COLUMNS = [
    "NiFe_G1", "NiFe_G2", "NiFe_G3", "NiFe_G4",
    "FeFe_GA", "FeFe_GB", "FeFe_GC", "Fe_G1",
]

DETAILED_SUBTYPE_COLUMNS = [
    "[FeFe]_Group_A1", "[FeFe]_Group_A2", "[FeFe]_Group_A3", "[FeFe]_Group_A4",
    "[FeFe]_Group_B", "[FeFe]_Group_C1", "[FeFe]_Group_C2", "[FeFe]_Group_C3",
    "[Fe]",
    "[NiFe]_Group_1a", "[NiFe]_Group_1b", "[NiFe]_Group_1c", "[NiFe]_Group_1d",
    "[NiFe]_Group_1e", "[NiFe]_Group_1f", "[NiFe]_Group_1g", "[NiFe]_Group_1h",
    "[NiFe]_Group_1i", "[NiFe]_Group_1k", "[NiFe]_Group_1l",
    "[NiFe]_Group_2a", "[NiFe]_Group_2b", "[NiFe]_Group_2c", "[NiFe]_Group_2d",
    "[NiFe]_Group_2e",
    "[NiFe]_Group_3a", "[NiFe]_Group_3b", "[NiFe]_Group_3c", "[NiFe]_Group_3d",
    "[NiFe]_Group_4a", "[NiFe]_Group_4b", "[NiFe]_Group_4c", "[NiFe]_Group_4d",
    "[NiFe]_Group_4e", "[NiFe]_Group_4f", "[NiFe]_Group_4g", "[NiFe]_Group_4h",
    "[NiFe]_Group_4i",
]

GENOME_OUTPUT_BASE_COLUMNS = [
    "bin", "completeness", "contamination", "GC", "classification", "warnings",
    "total_hydrogenase_genes",
]

DEFAULT_BROAD_SUBTYPE_TSV = "/dssg/home/acct-schhui/share/HyEn_assembly/gses_all_rawMAGs_fa_V0_final.tsv"

def sanitize_id(seq_id):
    """Must match the sanitization rule used in classify_subtypes_placement.py"""
    return re.sub(r'[^a-zA-Z0-9_\-\.]', '_', seq_id)

def parse_query_fasta_mapping(fasta_file):
    """
    Parse original query FASTA to map:
    sanitized_id -> (original_header, genome_id)
    Expects header format containing '[Genome:MAG_NAME]'
    """
    mapping = {}
    print(f"正在解析原始待测 FASTA 文件建立映射关系: {fasta_file} ...")
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                original_header = line.strip()[1:]
                # Extract first word as seq ID (which was sanitized in step 5)
                seq_id = original_header.split()[0]
                sanitized_id_key = sanitize_id(seq_id)
                
                # Extract host genome name from [Genome:xxx]
                genome_match = re.search(r'\[Genome:([^\]]+)\]', original_header)
                genome_id = genome_match.group(1).strip() if genome_match else "UnknownGenome"
                
                mapping[sanitized_id_key] = (original_header, genome_id)
    print(f"成功解析了 {len(mapping)} 条基因与基因组的映射关系。")
    return mapping

def main():
    parser = argparse.ArgumentParser(description="Merge phylogenetic subtyping results with host MAG metadata and GTDB taxonomy.")
    parser.add_argument("--query-fasta", required=True, help="Path to the original query FASTA file (containing [Genome:...] tags)")
    parser.add_argument("--gappa-tsv", required=True, help="Path to Gappa's subtype_classification_per_query.tsv")
    parser.add_argument("--mag-info-tsv", required=True, help="Path to merged_mag_info_gtdbtk.tsv (previously merged MAG metadata + GTDB)")
    parser.add_argument("--broad-subtype-tsv", default=DEFAULT_BROAD_SUBTYPE_TSV, help="Path to genome-level broad hydrogenase subtype table (default: %(default)s)")
    parser.add_argument("--out-gene", required=True, help="Output path for gene-level classification summary TSV")
    parser.add_argument("--out-genome", required=True, help="Output path for genome-level summary TSV with subtype counts")
    parser.add_argument("--lwr-cutoff", type=float, default=0.5, help="Likelihood Weight Ratio (LWR) cutoff for confident classification (default: 0.5)")
    
    args = parser.parse_args()
    
    # 1. Parse original FASTA mapping
    fasta_map = parse_query_fasta_mapping(args.query_fasta)
    
    # 2. Parse Gappa classifications
    print(f"正在读取并过滤 Gappa 亚型分类文件: {args.gappa_tsv} ...")
    df_gappa = pd.read_csv(args.gappa_tsv, sep="\t")
    
    # Filter to keep only the specific leaf node assignments (3 levels in taxopath, e.g. Root;Subtypes;Group)
    df_sub = df_gappa[df_gappa['taxopath'].str.count(';') == 2].copy()
    
    if df_sub.empty:
        print("Error: Gappa outputs do not contain specific subtype assignments (leaves). Please check your taxonomy map.")
        return
        
    # Get the best placement (highest LWR) for each query sequence
    best_hits = df_sub.loc[df_sub.groupby('name')['LWR'].idxmax()].copy()
    best_hits['subtype'] = best_hits['taxopath'].apply(lambda x: x.split(';')[-1])
    
    # 3. Create Gene-level dataframe
    print("正在构建基因水平的亚型分类汇总表...")
    gene_rows = []
    for idx, row in best_hits.iterrows():
        sanitized_name = row['name']
        lwr = row['LWR']
        subtype = row['subtype']
        
        # Resolve original header and host genome
        orig_header, genome_id = fasta_map.get(sanitized_name, (sanitized_name, "UnknownGenome"))
        
        # Determine classification status based on cutoff
        status = subtype if lwr >= args.lwr_cutoff else "Unclassified_Low_Confidence"
        
        gene_rows.append({
            "gene_name": sanitized_name,
            "original_header": orig_header,
            "genome": genome_id,
            "subtype": subtype,
            "confidence": lwr,
            "classification_status": status
        })
    df_genes = pd.DataFrame(gene_rows)
    
    # 4. Load the merged MAG info (already contains metadata and GTDB taxonomy)
    print(f"正在加载宿主基因组元数据表: {args.mag_info_tsv} ...")
    df_mags = pd.read_csv(args.mag_info_tsv, sep="\t", keep_default_na=False)

    print(f"正在加载上一级氢化酶 subtype 汇总表: {args.broad_subtype_tsv} ...")
    df_broad = pd.read_csv(args.broad_subtype_tsv, sep="\t", keep_default_na=False)
    missing_broad_cols = [col for col in ["Genome_Name"] + BROAD_SUBTYPE_COLUMNS if col not in df_broad.columns]
    if missing_broad_cols:
        raise ValueError(f"Broad subtype table is missing required columns: {', '.join(missing_broad_cols)}")
    df_broad = df_broad[["Genome_Name"] + BROAD_SUBTYPE_COLUMNS].copy()
    
    # Merge gene-level data with MAG metadata
    df_genes_merged = pd.merge(df_genes, df_mags, left_on="genome", right_on="bin", how="left")
    
    # Save Gene-level summary
    os.makedirs(os.path.dirname(args.out_gene), exist_ok=True)
    df_genes_merged.to_csv(args.out_gene, sep="\t", index=False)
    print(f"基因水平汇总文件成功保存至: {args.out_gene}")
    
    # 5. Build Genome-level subtype count table
    print("正在计算每个基因组内的各亚型拷贝数，构建基因组水平汇总表...")
    # Count occurrences of confident subtypes per genome
    df_genes_confident = df_genes[df_genes['classification_status'] != "Unclassified_Low_Confidence"]
    
    # Create pivot table of subtype counts per genome
    if not df_genes_confident.empty:
        df_counts = df_genes_confident.pivot_table(
            index="genome", 
            columns="subtype", 
            values="gene_name", 
            aggfunc="count", 
            fill_value=0
        ).reset_index()
        
        # Calculate total hydrogenase genes per genome (including low confidence ones)
        df_total_count = df_genes.groupby("genome")["gene_name"].count().reset_index().rename(columns={"gene_name": "total_hydrogenase_genes"})
        
        # Concatenate formatted list of subtypes for each genome, e.g., "[FeFe]_Group_A3:2, [NiFe]_Group_1a:1"
        def make_subtype_summary(subtypes):
            counts = subtypes.value_counts()
            return ", ".join([f"{sub}:{count}" for sub, count in counts.items()])
            
        df_summary_str = (
            df_genes_confident.groupby("genome")["subtype"]
            .agg(make_subtype_summary)
            .reset_index(name="hydrogenase_subtypes")
        )
        
        # Combine counts and summary string
        df_counts_merged = pd.merge(df_total_count, df_counts, on="genome", how="left")
        df_counts_merged = pd.merge(df_counts_merged, df_summary_str, on="genome", how="left")
        
        # Fill numeric columns with 0, and string columns with "None"
        for col in df_counts_merged.columns:
            if col == "hydrogenase_subtypes":
                df_counts_merged[col] = df_counts_merged[col].fillna("None")
            else:
                df_counts_merged[col] = df_counts_merged[col].fillna(0)
    else:
        # Fallback if no confident classifications
        df_counts_merged = pd.DataFrame(columns=["genome", "total_hydrogenase_genes", "hydrogenase_subtypes"])
        df_counts_merged["genome"] = df_genes["genome"].unique()
        df_counts_merged["total_hydrogenase_genes"] = 0
        df_counts_merged["hydrogenase_subtypes"] = ""
        
    # Merge detailed and broad subtype summaries back into original df_mags
    df_genome_final = pd.merge(df_mags, df_counts_merged, left_on="bin", right_on="genome", how="left")
    df_genome_final = pd.merge(df_genome_final, df_broad, left_on="bin", right_on="Genome_Name", how="left")
    # Clean columns: fill NaNs with 0 for gene counts
    df_genome_final["total_hydrogenase_genes"] = df_genome_final["total_hydrogenase_genes"].fillna(0)
    df_genome_final["hydrogenase_subtypes"] = df_genome_final["hydrogenase_subtypes"].fillna("None")
    
    # Drop temp merge columns if present
    for temp_col in ["genome", "Genome_Name"]:
        if temp_col in df_genome_final.columns:
            df_genome_final.drop(columns=[temp_col], inplace=True)
        
    # Ensure requested subtype columns always exist and are exported in a stable order.
    for col in BROAD_SUBTYPE_COLUMNS + DETAILED_SUBTYPE_COLUMNS:
        if col not in df_genome_final.columns:
            df_genome_final[col] = 0
        df_genome_final[col] = df_genome_final[col].fillna(0).astype(int)
    df_genome_final["total_hydrogenase_genes"] = df_genome_final["total_hydrogenase_genes"].astype(int)

    genome_output_columns = (
        GENOME_OUTPUT_BASE_COLUMNS
        + BROAD_SUBTYPE_COLUMNS
        + DETAILED_SUBTYPE_COLUMNS
        + ["hydrogenase_subtypes"]
    )
    missing_output_cols = [col for col in genome_output_columns if col not in df_genome_final.columns]
    if missing_output_cols:
        raise ValueError(f"MAG info table is missing required output columns: {', '.join(missing_output_cols)}")
    df_genome_final = df_genome_final[genome_output_columns]
    
    # Save Genome-level summary
    os.makedirs(os.path.dirname(args.out_genome), exist_ok=True)
    df_genome_final.to_csv(args.out_genome, sep="\t", index=False)
    print(f"基因组水平汇总文件成功保存至: {args.out_genome}")
    print("\n完成！")

if __name__ == "__main__":
    main()
