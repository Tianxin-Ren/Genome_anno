#!/usr/bin/env python3
"""
Summarize Mapping Abundance by Target Genes

This script:
1. For each sample & MAG: Create a table of abundance metrics (rows) x target genes (cols)
2. For each sample & MAG type: Aggregate abundance across MAGs of each type

Usage:
    python summarize_abundance_by_gene.py <base_dir> <output_dir>

Expected directory structure:
    base_dir/
    ├── abundance/           # {sample}_ratio.tsv files
    ├── annotation/          # {MAG}.unique_sampled.tsv files
    ├── all_lab.tsv          # Target genes (knum, Symbol)
    └── type.tsv             # MAG type mapping

Output:
    output_dir/
    ├── by_mag/
    │   └── {sample}/
    │       └── {MAG}.tsv    # Rows: metrics, Cols: target genes
    └── by_type/
        └── {sample}/
            └── {type}.tsv   # Rows: metrics, Cols: target genes (aggregated)
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict


# Abundance metrics to extract
ABUNDANCE_METRICS = [
    'DNA_Counts', 'DNA_RPK', 'DNA_TPM_Mapped', 'DNA_RPKM_Mapped', 
    'DNA_RPKM_Lib', 'DNA_TPM_Lib',
    'RNA_Counts', 'RNA_RPK', 'RNA_TPM_Mapped', 'RNA_RPKM_Mapped',
    'RNA_RPKM_Lib', 'RNA_TPM_Lib', 'Ratio'
]


def load_target_genes(all_lab_path):
    """Load target K-numbers from all_lab.tsv"""
    df = pd.read_csv(all_lab_path, sep='\t')
    # Create knum -> symbol mapping
    knum_to_symbol = dict(zip(df['knum'], df['Symbol']))
    target_knums = set(df['knum'])
    print(f"Loaded {len(target_knums)} target genes from {all_lab_path}")
    return target_knums, knum_to_symbol


def load_mag_types(type_path):
    """Load MAG -> type mapping from type.tsv"""
    df = pd.read_csv(type_path, sep='\t')
    # Assume columns: MAG, type (adjust if different)
    if 'MAG' not in df.columns:
        # Try first two columns
        df.columns = ['MAG', 'type'] + list(df.columns[2:])
    mag_to_type = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
    print(f"Loaded {len(mag_to_type)} MAG type mappings")
    return mag_to_type


def load_annotation(annotation_path, target_knums):
    """
    Load annotation file and filter for target K-numbers.
    Returns dict: LocalGene -> K-number
    """
    # Annotation format: K-number, ?, LocalGene, ...
    gene_to_knum = {}
    with open(annotation_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                knum = parts[0]
                local_gene = parts[2]
                if knum in target_knums:
                    gene_to_knum[local_gene] = knum
    return gene_to_knum


def load_abundance(abundance_path, mag_filter=None):
    """Load abundance file, optionally filter by MAG"""
    df = pd.read_csv(abundance_path, sep='\t')
    if mag_filter:
        df = df[df['MAG'] == mag_filter]
    return df


def summarize_mag_abundance(abundance_df, gene_to_knum, knum_to_symbol):
    """
    Create summary table: rows = metrics, cols = target genes
    Aggregates by K-number (sum for counts, mean for normalized metrics)
    """
    # Filter abundance rows that have target genes
    abundance_df = abundance_df[abundance_df['LocalGene'].isin(gene_to_knum.keys())].copy()
    
    if abundance_df.empty:
        return None
    
    # Add K-number column
    abundance_df['knum'] = abundance_df['LocalGene'].map(gene_to_knum)
    
    # Aggregate by K-number
    agg_funcs = {}
    for metric in ABUNDANCE_METRICS:
        if metric in abundance_df.columns and metric != 'Ratio':
            if 'Counts' in metric:
                agg_funcs[metric] = 'sum'
            else:
                agg_funcs[metric] = 'sum'  # Sum all normalized metrics too
    
    if not agg_funcs:
        return None
        
    grouped = abundance_df.groupby('knum').agg(agg_funcs).T
    
    # Recalculate Ratio for each K-number
    if 'RNA_TPM_Lib' in grouped.index and 'DNA_TPM_Lib' in grouped.index:
        dna_sums = grouped.loc['DNA_TPM_Lib']
        rna_sums = grouped.loc['RNA_TPM_Lib']
        
        # Use a vectorized way to set NA where DNA is 0
        ratios = []
        for d, r in zip(dna_sums, rna_sums):
            if d > 0:
                ratios.append(r / d)
            else:
                ratios.append("NA")
        grouped.loc['Ratio'] = ratios
    
    # Rename columns to include symbol
    new_cols = {}
    for knum in grouped.columns:
        symbol = knum_to_symbol.get(knum, '')
        new_cols[knum] = f"{knum}|{symbol}" if symbol else knum
    grouped = grouped.rename(columns=new_cols)
    
    return grouped


def process_sample(sample_name, abundance_path, annotation_dir, 
                   target_knums, knum_to_symbol, mag_to_type, output_dir):
    """Process a single sample, generate both by-MAG and by-type summaries"""
    
    print(f"\nProcessing sample: {sample_name}")
    
    # Load full abundance file
    abundance_df = load_abundance(abundance_path)
    mags_in_sample = abundance_df['MAG'].unique()
    print(f"  Found {len(mags_in_sample)} MAGs in abundance file")
    
    # Output directories
    mag_output_dir = output_dir / 'by_mag' / sample_name
    type_output_dir = output_dir / 'by_type' / sample_name
    mag_output_dir.mkdir(parents=True, exist_ok=True)
    type_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Store summaries for type aggregation
    type_summaries = defaultdict(list)
    
    # Process each MAG
    for mag in mags_in_sample:
        annotation_file = annotation_dir / f"{mag}.unique_sampled.tsv"
        
        if not annotation_file.exists():
            print(f"  Warning: Annotation file not found for {mag}")
            continue
        
        # Load annotation for this MAG
        gene_to_knum = load_annotation(annotation_file, target_knums)
        if not gene_to_knum:
            print(f"  No target genes found in annotation for {mag}")
            continue
        
        # Filter abundance for this MAG
        mag_abundance = abundance_df[abundance_df['MAG'] == mag]
        
        # Generate summary
        summary = summarize_mag_abundance(mag_abundance, gene_to_knum, knum_to_symbol)
        
        if summary is not None:
            # Save by-MAG summary
            output_path = mag_output_dir / f"{mag}.tsv"
            summary.to_csv(output_path, sep='\t')
            print(f"  Saved: {output_path}")
            
            # Store for type aggregation
            mag_type = mag_to_type.get(mag, 'unknown')
            type_summaries[mag_type].append(summary)
    
    # Aggregate by type
    for mag_type, summaries in type_summaries.items():
        if summaries:
            # Combine all MAG summaries of this type
            combined = pd.concat(summaries, axis=1)
            # Remove Ratio row before summing to avoid string type errors
            if 'Ratio' in combined.index:
                combined = combined.drop('Ratio')
            # Sum across same columns (same K-numbers)
            combined = combined.groupby(combined.columns, axis=1).sum()
            
            # Recalculate Ratio for the type
            if 'RNA_TPM_Lib' in combined.index and 'DNA_TPM_Lib' in combined.index:
                dna_sums = combined.loc['DNA_TPM_Lib']
                rna_sums = combined.loc['RNA_TPM_Lib']
                
                ratios = []
                for d, r in zip(dna_sums, rna_sums):
                    if d > 0:
                        ratios.append(r / d)
                    else:
                        ratios.append("NA")
                combined.loc['Ratio'] = ratios
            output_path = type_output_dir / f"{mag_type}.tsv"
            combined.to_csv(output_path, sep='\t')
            print(f"  Saved type summary: {output_path}")


def merge_all_results(output_dir):
    """
    Merge all individual results into master tables:
    1. master_by_mag.tsv - Long format: Sample, MAG, Metric, Gene, Value
    2. master_by_type.tsv - Long format: Sample, Type, Metric, Gene, Value
    """
    print("\n" + "="*50)
    print("Merging all results into master tables...")
    
    output_dir = Path(output_dir)
    
    # Merge by_mag results
    by_mag_dir = output_dir / 'by_mag'
    mag_records = []
    
    if by_mag_dir.exists():
        for sample_dir in by_mag_dir.iterdir():
            if sample_dir.is_dir():
                sample_name = sample_dir.name
                for mag_file in sample_dir.glob('*.tsv'):
                    mag_name = mag_file.stem
                    df = pd.read_csv(mag_file, sep='\t', index_col=0)
                    # Convert wide to long format
                    for metric in df.index:
                        for gene in df.columns:
                            mag_records.append({
                                'Sample': sample_name,
                                'MAG': mag_name,
                                'Metric': metric,
                                'Gene': gene,
                                'Value': df.loc[metric, gene]
                            })
    
    if mag_records:
        master_mag = pd.DataFrame(mag_records)
        master_mag_path = output_dir / 'master_by_mag.tsv'
        master_mag.to_csv(master_mag_path, sep='\t', index=False)
        print(f"  Saved: {master_mag_path} ({len(mag_records)} records)")
        
        # Also create wide format pivot table for easy viewing
        # Rows: Sample|MAG|Metric, Columns: Genes
        pivot_mag = master_mag.pivot_table(
            index=['Sample', 'MAG', 'Metric'],
            columns='Gene',
            values='Value',
            aggfunc='sum'
        )
        pivot_mag_path = output_dir / 'master_by_mag_wide.tsv'
        pivot_mag.to_csv(pivot_mag_path, sep='\t')
        print(f"  Saved: {pivot_mag_path}")
    
    # Merge by_type results
    by_type_dir = output_dir / 'by_type'
    type_records = []
    
    if by_type_dir.exists():
        for sample_dir in by_type_dir.iterdir():
            if sample_dir.is_dir():
                sample_name = sample_dir.name
                for type_file in sample_dir.glob('*.tsv'):
                    type_name = type_file.stem
                    df = pd.read_csv(type_file, sep='\t', index_col=0)
                    for metric in df.index:
                        for gene in df.columns:
                            type_records.append({
                                'Sample': sample_name,
                                'Type': type_name,
                                'Metric': metric,
                                'Gene': gene,
                                'Value': df.loc[metric, gene]
                            })
    
    if type_records:
        master_type = pd.DataFrame(type_records)
        master_type_path = output_dir / 'master_by_type.tsv'
        master_type.to_csv(master_type_path, sep='\t', index=False)
        print(f"  Saved: {master_type_path} ({len(type_records)} records)")
        
        # Wide format pivot
        pivot_type = master_type.pivot_table(
            index=['Sample', 'Type', 'Metric'],
            columns='Gene',
            values='Value',
            aggfunc='sum'
        )
        pivot_type_path = output_dir / 'master_by_type_wide.tsv'
        pivot_type.to_csv(pivot_type_path, sep='\t')
        print(f"  Saved: {pivot_type_path}")


def main():
    parser = argparse.ArgumentParser(description='Summarize abundance by target genes.')
    parser.add_argument('--abundance-dir', required=True, help='Directory containing ratio.tsv files')
    parser.add_argument('--annotation-dir', required=True, help='Directory containing unique_sampled.tsv files')
    parser.add_argument('--target-genes', required=True, help='Path to target_gene.tsv')
    parser.add_argument('--mag-types', required=True, help='Path to type.tsv')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    args = parser.parse_args()
    
    abundance_dir = Path(args.abundance_dir)
    annotation_dir = Path(args.annotation_dir)
    target_genes_path = Path(args.target_genes)
    type_path = Path(args.mag_types)
    output_dir = Path(args.output_dir)
    
    # Validate paths
    for path in [abundance_dir, annotation_dir, target_genes_path, type_path]:
        if not path.exists():
            print(f"Error: {path} not found")
            sys.exit(1)
    
    # Load reference data
    target_knums, knum_to_symbol = load_target_genes(target_genes_path)
    mag_to_type = load_mag_types(type_path)
    
    # Process each sample
    # Search for both .ratio.tsv and _ratio.tsv for compatibility
    for ratio_file in list(abundance_dir.glob('*.ratio.tsv')) + list(abundance_dir.glob('*_ratio.tsv')):
        sample_name = ratio_file.stem.replace('.ratio', '').replace('_ratio', '')
        process_sample(
            sample_name, ratio_file, annotation_dir,
            target_knums, knum_to_symbol, mag_to_type, output_dir
        )
    
    # Merge all results into master tables
    merge_all_results(output_dir)
    
    print("\n" + "="*50)
    print("Processing complete!")
    print(f"Results saved to: {output_dir}")
    print("\nMaster files:")
    print(f"  - {output_dir}/master_by_mag.tsv (long format)")
    print(f"  - {output_dir}/master_by_mag_wide.tsv (wide format)")
    print(f"  - {output_dir}/master_by_type.tsv (long format)")
    print(f"  - {output_dir}/master_by_type_wide.tsv (wide format)")


if __name__ == '__main__':
    main()

