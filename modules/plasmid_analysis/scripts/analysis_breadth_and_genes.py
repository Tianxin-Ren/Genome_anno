import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt
import os
import json

def generate_interactive_html(fig, output_path, title="Taxonomic Analysis"):
    """
    Helper to wrap a Plotly figure with custom HTML/JS for filtering, depth control, and navigation.
    """
    trace_data = fig.data[0]
    data_dict = {
        'ids': list(trace_data.ids),
        'labels': list(trace_data.labels),
        'parents': list(trace_data.parents),
        'values': [float(v) for v in trace_data.values],
        'colors': list(trace_data.marker.colors)
    }
    plotly_html = fig.to_html(include_plotlyjs='cdn', full_html=False, div_id='plotly-div')
    custom_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>{{title}}</title>
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 20px; background-color: #f4f7f9; color: #333; }}
            .control-panel {{ 
                margin-bottom: 25px; padding: 20px; background: white; 
                border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.05);
                display: flex; align-items: center; gap: 20px; flex-wrap: wrap;
            }}
            .control-group {{ display: flex; align-items: center; gap: 10px; border-right: 1px solid #eee; padding-right: 20px; }}
            .control-group:last-child {{ border-right: none; }}
            input, select {{ padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; outline: none; }}
            button {{ 
                padding: 8px 18px; cursor: pointer; border: none; 
                border-radius: 6px; background-color: #007bff; color: white; 
                font-weight: 600; transition: all 0.2s;
            }}
            button:hover {{ background-color: #0056b3; transform: translateY(-1px); }}
            .reset-btn {{ background-color: #6c757d; }}
            .reset-btn:hover {{ background-color: #5a6268; }}
            #nodes-info {{ margin-left: auto; font-weight: 500; color: #555; }}
        </style>
    </head>
    <body>
        <h2 style="margin-top: 0; color: #2c3e50;">{{title}} - Advanced Dashboard</h2>
        <div class="control-panel">
            <div class="control-group">
                <strong>Filtering:</strong>
                <label>Min Genomes (N) ≥</label>
                <input type="number" id="min-genomes-input" value="1" min="1">
                <button onclick="applyFilter()">Apply</button>
            </div>
            <div class="control-group">
                <strong>View Depth:</strong>
                <select id="depth-select" onchange="changeDepth(this.value)">
                    <option value="1">1 Level (Domain Only)</option>
                    <option value="2">2 Levels (Phylum)</option>
                    <option value="3" selected>3 Levels (Class)</option>
                    <option value="4">4 Levels (Order)</option>
                    <option value="5">5 Levels (Family)</option>
                    <option value="all">All Levels (Full Genus)</option>
                </select>
            </div>
            <div class="control-group">
                <button class="reset-btn" onclick="resetView()">Reset Top / Back</button>
            </div>
            <div id="nodes-info">Visible Taxa: <span id="nodes-count">-</span></div>
        </div>
        <div style="background: white; padding: 10px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.05);">
            {plotly_html}
        </div>
        <script>
            const masterData = {json.dumps(data_dict)};
            const gd = document.getElementById('plotly-div');
            function applyFilter() {{
                const N = parseInt(document.getElementById('min-genomes-input').value) || 1;
                const validIds = new Set();
                for (let i = 0; i < masterData.ids.length; i++) {{
                    if (masterData.values[i] >= N) {{
                        const id = masterData.ids[i];
                        validIds.add(id);
                        let currentId = id;
                        while(currentId.includes('/')) {{
                            const parentId = currentId.substring(0, currentId.lastIndexOf('/'));
                            validIds.add(parentId);
                            currentId = parentId;
                        }}
                        validIds.add(currentId);
                    }}
                }}
                const f = {{ ids: [], labels: [], parents: [], values: [], colors: [] }};
                for (let i = 0; i < masterData.ids.length; i++) {{
                    if (validIds.has(masterData.ids[i])) {{
                        f.ids.push(masterData.ids[i]);
                        f.labels.push(masterData.labels[i]);
                        f.parents.push(masterData.parents[i]);
                        f.values.push(masterData.values[i]);
                        f.colors.push(masterData.colors[i]);
                    }}
                }}
                Plotly.restyle(gd, {{
                    ids: [f.ids], labels: [f.labels], parents: [f.parents],
                    values: [f.values], 'marker.colors': [f.colors]
                }});
                updateCounter(f.ids.length);
            }}
            function changeDepth(val) {{
                const depth = val === 'all' ? -1 : parseInt(val);
                Plotly.restyle(gd, 'maxdepth', [depth]);
            }}
            function resetView() {{
                Plotly.restyle(gd, {{'level': ''}});
            }}
            function updateCounter(count) {{
                document.getElementById('nodes-count').innerText = count;
            }}
            window.onload = () => updateCounter(masterData.ids.length);
        </script>
    </body>
    </html>
    """
    with open(output_path, 'w') as f:
        f.write(custom_html)

def process_and_visualize(input_file, output_dir="analysis_results", min_genomes=2):
    """
    Process genome annotation results and generate interactive visualizations.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    print(f"Loading data from {input_file}...")
    df = pd.read_csv(input_file, sep='\t')
    
    # --- 2. Taxonomy Processing ---
    print("Processing taxonomy...")
    tax_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    tax_prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    tax_col = 'gtdb_taxonomy' if 'gtdb_taxonomy' in df.columns else 'taxonomy'
    df_tax = df[tax_col].str.split(';', expand=True)
    for i in range(len(tax_levels)):
        if i < df_tax.shape[1]:
            df[tax_levels[i]] = df_tax[i].str.replace(r'^[a-z]__', '', regex=True)
        else:
            df[tax_levels[i]] = 'Unassigned'

    # Identify metadata & gene columns
    metadata_cols = [
        'genome_id', 'gtdb_taxonomy', 'ncbi_genbank_assembly_accession', 
        'sample', 'MAG_id', 'taxonomy', 'genome_size', 'effective_genome_size', 
        'CP', 'CT', 'CP2', 'CT2', 'quality', 'N50', 'num_contigs', 'GC', 
        'Coding_density', 'tRNA_count', '16srRNA_count', '23srRNA_count', 
        '5srRNA_count', 'gs_oxygen', 'gs_pH', 'gs_salinity', 'gs_temperature'
    ]
    try:
        anchor_idx = list(df.columns).index('gs_temperature')
        gene_cols = [c for c in list(df.columns)[anchor_idx + 1:] if c not in tax_levels and c not in metadata_cols]
    except ValueError:
        gene_cols = [c for c in df.columns if c not in metadata_cols and c not in tax_levels]
    
    for col in gene_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    
    df['total_target_genes'] = df[gene_cols].sum(axis=1)
    df['is_zero_hit'] = (df['total_target_genes'] == 0).astype(int)

    # --- 3. Integrated Multi-level Taxonomic Analysis ---
    print("\n--- Performing Integrated Multi-level Analysis ---")
    integration_output_dir = os.path.join(output_dir, "integrated_taxonomy_summary")
    if not os.path.exists(integration_output_dir): os.makedirs(integration_output_dir)

    all_level_stats = []
    for i, level in enumerate(tax_levels):
        current_hierarchy = tax_levels[:i+1]
        level_df = df.groupby(current_hierarchy).agg(
            Total_Genomes=('total_target_genes', 'count'),
            Zero_Hit_Genomes=('is_zero_hit', 'sum')
        ).reset_index()
        
        # Merge taxonomy into single string, stopping at the currently analyzed level
        def format_taxonomy(row):
            parts = []
            for j, lvl in enumerate(current_hierarchy):
                parts.append(f"{tax_prefixes[j]}{row[lvl]}")
            return ";".join(parts)
        
        level_df['Full_Taxonomy'] = level_df.apply(format_taxonomy, axis=1)
        level_df['Analyzed_Level'] = level
        level_df['Zero_Hit_Percentage'] = (level_df['Zero_Hit_Genomes'] / level_df['Total_Genomes']) * 100
        
        # Keep only essential columns
        all_level_stats.append(level_df[['Full_Taxonomy', 'Analyzed_Level', 'Total_Genomes', 'Zero_Hit_Genomes', 'Zero_Hit_Percentage']])

    pd.concat(all_level_stats, ignore_index=True).to_csv(os.path.join(integration_output_dir, "taxonomic_zero_hit_report.tsv"), sep='\t', index=False)

    # (A) Sunburst & Treemap
    sunburst_hierarchy = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
    df_agg = df.groupby(sunburst_hierarchy).agg(Count=('is_zero_hit', 'count'), Zero_Hit_Rate=('is_zero_hit', 'mean')).reset_index()

    for plot_type in ['sunburst', 'treemap']:
        if plot_type == 'sunburst':
            fig = px.sunburst(df_agg, path=sunburst_hierarchy, values='Count', color='Zero_Hit_Rate', color_continuous_scale='RdBu_r')
            fig.update_traces(maxdepth=3)
        else:
            fig = px.treemap(df_agg, path=sunburst_hierarchy, values='Count', color='Zero_Hit_Rate', color_continuous_scale='RdBu_r')
        fig.update_layout(coloraxis_colorbar=dict(title="Zero-Hit Prevalence"))
        generate_interactive_html(fig, os.path.join(integration_output_dir, f"integrated_taxonomic_{plot_type}.html"), title=f"Taxonomic {plot_type.capitalize()}")

    # --- 4. Per-Domain Analysis (Bacteria/Archaea) ---
    domains = [d for d in df['Domain'].unique() if d != 'Unassigned']
    for domain in domains:
        print(f"\nProcessing Domain: {domain}...")
        domain_df = df[df['Domain'] == domain].copy()
        domain_output_dir = os.path.join(output_dir, domain.lower().replace(" ", "_"))
        if not os.path.exists(domain_output_dir): os.makedirs(domain_output_dir)
        fig = px.sunburst(domain_df, path=['Phylum', 'Class', 'Order'], title=f"Taxonomic Breadth: {domain} (N={len(domain_df)})", color_discrete_sequence=px.colors.qualitative.Pastel if domain == 'Bacteria' else px.colors.qualitative.Set2)
        fig.write_html(os.path.join(domain_output_dir, "taxonomic_breadth.html"))
        domain_genes = (domain_df[gene_cols] > 0).astype(int); domain_genes['Phylum'] = domain_df['Phylum']
        phylum_stats = domain_genes.groupby('Phylum').agg(['mean', 'count'])
        keep_phyla_stats = phylum_stats[phylum_stats[gene_cols[0]]['count'] >= min_genomes].copy()
        keep_phyla_stats = keep_phyla_stats.sort_values((gene_cols[0], 'count'), ascending=False)
        if not keep_phyla_stats.empty:
            freq_data = keep_phyla_stats.loc[:, (slice(None), 'mean')] * 100; freq_data.index = [f"{p} (n={int(keep_phyla_stats.loc[p, (gene_cols[0], 'count')])})" for p in freq_data.index]
            try:
                if freq_data.shape[0] >= 2 and freq_data.shape[1] >= 2:
                    sns.clustermap(freq_data, cmap="YlGnBu", figsize=(12, 8), row_cluster=False).savefig(os.path.join(domain_output_dir, "gene_distribution_heatmap.pdf"), bbox_inches='tight')
                else:
                    plt.figure(figsize=(12, 8)); sns.heatmap(freq_data, cmap="YlGnBu"); plt.savefig(os.path.join(domain_output_dir, "gene_distribution_heatmap.pdf"), bbox_inches='tight')
                plt.close('all')
            except Exception: pass
        phylum_zero = domain_df.groupby('Phylum')['is_zero_hit'].agg(['mean', 'count']); phylum_zero = phylum_zero[phylum_zero['count'] >= min_genomes].sort_values('count', ascending=False)
        if not phylum_zero.empty:
            phylum_zero['percent'] = phylum_zero['mean'] * 100
            phylum_zero['label'] = [f"{p} (n={int(c)})" for p, c in zip(phylum_zero.index, phylum_zero['count'])]
            plt.figure(figsize=(10, 8)); sns.barplot(data=phylum_zero.head(20), x='percent', y='label', palette="Reds_r"); plt.title(f"{domain}: Zero-Hit %"); plt.savefig(os.path.join(domain_output_dir, "zero_hit_prevalence.pdf")); plt.close()
        if domain == 'Archaea':
            bathy_df = domain_df[domain_df['Class'] == 'Bathyarchaeia'].copy()
            if not bathy_df.empty:
                bathy_dir = os.path.join(domain_output_dir, "bathyarchaeia_focus"); os.makedirs(bathy_dir, exist_ok=True)
                ba_genes = (bathy_df[gene_cols] > 0).astype(int); ba_genes['Family'] = bathy_df['Family']
                fam_stats = ba_genes.groupby('Family').agg(['mean', 'count']).sort_values((gene_cols[0], 'count'), ascending=False)
                if not fam_stats.empty:
                    f_freq = fam_stats.loc[:, (slice(None), 'mean')] * 100; f_freq.index = [f"{f} (n={int(fam_stats.loc[f, (gene_cols[0], 'count')])})" for f in f_freq.index]
                    plt.figure(figsize=(15, 10)); sns.heatmap(f_freq, cmap="YlGnBu"); plt.title("Bathyarchaeia ATPase Distribution"); plt.savefig(os.path.join(bathy_dir, "family_level_distribution.pdf")); plt.close()
                f_zero = bathy_df.groupby('Family')['is_zero_hit'].agg(['mean', 'count']).sort_values('count', ascending=False)
                if not f_zero.empty:
                    f_zero['percent'] = f_zero['mean'] * 100; f_zero['label'] = [f"{f} (n={int(c)})" for f, c in zip(f_zero.index, f_zero['count'])]
                    plt.figure(figsize=(10, 8)); sns.barplot(data=f_zero.head(20), x='percent', y='label', palette="Reds_r"); plt.title("Bathyarchaeia Family Zero-Hit %"); plt.savefig(os.path.join(bathy_dir, "zero_hit_by_family.pdf")); plt.close()
                bathy_df[bathy_df['is_zero_hit'] == 1][['genome_id', 'gtdb_taxonomy', 'quality']].to_csv(os.path.join(bathy_dir, "zero_hit_genomes.tsv"), sep='\t', index=False)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input"); parser.add_argument("-o", "--output", default="analysis_results"); parser.add_argument("-m", "--min-count", type=int, default=2)
    args = parser.parse_args()
    if os.path.exists(args.input): process_and_visualize(args.input, args.output, args.min_count)
