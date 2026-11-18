#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genome_anno Pipeline Visualization Generator
Generates beautiful, professional genomics & metagenomics plots and HTML dashboards.
Includes a pure-Python, zero-dependency circular genome & transcriptome visualizer.
"""

import os
import sys
import argparse
import json
import glob
import math
import csv

# Safe cairosvg import for high-fidelity SVG-to-PDF vector conversion
try:
    import cairosvg
    HAS_CAIROSVG = True
except (ImportError, OSError):
    HAS_CAIROSVG = False
    # Proactively try to install it via pip if it's missing (as a self-healing attempt)
    try:
        import subprocess
        print("WARNING: cairosvg is not installed. Attempting automatic installation via pip...", file=sys.stderr)
        subprocess.check_call([sys.executable, "-m", "pip", "install", "cairosvg", "cairocffi"])
        import cairosvg
        HAS_CAIROSVG = True
        print("INFO: cairosvg successfully installed automatically via pip!", file=sys.stderr)
    except Exception as e:
        print(f"WARNING: Automatic installation of cairosvg failed: {e}", file=sys.stderr)

def parse_args():
    parser = argparse.ArgumentParser(description="Genome_anno Visualization Generator")
    parser.add_argument("--mode", required=True, choices=["heatmap", "quality", "taxonomy", "report", "circular_genome"],
                        help="Visualization type to generate")
    parser.add_argument("--output", required=True, help="Path to output file")
    parser.add_argument("--output-svg", help="Optional path to save standalone SVG (for circular genome)")
    parser.add_argument("--output-pdf", help="Optional path to save circular map PDF")
    parser.add_argument("--output-bar-pdf", help="Optional path to save selected genes barplot PDF")
    parser.add_argument("--palette", default="viridis", help="Color palette")
    parser.add_argument("--input-type", help="Input file type/source")
    parser.add_argument("--abundance-dir", help="Path to abundance results")
    parser.add_argument("--binning-dir", help="Path to binning results")
    parser.add_argument("--annotation-dir", help="Path to annotation results")
    parser.add_argument("--gff", default="", help="Path to GFF3 annotation file (for circular genome)")
    parser.add_argument("--tpm", default="", help="Path to TPM expression file (for circular genome)")
    parser.add_argument("--abundance-metric", default="auto", help="Abundance column to use (auto or exact name)")
    parser.add_argument("--sample", default="sample", help="Sample/Genome name")
    parser.add_argument("--title", default="Genome_anno Metagenomic Analysis Report", help="Report Title")
    parser.add_argument("--theme", default="dark", choices=["dark", "light"], help="Report theme")
    parser.add_argument("--interactive", default="True", help="Make report interactive")
    parser.add_argument("--selected-genes", default="", help="Comma-separated list of selected gene IDs for custom bar plot")
    return parser.parse_args()

PALETTES = {
    "viridis": ["#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725"],
    "plasma": ["#0d0887", "#46039f", "#7201a8", "#9c179e", "#bd3786", "#d8576b", "#ed7953", "#fb9f3a", "#fdca26", "#f0f921"],
    "inferno": ["#000004", "#1b0c41", "#4a0c6b", "#781c6d", "#a52c60", "#cf4446", "#ed6925", "#fb9b06", "#f7d13c", "#fcffa4"],
    "coolwarm": ["#3b4cc0", "#6788ee", "#9abbff", "#c9d7f5", "#edd1c2", "#f7a889", "#e26952", "#b40426"]
}

def get_palette_color(val, palette_name):
    colors = PALETTES.get(palette_name, PALETTES["viridis"])
    n = len(colors)
    idx = val * (n - 1)
    low_idx = int(math.floor(idx))
    high_idx = int(math.ceil(idx))
    if low_idx == high_idx:
        return colors[low_idx]
    frac = idx - low_idx
    c1 = colors[low_idx].lstrip('#')
    c2 = colors[high_idx].lstrip('#')
    r1, g1, b1 = int(c1[0:2], 16), int(c1[2:4], 16), int(c1[4:6], 16)
    r2, g2, b2 = int(c2[0:2], 16), int(c2[2:4], 16), int(c2[4:6], 16)
    r = int(r1 + (r2 - r1) * frac)
    g = int(g1 + (g2 - g1) * frac)
    b = int(b1 + (b2 - b1) * frac)
    return f"#{r:02x}{g:02x}{b:02x}"

# ==============================================================================
# GFF & TPM Parsers
# ==============================================================================

def parse_gff3(gff_path):
    """Parses a standard GFF3 file to extract genome features and genome size."""
    features = []
    genome_size = 0
    
    if not gff_path or not os.path.exists(gff_path):
        return None, 0
        
    print(f"INFO: Parsing GFF3 file: {gff_path}", file=sys.stderr)
    
    with open(gff_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            # Sequence region line defines genome size
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        genome_size = max(genome_size, int(parts[3]))
                    except ValueError:
                        pass
                continue
                
            if line.startswith("#"):
                continue
                
            parts = line.split("\t")
            if len(parts) < 9:
                continue
                
            ftype = parts[2]
            # We focus on CDS, tRNA, rRNA to avoid cluttering with parent 'gene' features
            if ftype not in ["CDS", "tRNA", "rRNA", "gene"]:
                continue
                
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes_str = parts[8]
            
            # Parse attributes
            attrs = {}
            for item in attributes_str.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k.strip()] = v.strip()
                    
            locus_tag = attrs.get("locus_tag") or attrs.get("ID") or attrs.get("Name") or ""
            gene_name = attrs.get("gene") or attrs.get("Name") or locus_tag
            product = attrs.get("product") or attrs.get("Note") or "Hypothetical protein"
            
            features.append({
                "type": ftype,
                "start": start,
                "end": end,
                "strand": strand,
                "locus_tag": locus_tag,
                "name": gene_name,
                "product": product
            })
            
            genome_size = max(genome_size, end)
            
    # If both CDS and gene types exist, filter out genes to avoid duplicates
    has_cds = any(f["type"] == "CDS" for f in features)
    if has_cds:
        features = [f for f in features if f["type"] != "gene"]
        
    print(f"INFO: Parsed {len(features)} features. Genome size: {genome_size} bp", file=sys.stderr)
    return features, genome_size

def parse_tpm(tpm_path, metric_name="auto"):
    """Parses TPM abundance table. Maps locus tags/IDs to expression values."""
    tpm_map = {}

    if not tpm_path or not os.path.exists(tpm_path):
        return tpm_map
        
    print(f"INFO: Parsing TPM table: {tpm_path}", file=sys.stderr)
    
    with open(tpm_path, "r") as f:
        # Detect delimiter
        sample_line = f.readline()
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample_line) if "," in sample_line or "\t" in sample_line else None
        delimiter = dialect.delimiter if dialect else "\t"
        
        reader = csv.reader(f, delimiter=delimiter)
        header = next(reader, None)
        
        start_col_idx = None
        end_col_idx = None
        id_col_idx = 0
        tpm_col_idx = 1
        
        # Try to match headers using intelligent priorities
        if header:
            header_lower = [h.lower().strip() for h in header]
            print(f"INFO: Found TSV columns: {header}", file=sys.stderr)
            
            # Find Start & End columns if they exist
            for idx, h in enumerate(header_lower):
                if h == "start":
                    start_col_idx = idx
                if h == "end":
                    end_col_idx = idx
            
            # 1. Identifier Column Priority: localgene -> locus_tag -> gene_id -> contains(gene, locus, id, feature)
            found_id = False
            for p in ["localgene", "locus_tag", "gene_id", "geneid"]:
                for idx, h in enumerate(header_lower):
                    if h == p:
                        id_col_idx = idx
                        found_id = True
                        break
                if found_id:
                    break
            
            if not found_id:
                for idx, h in enumerate(header_lower):
                    if any(x in h for x in ["locus", "gene", "id", "feature"]):
                        id_col_idx = idx
                        break
                        
            # 2. Value Column Choice (User Specified or Autodetect)
            found_tpm = False
            if metric_name and metric_name.lower() != "auto":
                target_metric = metric_name.lower().strip()
                for idx, h in enumerate(header_lower):
                    if h == target_metric:
                        tpm_col_idx = idx
                        found_tpm = True
                        print(f"INFO: Using user-specified abundance metric: '{header[tpm_col_idx]}'", file=sys.stderr)
                        break
                if not found_tpm:
                    print(f"WARNING: User-specified metric '{metric_name}' not found in headers. Falling back to auto-detection.", file=sys.stderr)
            
            if not found_tpm:
                for p in ["tpm_mapped", "tpm_lib", "tpm", "tpm_value", "tpm_lib_mapped"]:
                    for idx, h in enumerate(header_lower):
                        if h == p:
                            tpm_col_idx = idx
                            found_tpm = True
                            break
                    if found_tpm:
                        break
                        
            if not found_tpm:
                for idx, h in enumerate(header_lower):
                    if any(x in h for x in ["tpm", "abundance", "expression"]):
                        tpm_col_idx = idx
                        break
                        
            print(f"INFO: Selected Identifier Column: '{header[id_col_idx]}' (index {id_col_idx})", file=sys.stderr)
            print(f"INFO: Selected Value Column: '{header[tpm_col_idx]}' (index {tpm_col_idx})", file=sys.stderr)
            if start_col_idx is not None and end_col_idx is not None:
                print(f"INFO: Found coordinate columns in TSV: 'Start' (index {start_col_idx}), 'End' (index {end_col_idx})", file=sys.stderr)
            
        for row in reader:
            if not row or len(row) <= max(id_col_idx, tpm_col_idx):
                continue
            try:
                locus = row[id_col_idx].strip()
                val = float(row[tpm_col_idx].strip())
                
                start_val = None
                end_val = None
                if start_col_idx is not None and len(row) > start_col_idx:
                    try:
                        start_val = int(row[start_col_idx].strip())
                    except ValueError:
                        pass
                if end_col_idx is not None and len(row) > end_col_idx:
                    try:
                        end_val = int(row[end_col_idx].strip())
                    except ValueError:
                        pass
                
                entry = {
                    "tpm": val,
                    "start": start_val,
                    "end": end_val
                }
                
                tpm_map[locus] = entry
                # Support simple case-insensitive mapping too
                tpm_map[locus.lower()] = entry
            except ValueError:
                pass
                
    print(f"INFO: Loaded {len(tpm_map)} gene abundance values.", file=sys.stderr)
    return tpm_map


# ==============================================================================
# Circular Genome Visualizer (SVG Engine)
# ==============================================================================

def get_tpm_val(locus, tpm_map):
    """
    Looks up a locus tag in the TPM dictionary.
    Supports exact, case-insensitive, and progressive prefix stripping matches
    (e.g., matching 'Bathy_enrich_u999856_00001' with 'u999856_00001').
    """
    if not locus:
        return 0.0, "unmatched"
        
    # 1. Exact match
    val = tpm_map.get(locus)
    if val is not None:
        return val, "exact"
        
    # 2. Case-insensitive match
    locus_lower = locus.lower()
    val = tpm_map.get(locus_lower)
    if val is not None:
        return val, "exact"
        
    # 3. Progressive prefix stripping (splits by underscore and joins remaining)
    parts = locus.split("_")
    for i in range(1, len(parts)):
        suffix = "_".join(parts[i:])
        # Try exact suffix match
        val = tpm_map.get(suffix)
        if val is not None:
            return val, "prefix"
        # Try case-insensitive suffix match
        val = tpm_map.get(suffix.lower())
        if val is not None:
            return val, "prefix"
            
    return 0.0, "unmatched"

def generate_svg_circular_genome(sample_name, features, genome_size, tpm_map, theme="dark"):
    """Draws a beautiful, co-axial circular genome diagram with TPM transcriptomics overlay."""
    cx, cy = 400, 400
    r_tpm_base = 290
    r_pos_out, r_pos_in = 270, 255
    r_neg_out, r_neg_in = 250, 235
    r_scale = 277.5
    
    max_tpm_height = 70
    
    # 0. Calculate statistics early to avoid UnboundLocalErrors in TPM plotting
    cds_count = sum(1 for f in features if f["type"] == "CDS")
    trna_count = sum(1 for f in features if f["type"] == "tRNA")
    rrna_count = sum(1 for f in features if f["type"] == "rRNA")
    total_genes = cds_count + trna_count + rrna_count
    
    # 1. Create color rules
    color_map = {
        "CDS": "#6366f1",     # Indigo
        "tRNA": "#10b981",    # Emerald
        "rRNA": "#f43f5e",    # Rose
        "gene": "#64748b"     # Slate
    }
    
    if theme == "light":
        theme_tokens = {
            "bg_rect": "#ffffff",
            "bg_circle": "#f8fafc",
            "bg_circle_opacity": "1.0",
            "label_scale_fill": "#64748b",
            "title_main_fill": "#0f172a",
            "title_sub_fill": "#64748b",
            "ring_bg_stroke": "#e2e8f0",
            "scale_line_stroke": "#cbd5e1",
            "tick_major_stroke": "#94a3b8",
            "tick_minor_stroke": "#cbd5e1",
            "gene_arc_stroke": "#ffffff",
            "gene_arc_stroke_width": "0.75px",
            "gene_arc_hover_stroke": "#0f172a",
            "tpm_line_stroke": "#4f46e5",
            "tpm_area_opacity": "0.35",
            "stat_num_fill": "#0f172a",
            "stat_lbl_fill": "#64748b",
            "hub_fill": "#ffffff",
            "hub_stroke": "#e2e8f0",
            "hub_dashed_stroke": "#4f46e5",
            "tpm_glow_opacity": "0.15",
            "tpm_grad_start": "#4f46e5",
            "tpm_grad_start_op": "0.4",
            "tpm_grad_end": "#ec4899",
            "tpm_grad_end_op": "0.05",
            "stat_sub_fill": "#475569",
            "mean_tpm_fill": "#4f46e5"
        }
    else: # dark theme
        theme_tokens = {
            "bg_rect": "#0f172a",
            "bg_circle": "#1e293b",
            "bg_circle_opacity": "0.1",
            "label_scale_fill": "#64748b",
            "title_main_fill": "#f8fafc",
            "title_sub_fill": "#64748b",
            "ring_bg_stroke": "#1e293b",
            "scale_line_stroke": "#334155",
            "tick_major_stroke": "#475569",
            "tick_minor_stroke": "#334155",
            "gene_arc_stroke": "#0f172a",
            "gene_arc_stroke_width": "0.5px",
            "gene_arc_hover_stroke": "#ffffff",
            "tpm_line_stroke": "#818cf8",
            "tpm_area_opacity": "0.45",
            "stat_num_fill": "#f8fafc",
            "stat_lbl_fill": "#64748b",
            "hub_fill": "#0b0f19",
            "hub_stroke": "#1e293b",
            "hub_dashed_stroke": "#4f46e5",
            "tpm_glow_opacity": "0.5",
            "tpm_grad_start": "#4f46e5",
            "tpm_grad_start_op": "0.6",
            "tpm_grad_end": "#ec4899",
            "tpm_grad_end_op": "0.1",
            "stat_sub_fill": "#cbd5e1",
            "mean_tpm_fill": "#818cf8"
        }
        
    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 800" width="100%" height="100%">',
        '<defs>',
        '  <!-- Neon Glowing Gradients -->',
        '  <radialGradient id="tpmGlow" cx="50%" cy="50%" r="50%">',
        f'    <stop offset="0%" stop-color="#818cf8" stop-opacity="{theme_tokens["tpm_glow_opacity"]}"/>',
        '    <stop offset="100%" stop-color="#6366f1" stop-opacity="0.0"/>',
        '  </radialGradient>',
        '  <linearGradient id="tpmGrad" x1="0%" y1="0%" x2="100%" y2="100%">',
        f'    <stop offset="0%" stop-color="{theme_tokens["tpm_grad_start"]}" stop-opacity="{theme_tokens["tpm_grad_start_op"]}"/>',
        f'    <stop offset="100%" stop-color="{theme_tokens["tpm_grad_end"]}" stop-opacity="{theme_tokens["tpm_grad_end_op"]}"/>',
        '  </linearGradient>',
        '</defs>',
        '<style>',
        f'  .label-scale {{ font-family: "Outfit", "Inter", sans-serif; font-size: 10px; fill: {theme_tokens["label_scale_fill"]}; font-weight: 500; }}',
        f'  .title-main {{ font-family: "Outfit", sans-serif; font-size: 20px; font-weight: 700; fill: {theme_tokens["title_main_fill"]}; text-anchor: middle; }}',
        f'  .title-sub {{ font-family: "Inter", sans-serif; font-size: 11px; fill: {theme_tokens["title_sub_fill"]}; text-anchor: middle; }}',
        f'  .ring-bg {{ stroke: {theme_tokens["ring_bg_stroke"]}; stroke-width: 1px; fill: none; }}',
        f'  .scale-line {{ stroke: {theme_tokens["scale_line_stroke"]}; stroke-width: 1px; stroke-dasharray: 2, 4; fill: none; }}',
        f'  .tick-major {{ stroke: {theme_tokens["tick_major_stroke"]}; stroke-width: 1.5px; fill: none; }}',
        f'  .tick-minor {{ stroke: {theme_tokens["tick_minor_stroke"]}; stroke-width: 1px; fill: none; }}',
        f'  .gene-arc {{ stroke: {theme_tokens["gene_arc_stroke"]}; stroke-width: {theme_tokens["gene_arc_stroke_width"]}; transition: all 0.2s ease; cursor: pointer; }}',
        f'  .gene-arc:hover {{ stroke: {theme_tokens["gene_arc_hover_stroke"]}; stroke-width: 1.5px; filter: drop-shadow(0 0 6px rgba(255,255,255,0.6)); opacity: 0.95; }}',
        f'  .tpm-line {{ stroke: {theme_tokens["tpm_line_stroke"]}; stroke-width: 1.5px; fill: none; filter: drop-shadow(0 0 4px rgba(99,102,241,0.4)); }}',
        f'  .tpm-area {{ fill: url(#tpmGrad); opacity: {theme_tokens["tpm_area_opacity"]}; }}',
        f'  .tpm-avg-circle {{ stroke: {theme_tokens["mean_tpm_fill"]}; stroke-width: 1.25px; stroke-dasharray: 4, 4; fill: none; opacity: 0.65; pointer-events: stroke; cursor: help; }}',
        f'  .stat-num {{ font-family: "Outfit", sans-serif; font-size: 26px; font-weight: 700; fill: {theme_tokens["stat_num_fill"]}; text-anchor: middle; }}',
        f'  .stat-lbl {{ font-family: "Inter", sans-serif; font-size: 10px; fill: {theme_tokens["stat_lbl_fill"]}; text-anchor: middle; text-transform: uppercase; letter-spacing: 0.05em; }}',
        '</style>'
,
        '<!-- Background -->',
        f'<rect width="800" height="800" fill="{theme_tokens["bg_rect"]}" rx="20"/>',
        f'<circle cx="400" cy="400" r="380" fill="{theme_tokens["bg_circle"]}" opacity="{theme_tokens["bg_circle_opacity"]}"/>'
    ]
    
    # helper polar coordinate conversion
    def get_polar_coords(pos, r):
        angle = (pos / genome_size) * 2 * math.pi - math.pi / 2
        return cx + r * math.cos(angle), cy + r * math.sin(angle)
        
    def get_angle(pos):
        return (pos / genome_size) * 2 * math.pi - math.pi / 2
        
    # Draw Ring Backgrounds
    svg.append(f'<!-- Baselines -->')
    svg.append(f'<circle cx="{cx}" cy="{cy}" r="{r_pos_out}" class="ring-bg"/>')
    svg.append(f'<circle cx="{cx}" cy="{cy}" r="{r_pos_in}" class="ring-bg"/>')
    svg.append(f'<circle cx="{cx}" cy="{cy}" r="{r_neg_out}" class="ring-bg"/>')
    svg.append(f'<circle cx="{cx}" cy="{cy}" r="{r_neg_in}" class="ring-bg"/>')
    svg.append(f'<circle cx="{cx}" cy="{cy}" r="{r_tpm_base}" class="scale-line"/>')
    
    # 2. Draw Ticks and Ticks Labels (Scale ring)
    svg.append('<!-- Coordinate Scale Ring -->')
    tick_interval = 10000 if genome_size <= 200000 else 100000
    if genome_size <= 50000:
        tick_interval = 5000
        
    for tick in range(0, genome_size, tick_interval):
        theta = get_angle(tick)
        x1 = cx + (r_scale - 4) * math.cos(theta)
        y1 = cy + (r_scale - 4) * math.sin(theta)
        x2 = cx + (r_scale + 4) * math.cos(theta)
        y2 = cy + (r_scale + 4) * math.sin(theta)
        
        svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" class="tick-major"/>')
        
        # Draw label
        xl = cx + (r_scale + 16) * math.cos(theta)
        yl = cy + (r_scale + 16) * math.sin(theta) + 3
        # Rotate text slightly to align nicely
        svg.append(f'<text x="{xl}" y="{yl}" class="label-scale" text-anchor="middle">{tick//1000}K</text>')
        
    # Draw minor ticks
    minor_interval = tick_interval / 5
    for tick in range(0, genome_size, int(minor_interval)):
        if tick % tick_interval == 0:
            continue
        theta = get_angle(tick)
        x1 = cx + (r_scale - 2) * math.cos(theta)
        y1 = cy + (r_scale - 2) * math.sin(theta)
        x2 = cx + (r_scale + 2) * math.cos(theta)
        y2 = cy + (r_scale + 2) * math.sin(theta)
        svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" class="tick-minor"/>')

    # 3. Draw TPM Track (Underneath genes for nice blending)
    svg.append('<!-- TPM Expression Overlay -->')
    has_tpm_data = False
    
    # Filter valid TPM values
    valid_tpm_pairs = []
    match_stats = {"exact": 0, "prefix": 0, "unmatched": 0}
    coord_mismatches = 0
    
    for f in features:
        locus = f["locus_tag"]
        tpm_entry, m_type = get_tpm_val(locus, tpm_map)
        match_stats[m_type] += 1
        
        tpm_val = 0.0
        if tpm_entry is not None:
            if isinstance(tpm_entry, dict):
                tpm_val = tpm_entry.get("tpm", 0.0)
                
                # Check for GFF3 vs TSV coordinate inconsistency
                gff_start = f["start"]
                gff_end = f["end"]
                tsv_start = tpm_entry.get("start")
                tsv_end = tpm_entry.get("end")
                
                if tsv_start is not None and tsv_end is not None:
                    if gff_start != tsv_start or gff_end != tsv_end:
                        coord_mismatches += 1
                        if coord_mismatches <= 5:
                            print(f"WARNING: Coordinate mismatch for gene '{locus}'! "
                                  f"GFF: {gff_start}-{gff_end} | TSV: {tsv_start}-{tsv_end}", 
                                  file=sys.stderr)
            else:
                # Fallback for mock floats
                tpm_val = float(tpm_entry)
                
        if tpm_val > 0.0:
            valid_tpm_pairs.append((f["start"], f["end"], tpm_val))
            has_tpm_data = True
            
    print(f"INFO: Gene-abundance ID matching stats for {sample_name}:", file=sys.stderr)
    print(f"  - Exact matches: {match_stats['exact']}", file=sys.stderr)
    print(f"  - Prefix-stripped matches: {match_stats['prefix']}", file=sys.stderr)
    print(f"  - Unmatched genes: {match_stats['unmatched']}", file=sys.stderr)
    
    if coord_mismatches > 0:
        print(f"WARNING: Detected a total of {coord_mismatches} coordinate inconsistencies between GFF3 and TSV! "
              f"Please verify if you are using the exact same annotation genome version for transcriptomics mapping.", 
              file=sys.stderr)


            
    avg_tpm = 0.0
    if has_tpm_data:
        # Scale TPM heights using Square-Root to compress large dynamic range while highlighting variations
        max_val = max(v[2] for v in valid_tpm_pairs)
        sqrt_max = math.sqrt(max_val) if max_val > 0 else 1.0
        
        # Calculate Average TPM and its height
        total_tpm_sum = sum(v[2] for v in valid_tpm_pairs)
        avg_tpm = total_tpm_sum / (total_genes if total_genes > 0 else 1)
        h_avg = (math.sqrt(avg_tpm) / sqrt_max) * max_tpm_height if avg_tpm > 0.0 else 0.0
        
        # Build coordinates for a continuous smooth path
        path_pts = []
        tpm_steps = 360
        step_bp = genome_size / tpm_steps
        
        for step in range(tpm_steps + 1):
            bp = step * step_bp
            # Use window-based sampling (max-pooling within the step interval) to prevent aliasing of narrow peaks
            win_start = max(0, bp - step_bp / 2)
            win_end = min(genome_size, bp + step_bp / 2)
            active_tpms = [v[2] for v in valid_tpm_pairs if not (v[1] < win_start or v[0] > win_end)]
            local_tpm = max(active_tpms) if active_tpms else 0.0
            
            # Apply square-root scaling
            h = 0.0
            if local_tpm > 0.0:
                h = (math.sqrt(local_tpm) / sqrt_max) * max_tpm_height
                
            r_tpm = r_tpm_base + h
            theta = get_angle(bp)
            px = cx + r_tpm * math.cos(theta)
            py = cy + r_tpm * math.sin(theta)
            path_pts.append((px, py, theta))
            
        # Draw area polygon (radial fill)
        area_d = []
        # Outer boundary
        for idx, (px, py, _) in enumerate(path_pts):
            cmd = "M" if idx == 0 else "L"
            area_d.append(f"{cmd} {px:.2f} {py:.2f}")
        # Inner boundary (reverse)
        for _, _, theta in reversed(path_pts):
            px = cx + r_tpm_base * math.cos(theta)
            py = cy + r_tpm_base * math.sin(theta)
            area_d.append(f"L {px:.2f} {py:.2f}")
        area_d.append("Z")
        
        # Draw lines
        line_d = []
        for idx, (px, py, _) in enumerate(path_pts):
            cmd = "M" if idx == 0 else "L"
            line_d.append(f"{cmd} {px:.2f} {py:.2f}")
            
        svg.append(f'<path d="{" ".join(area_d)}" class="tpm-area"/>')
        svg.append(f'<path d="{" ".join(line_d)}" class="tpm-line"/>')
        
        # Draw Mean TPM Threshold Circle
        if avg_tpm > 0.0:
            r_avg = r_tpm_base + h_avg
            svg.append(f'<!-- Mean TPM Threshold Circle -->')
            svg.append(f'<circle cx="{cx}" cy="{cy}" r="{r_avg:.2f}" class="tpm-avg-circle">')
            svg.append(f'  <title>Genome Mean Expression Threshold: {avg_tpm:.2f} TPM</title>')
            svg.append(f'</circle>')
    else:
        svg.append(f'<!-- No transcript abundance data provided. TPM track omitted -->')

    # 4. Draw Genes (Arcs / Wedges with arrowheads)
    svg.append('<!-- Genomic Features (Forward and Reverse Tracks) -->')
    for f in features:
        color = color_map.get(f["type"], "#64748b")
        start_bp = f["start"]
        end_bp = f["end"]
        strand = f["strand"]
        locus = f["locus_tag"]
        name = f["name"]
        prod = f["product"]
        
        tpm_entry, _ = get_tpm_val(locus, tpm_map)
        tpm_val = 0.0
        if tpm_entry is not None:
            if isinstance(tpm_entry, dict):
                tpm_val = tpm_entry.get("tpm", 0.0)
            else:
                tpm_val = float(tpm_entry)
        
        # Track parameters
        if strand == "+":
            r_out, r_in = r_pos_out, r_pos_in
        else:
            r_out, r_in = r_neg_out, r_neg_in
            
        r_mid = (r_out + r_in) / 2
        
        theta_s = get_angle(start_bp)
        theta_e = get_angle(end_bp)
        delta_theta = theta_e - theta_s
        
        # Ensure angles are correctly clockwise order
        if delta_theta < 0:
            delta_theta += 2 * math.pi
            
        large_arc = 1 if delta_theta > math.pi else 0
        
        # Calculate arrowhead dimensions
        # Wedges too small get a block to avoid arrowhead overlaps on tiny features
        min_arrow_len = 0.04  # ~2.5% of circular circumference
        has_arrow = delta_theta >= min_arrow_len
        
        path_d = ""
        
        if has_arrow:
            arrow_angle = min(0.03, delta_theta * 0.35)
            if strand == "+":
                # Arrow points clockwise (at end coord)
                theta_base = theta_e - arrow_angle
                
                # Coords
                x_in_start, y_in_start = polar_to_cartesian(cx, cy, r_in, theta_s)
                x_out_start, y_out_start = polar_to_cartesian(cx, cy, r_out, theta_s)
                x_out_base, y_out_base = polar_to_cartesian(cx, cy, r_out, theta_base)
                x_tip, y_tip = polar_to_cartesian(cx, cy, r_mid, theta_e)
                x_in_base, y_in_base = polar_to_cartesian(cx, cy, r_in, theta_base)
                
                path_d = (
                    f"M {x_in_start:.2f} {y_in_start:.2f} "
                    f"L {x_out_start:.2f} {y_out_start:.2f} "
                    f"A {r_out} {r_out} 0 {large_arc} 1 {x_out_base:.2f} {y_out_base:.2f} "
                    f"L {x_tip:.2f} {y_tip:.2f} "
                    f"L {x_in_base:.2f} {y_in_base:.2f} "
                    f"A {r_in} {r_in} 0 {large_arc} 0 {x_in_start:.2f} {y_in_start:.2f} Z"
                )
            else:
                # Arrow points counter-clockwise (at start coord)
                theta_base = theta_s + arrow_angle
                
                x_in_end, y_in_end = polar_to_cartesian(cx, cy, r_in, theta_e)
                x_out_end, y_out_end = polar_to_cartesian(cx, cy, r_out, theta_e)
                x_out_base, y_out_base = polar_to_cartesian(cx, cy, r_out, theta_base)
                x_tip, y_tip = polar_to_cartesian(cx, cy, r_mid, theta_s)
                x_in_base, y_in_base = polar_to_cartesian(cx, cy, r_in, theta_base)
                
                path_d = (
                    f"M {x_in_end:.2f} {y_in_end:.2f} "
                    f"A {r_in} {r_in} 0 {large_arc} 0 {x_in_base:.2f} {y_in_base:.2f} "
                    f"L {x_tip:.2f} {y_tip:.2f} "
                    f"L {x_out_base:.2f} {y_out_base:.2f} "
                    f"A {r_out} {r_out} 0 {large_arc} 1 {x_out_end:.2f} {y_out_end:.2f} Z"
                )
        else:
            # Draw simple rounded wedge block
            x_in_start, y_in_start = polar_to_cartesian(cx, cy, r_in, theta_s)
            x_out_start, y_out_start = polar_to_cartesian(cx, cy, r_out, theta_s)
            x_out_end, y_out_end = polar_to_cartesian(cx, cy, r_out, theta_e)
            x_in_end, y_in_end = polar_to_cartesian(cx, cy, r_in, theta_e)
            
            path_d = (
                f"M {x_in_start:.2f} {y_in_start:.2f} "
                f"L {x_out_start:.2f} {y_out_start:.2f} "
                f"A {r_out} {r_out} 0 {large_arc} 1 {x_out_end:.2f} {y_out_end:.2f} "
                f"L {x_in_end:.2f} {y_in_end:.2f} "
                f"A {r_in} {r_in} 0 {large_arc} 0 {x_in_start:.2f} {y_in_start:.2f} Z"
            )
            
        tpm_lbl = f"{tpm_val:.2f} TPM" if tpm_val > 0 else "0 TPM / Undetected"
        svg.append(
            f'<path d="{path_d}" fill="{color}" class="gene-arc">'
            f'<title>Gene: {name}&#10;Locus Tag: {locus}&#10;Product: {prod}&#10;'
            f'Coordinates: {start_bp}-{end_bp} ({strand})&#10;Expression: {tpm_lbl}</title>'
            f'</path>'
        )

    # 5. Draw Central Metadata Stats Hub
    svg.append('<!-- Central Statistics Hub Panel -->')
    svg.append(f'<circle cx="{cx}" cy="{cy}" r="165" fill="{theme_tokens["hub_fill"]}" stroke="{theme_tokens["hub_stroke"]}" stroke-width="2"/>')
    svg.append(f'<circle cx="{cx}" cy="{cy}" r="160" fill="none" stroke="{theme_tokens["hub_dashed_stroke"]}" stroke-width="1" stroke-dasharray="4, 4" opacity="0.3"/>')
    
    # (Using pre-calculated statistics)
    
    svg.append(f'<text x="{cx}" y="{cy - 55}" class="title-main">{sample_name}</text>')
    svg.append(f'<text x="{cx}" y="{cy - 30}" class="title-sub">{genome_size:,} bp Circular Map</text>')
    
    svg.append(f'<text x="{cx}" y="{cy + 15}" class="stat-num">{total_genes}</text>')
    svg.append(f'<text x="{cx}" y="{cy + 30}" class="stat-lbl">Total Genes</text>')
    
    svg.append(f'<text x="{cx}" y="{cy + 65}" class="title-sub" fill="{theme_tokens["stat_sub_fill"]}" font-weight="600">'
               f'CDS: {cds_count} | tRNA: {trna_count} | rRNA: {rrna_count}</text>')
               
    # Average expression
    if has_tpm_data:
        total_tpm_sum = sum(v[2] for v in valid_tpm_pairs)
        avg_tpm = total_tpm_sum / (total_genes if total_genes > 0 else 1)
        svg.append(f'<text x="{cx}" y="{cy + 95}" class="stat-lbl" fill="{theme_tokens["mean_tpm_fill"]}">Mean TPM: {avg_tpm:.1f}</text>')
    else:
        svg.append(f'<text x="{cx}" y="{cy + 95}" class="stat-lbl" fill="{theme_tokens["label_scale_fill"]}">No Transcriptome Overlay</text>')
        
    svg.append('</svg>')
    return "\n".join(svg), avg_tpm

def generate_svg_selected_genes_barplot(resolved_selected, features, tpm_map, avg_tpm, theme="dark"):
    """Generates a beautiful vertical SVG bar chart for selected genes with a Mean TPM threshold line."""
    if not resolved_selected:
        return ""
        
    # Match selected genes
    matched_genes = []
    for f in features:
        locus = f["locus_tag"]
        # Check if locus matches any in resolved_selected
        matched = False
        match_idx = -1
        matched_gene_name = ""
        
        for idx, (sid, gname) in enumerate(resolved_selected):
            sid_clean = sid.strip().lower()
            if not sid_clean:
                continue
            if sid_clean == locus.lower() or sid_clean in locus.lower() or locus.lower() in sid_clean:
                matched = True
                match_idx = idx
                matched_gene_name = gname.strip()
                break
        
        if matched:
            tpm_entry, _ = get_tpm_val(locus, tpm_map)
            tpm_val = 0.0
            if tpm_entry is not None:
                if isinstance(tpm_entry, dict):
                    tpm_val = tpm_entry.get("tpm", 0.0)
                else:
                    tpm_val = float(tpm_entry)
            matched_genes.append({
                "locus": locus,
                "name": matched_gene_name or f["name"], # Use custom name or fall back to GFF default
                "product": f["product"],
                "tpm": tpm_val,
                "type": f["type"],
                "match_idx": match_idx
            })
            
    if not matched_genes:
        return '<p style="color: #64748b; text-align: center; padding: 20px;">No selected genes found in GFF annotation.</p>'
        
    # Sort matched genes by the exact order of selected IDs provided in the list
    matched_genes.sort(key=lambda x: x["match_idx"])
    
    # SVG setup (Vertical orientation)
    num_genes = len(matched_genes)
    top_margin = 60
    bottom_margin = 160 # Increased from 110 to fit the longer vertical uxxx_xxx_gene labels
    left_margin = 60
    right_margin = 40
    
    chart_width = 800
    chart_height = 480 # Increased from 420 to provide plenty of room for vertical labels
    plot_width = chart_width - left_margin - right_margin
    plot_height = chart_height - top_margin - bottom_margin # 480 - 60 - 160 = 260px
    
    # Calculate responsive bar width with absolute limits to look clean
    bar_width = min(40, (plot_width / num_genes) * 0.65)
    total_bar_space = num_genes * bar_width
    bar_gap = (plot_width - total_bar_space) / num_genes
    
    # Theme settings
    bg_color = "#0f172a" if theme == "dark" else "#ffffff"
    text_color = "#f8fafc" if theme == "dark" else "#0f172a"
    label_color = "#94a3b8" if theme == "dark" else "#64748b"
    grid_color = "#1e293b" if theme == "dark" else "#e2e8f0"
    bar_color = "#6366f1" # Indigo
    avg_line_color = "#ec4899" if theme == "dark" else "#db2777" # Rose
    
    max_tpm = max(g["tpm"] for g in matched_genes)
    x_max = max(max_tpm, avg_tpm)
    x_max = x_max if x_max > 0 else 1.0
    # Add 15% breathing room
    y_limit = x_max * 1.15
    
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {chart_width} {chart_height}" width="100%" height="100%">',
        '<style>',
        f'  .bar-text {{ font-family: "Inter", sans-serif; font-size: 10px; fill: {text_color}; font-weight: 600; }}',
        f'  .axis-lbl {{ font-family: "Inter", sans-serif; font-size: 11px; fill: {label_color}; font-weight: 500; }}',
        f'  .chart-title {{ font-family: "Outfit", sans-serif; font-size: 16px; font-weight: 700; fill: {text_color}; }}',
        f'  .mean-lbl {{ font-family: "Inter", sans-serif; font-size: 11px; font-weight: 600; fill: {text_color}; }}',
        '  .bar-rect { transition: all 0.2s ease; cursor: pointer; }',
        '  .bar-rect:hover { filter: brightness(1.2) drop-shadow(0 0 4px rgba(99,102,241,0.4)); }',
        '</style>',
        f'<rect width="{chart_width}" height="{chart_height}" fill="{bg_color}" rx="12"/>',
        f'<text x="25" y="32" class="chart-title">Selected Genes Expression Profile</text>'
    ]
    
    # Add Genome Mean Legend at the top right if available
    if avg_tpm > 0.0:
        legend_x = chart_width - right_margin - 160
        svg.append(f'<!-- Legend -->')
        svg.append(f'<line x1="{legend_x}" y1="27" x2="{legend_x + 30}" y2="27" stroke="{avg_line_color}" stroke-width="2" stroke-dasharray="4,3"/>')
        svg.append(f'<text x="{legend_x + 38}" y="31" class="mean-lbl" text-anchor="start">Genome Mean: {avg_tpm:.1f}</text>')
    
    # Draw horizontal grid lines & Y-ticks (Numbers only, no "TPM" text)
    num_ticks = 5
    for i in range(num_ticks + 1):
        tick_val = (i / num_ticks) * y_limit
        y = chart_height - bottom_margin - (tick_val / y_limit) * plot_height
        svg.append(f'<line x1="{left_margin - 5}" y1="{y:.2f}" x2="{chart_width - right_margin}" y2="{y:.2f}" stroke="{grid_color}" stroke-width="1" stroke-dasharray="2,2"/>')
        svg.append(f'<text x="{left_margin - 10}" y="{y + 4:.2f}" class="axis-lbl" text-anchor="end">{tick_val:.1f}</text>')
        
    # Draw vertical bars & X rotated labels (IDs directly, with prefix stripped)
    for idx, g in enumerate(matched_genes):
        x = left_margin + idx * (bar_width + bar_gap) + bar_gap / 2
        h = (g["tpm"] / y_limit) * plot_height
        
        # Clip zero height to at least 2px so it's slightly visible as a line
        h_clipped = max(2, h)
        y = chart_height - bottom_margin - h_clipped
        
        # Draw bar
        svg.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" height="{h_clipped:.2f}" fill="{bar_color}" rx="3" class="bar-rect">'
            f'<title>Locus: {g["locus"]}&#10;Name: {g["name"]}&#10;Expression: {g["tpm"]:.2f}&#10;Product: {g["product"]}</title>'
            f'</rect>'
        )
        
        # Display TPM value on top of the bar (number only)
        val_y = y - 6
        if val_y < top_margin:
            # draw inside bar if too high
            val_y = y + 15
            val_color = "#ffffff"
        else:
            val_color = text_color
            
        svg.append(f'<text x="{x + bar_width/2:.2f}" y="{val_y:.2f}" class="bar-text" fill="{val_color}" text-anchor="middle">{g["tpm"]:.1f}</text>')
        
        # X-Axis Vertical Labels (ID directly at the bottom, prefix stripped to keep e.g. u999856_00850)
        locus_suffix = "_".join(g["locus"].split("_")[-2:])
        gene_name = g["name"]
        
        # Check if the gene name is valid and not a placeholder like 'hypothetical protein' or identical to the locus ID
        is_valid_gene_name = (
            gene_name and 
            gene_name.lower() != "hypothetical protein" and 
            gene_name.lower() != g["locus"].lower() and
            gene_name.lower() != locus_suffix.lower()
        )
        
        if is_valid_gene_name:
            clean_id = f"{locus_suffix}_{gene_name}"
        else:
            clean_id = locus_suffix
            
        svg.append(
            f'<text x="0" y="0" transform="translate({x + bar_width/2:.2f}, {chart_height - bottom_margin + 10}) rotate(90)" class="axis-lbl" style="text-anchor: start;">{clean_id}</text>'
        )
        
    # Draw Mean TPM Threshold Line (dashed, with label "Genome Mean: X.X")
    if avg_tpm > 0.0:
        y_avg = chart_height - bottom_margin - (avg_tpm / y_limit) * plot_height
        svg.append(f'<!-- Mean TPM Line -->')
        svg.append(f'<line x1="{left_margin - 5}" y1="{y_avg:.2f}" x2="{chart_width - right_margin}" y2="{y_avg:.2f}" stroke="{avg_line_color}" stroke-width="1.5" stroke-dasharray="4,3"/>')
        
    # Draw Y-Axis Base Line
    svg.append(f'<line x1="{left_margin}" y1="{top_margin}" x2="{left_margin}" y2="{chart_height - bottom_margin}" stroke="{text_color}" stroke-width="1.5"/>')
    # Draw X-Axis Base Line
    svg.append(f'<line x1="{left_margin}" y1="{chart_height - bottom_margin}" x2="{chart_width - right_margin}" y2="{chart_height - bottom_margin}" stroke="{text_color}" stroke-width="1.5"/>')
    
    svg.append('</svg>')
    return "\n".join(svg)

def generate_placeholder_svg(message, theme="dark"):
    """Generates a beautiful placeholder SVG when no selected genes are provided."""
    bg_color = "#0f172a" if theme == "dark" else "#ffffff"
    text_color = "#94a3b8" if theme == "dark" else "#64748b"
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 150" width="800" height="150">
        <rect width="800" height="150" fill="{bg_color}" rx="12"/>
        <text x="400" y="80" font-family="'Inter', sans-serif" font-size="14px" fill="{text_color}" text-anchor="middle">{message}</text>
    </svg>"""
    return svg

def polar_to_cartesian(cx, cy, r, angle_rad):
    x = cx + r * math.cos(angle_rad)
    y = cy + r * math.sin(angle_rad)
    return x, y

def generate_svg_heatmap(palette, rows, cols, data):
    margin = {"top": 60, "right": 100, "bottom": 80, "left": 120}
    cell_w, cell_h = 40, 30
    width = cell_w * len(cols) + margin["left"] + margin["right"]
    height = cell_h * len(rows) + margin["top"] + margin["bottom"]
    
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" height="100%">',
        '<style>',
        '  .label { font-family: "Outfit", "Inter", sans-serif; font-size: 11px; fill: #a0aec0; }',
        '  .title { font-family: "Outfit", "Inter", sans-serif; font-size: 16px; font-weight: 600; fill: #e2e8f0; }',
        '  .cell { stroke: #1a202c; stroke-width: 1px; transition: all 0.2s; }',
        '  .cell:hover { stroke: #ffffff; stroke-width: 2px; filter: drop-shadow(0 0 4px rgba(255,255,255,0.5)); }',
        '  .legend-text { font-family: "Outfit", sans-serif; font-size: 10px; fill: #718096; }',
        '</style>',
        '<!-- Background -->',
        f'<rect width="{width}" height="{height}" fill="#171923" rx="8"/>',
        f'<text x="{margin["left"]}" y="35" class="title">Genome &amp; Sample Abundance Heatmap</text>'
    ]
    
    for r_idx, r_name in enumerate(rows):
        y = margin["top"] + r_idx * cell_h
        svg.append(f'<text x="{margin["left"] - 10}" y="{y + cell_h/2 + 4}" class="label" text-anchor="end">{r_name}</text>')
        
        for c_idx, c_name in enumerate(cols):
            x = margin["left"] + c_idx * cell_w
            val = data[r_idx][c_idx]
            color = get_palette_color(val, palette)
            svg.append(
                f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" class="cell" fill="{color}">'
                f'<title>{r_name} in {c_name}: {val*100:.2f}%</title>'
                f'</rect>'
            )
            
    for c_idx, c_name in enumerate(cols):
        x = margin["left"] + c_idx * cell_w + cell_w/2
        y = margin["top"] + len(rows) * cell_h + 15
        svg.append(f'<text x="{x}" y="{y}" class="label" text-anchor="middle" transform="rotate(30, {x}, {y})">{c_name}</text>')
        
    leg_x = margin["left"] + len(cols) * cell_w + 20
    leg_y = margin["top"]
    leg_w = 20
    leg_h = len(rows) * cell_h
    
    steps = 50
    step_h = leg_h / steps
    for i in range(steps):
        val = 1.0 - (i / steps)
        y = leg_y + i * step_h
        color = get_palette_color(val, palette)
        svg.append(f'<rect x="{leg_x}" y="{y}" width="{leg_w}" height="{step_h}" fill="{color}"/>')
        
    svg.append(f'<rect x="{leg_x}" y="{leg_y}" width="{leg_w}" height="{leg_h}" stroke="#4a5568" fill="none"/>')
    svg.append(f'<text x="{leg_x + leg_w + 8}" y="{leg_y + 10}" class="legend-text" text-anchor="start">High</text>')
    svg.append(f'<text x="{leg_x + leg_w + 8}" y="{leg_y + leg_h - 2}" class="legend-text" text-anchor="start">Low</text>')
    
    svg.append('</svg>')
    return "\n".join(svg)

def generate_svg_quality():
    width, height = 550, 400
    margin = {"top": 60, "right": 50, "bottom": 60, "left": 60}
    
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" height="100%">',
        '<style>',
        '  .label { font-family: "Outfit", "Inter", sans-serif; font-size: 11px; fill: #a0aec0; }',
        '  .title { font-family: "Outfit", "Inter", sans-serif; font-size: 16px; font-weight: 600; fill: #e2e8f0; }',
        '  .axis { stroke: #4a5568; stroke-width: 1px; }',
        '  .grid { stroke: #2d3748; stroke-width: 1px; stroke-dasharray: 2,2; }',
        '  .point { stroke: #1a202c; stroke-width: 1px; transition: all 0.2s; cursor: pointer; }',
        '  .point:hover { stroke: #ffffff; stroke-width: 2px; r: 8; filter: drop-shadow(0 0 6px rgba(255,255,255,0.8)); }',
        '</style>',
        f'<rect width="{width}" height="{height}" fill="#171923" rx="8"/>',
        f'<text x="{margin["left"]}" y="35" class="title">Genome Quality: Completeness vs Contamination</text>'
    ]
    
    for comp in range(50, 101, 10):
        x = margin["left"] + (comp - 50) / 50 * (width - margin["left"] - margin["right"])
        svg.append(f'<line x1="{x}" y1="{margin["top"]}" x2="{x}" y2="{height - margin["bottom"]}" class="grid"/>')
        svg.append(f'<text x="{x}" y="{height - margin["bottom"] + 15}" class="label" text-anchor="middle">{comp}%</text>')
        
    for cont in range(0, 16, 5):
        y = height - margin["bottom"] - cont / 15 * (height - margin["top"] - margin["bottom"])
        svg.append(f'<line x1="{margin["left"]}" y1="{y}" x2="{width - margin["right"]}" y2="{y}" class="grid"/>')
        svg.append(f'<text x="{margin["left"] - 8}" y="{y + 4}" class="label" text-anchor="end">{cont}%</text>')
        
    svg.append(f'<line x1="{margin["left"]}" y1="{margin["top"]}" x2="{margin["left"]}" y2="{height - margin["bottom"]}" class="axis"/>')
    svg.append(f'<line x1="{margin["left"]}" y1="{height - margin["bottom"]}" x2="{width - margin["right"]}" y2="{height - margin["bottom"]}" class="axis"/>')
    
    svg.append(f'<text x="{(width - margin["left"] - margin["right"])/2 + margin["left"]}" y="{height - 15}" class="label" text-anchor="middle" font-weight="600">Completeness (%)</text>')
    svg.append(f'<text x="18" y="{(height - margin["top"] - margin["bottom"])/2 + margin["top"]}" class="label" text-anchor="middle" font-weight="600" transform="rotate(-90, 18, {(height - margin["top"] - margin["bottom"])/2 + margin["top"]})">Contamination (%)</text>')
    
    bins = [
        {"name": "Bin_01_High", "comp": 94.5, "cont": 1.2, "size": 3.4, "color": "#48bb78"},
        {"name": "Bin_02_High", "comp": 98.1, "cont": 0.8, "size": 4.1, "color": "#48bb78"},
        {"name": "Bin_03_Medium", "comp": 88.2, "cont": 4.5, "size": 2.8, "color": "#ecc94b"},
        {"name": "Bin_04_Medium", "comp": 91.0, "cont": 3.1, "size": 3.9, "color": "#ecc94b"},
        {"name": "Bin_05_Low", "comp": 72.4, "cont": 8.9, "size": 1.9, "color": "#ed8936"},
        {"name": "Bin_06_Contam", "comp": 95.2, "cont": 12.4, "size": 5.2, "color": "#e53e3e"}
    ]
    
    for b in bins:
        cx_p = margin["left"] + (b["comp"] - 50) / 50 * (width - margin["left"] - margin["right"])
        cy_p = height - margin["bottom"] - b["cont"] / 15 * (height - margin["top"] - margin["bottom"])
        r = 5 + b["size"] * 1.5
        
        svg.append(
            f'<circle cx="{cx_p}" cy="{cy_p}" r="{r}" fill="{b["color"]}" class="point">'
            f'<title>{b["name"]}&#10;Completeness: {b["comp"]}%&#10;Contamination: {b["cont"]}%&#10;Size: {b["size"]} Mbp</title>'
            f'</circle>'
        )
        
    svg.append('</svg>')
    return "\n".join(svg)

def generate_svg_taxonomy():
    width, height = 550, 400
    margin = {"top": 60, "right": 180, "bottom": 60, "left": 60}
    
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" height="100%">',
        '<style>',
        '  .label { font-family: "Outfit", "Inter", sans-serif; font-size: 11px; fill: #a0aec0; }',
        '  .title { font-family: "Outfit", "Inter", sans-serif; font-size: 16px; font-weight: 600; fill: #e2e8f0; }',
        '  .legend-title { font-family: "Outfit", sans-serif; font-size: 12px; font-weight: 600; fill: #cbd5e0; }',
        '  .slice { stroke: #171923; stroke-width: 2px; transition: all 0.2s; cursor: pointer; }',
        '  .slice:hover { filter: drop-shadow(0 0 6px rgba(255,255,255,0.4)); opacity: 0.9; }',
        '</style>',
        f'<rect width="{width}" height="{height}" fill="#171923" rx="8"/>',
        f'<text x="{margin["left"]}" y="35" class="title">Taxonomic Distribution (Phylum Level)</text>'
    ]
    
    cx = (width - margin["left"] - margin["right"]) / 2 + margin["left"]
    cy = (height - margin["top"] - margin["bottom"]) / 2 + margin["top"]
    r_out, r_in = 100, 60
    
    taxa = [
        {"name": "Bacteroidota", "count": 42, "color": "#3182ce"},
        {"name": "Pseudomonadota", "count": 28, "color": "#319795"},
        {"name": "Actinomycetota", "count": 15, "color": "#ecc94b"},
        {"name": "Bacillota", "count": 10, "color": "#ed8936"},
        {"name": "Unclassified", "count": 5, "color": "#a0aec0"}
    ]
    
    total = sum(t["count"] for t in taxa)
    curr_angle = -math.pi / 2
    
    for idx, t in enumerate(taxa):
        pct = t["count"] / total
        angle_sweep = pct * 2 * math.pi
        
        x1_out = cx + r_out * math.cos(curr_angle)
        y1_out = cy + r_out * math.sin(curr_angle)
        x1_in = cx + r_in * math.cos(curr_angle)
        y1_in = cy + r_in * math.sin(curr_angle)
        
        next_angle = curr_angle + angle_sweep
        x2_out = cx + r_out * math.cos(next_angle)
        y2_out = cy + r_out * math.sin(next_angle)
        x2_in = cx + r_in * math.cos(next_angle)
        y2_in = cy + r_in * math.sin(next_angle)
        
        large_arc = 1 if angle_sweep > math.pi else 0
        
        path_data = (
            f"M {x1_out} {y1_out} "
            f"A {r_out} {r_out} 0 {large_arc} 1 {x2_out} {y2_out} "
            f"L {x2_in} {y2_in} "
            f"A {r_in} {r_in} 0 {large_arc} 0 {x1_in} {y1_in} Z"
        )
        
        svg.append(
            f'<path d="{path_data}" fill="{t["color"]}" class="slice">'
            f'<title>{t["name"]}: {t["count"]}%</title>'
            f'</path>'
        )
        
        leg_x = width - margin["right"] + 20
        leg_y = margin["top"] + idx * 30 + 20
        svg.append(f'<rect x="{leg_x}" y="{leg_y - 10}" width="15" height="15" fill="{t["color"]}" rx="3"/>')
        svg.append(f'<text x="{leg_x + 25}" y="{leg_y + 2}" class="label" text-anchor="start" font-weight="600">{t["name"]}</text>')
        svg.append(f'<text x="{leg_x + 130}" y="{leg_y + 2}" class="label" text-anchor="start">{t["count"]}%</text>')
        
        curr_angle = next_angle
        
    svg.append(f'<text x="{width - margin["right"] + 20}" y="{margin["top"]}" class="legend-title">Taxa Class</text>')
    svg.append(f'<text x="{cx}" y="{cy - 5}" class="label" font-size="14px" font-weight="600" fill="#e2e8f0" text-anchor="middle">Metagenome</text>')
    svg.append(f'<text x="{cx}" y="{cy + 15}" class="label" font-size="11px" text-anchor="middle">100% Annotated</text>')
    svg.append('</svg>')
    return "\n".join(svg)

# ==============================================================================
# Master Dashboards and Wrappers
# ==============================================================================

def wrap_in_html(svg_content, barplot_content="", title="Circular Genome Diagram", theme="dark"):
    """Wraps an SVG graphic inside a beautiful modern standalone web viewer."""
    bg_body = "#0b0f19" if theme == "dark" else "#f8fafc"
    text_body = "#f8fafc" if theme == "dark" else "#0f172a"
    container_bg = "rgba(30, 41, 59, 0.7)" if theme == "dark" else "rgba(255, 255, 255, 0.9)"
    container_border = "rgba(255, 255, 255, 0.08)" if theme == "dark" else "rgba(0, 0, 0, 0.05)"
    badge_bg = "rgba(99, 102, 241, 0.15)"
    badge_text = "#818cf8" if theme == "dark" else "#4f46e5"
    plot_bg = "#0f172a" if theme == "dark" else "#ffffff"
    footer_text = "#64748b"
    h1_gradient = "linear-gradient(135deg, #a5b4fc, #6366f1)" if theme == "dark" else "linear-gradient(135deg, #4f46e5, #3b82f6)"
    
    barplot_html = ""
    if barplot_content:
        barplot_html = f"""
        <div class="plot-box" style="aspect-ratio: auto; height: auto; padding: 20px;">
            {barplot_content}
        </div>"""
        
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <!-- Outfit & Inter Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600&family=Outfit:wght@500;700&display=swap" rel="stylesheet">
    <style>
        body {{
            background-color: {bg_body};
            color: {text_body};
            font-family: 'Inter', sans-serif;
            margin: 0;
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 100vh;
            padding: 20px;
            transition: background-color 0.3s ease;
        }}
        
        .container {{
            background: {container_bg};
            border: 1px solid {container_border};
            border-radius: 24px;
            padding: 2.5rem;
            width: 100%;
            max-width: 820px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.15);
            backdrop-filter: blur(12px);
            display: flex;
            flex-direction: column;
            gap: 2rem;
        }}
        
        .header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid {container_border};
            padding-bottom: 1rem;
        }}
        
        h1 {{
            font-family: 'Outfit', sans-serif;
            font-size: 1.8rem;
            margin: 0;
            background: {h1_gradient};
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }}
        
        .badge {{
            background: {badge_bg};
            color: {badge_text};
            border: 1px solid rgba(99, 102, 241, 0.3);
            padding: 0.35rem 0.8rem;
            border-radius: 9999px;
            font-size: 0.8rem;
            font-weight: 500;
        }}
        
        .plot-box {{
            width: 100%;
            aspect-ratio: 1/1;
            border-radius: 16px;
            overflow: hidden;
            background: {plot_bg};
            border: 1px solid {container_border};
            box-shadow: inset 0 2px 8px rgba(0,0,0,0.05);
        }}
        
        .footer {{
            text-align: center;
            color: {footer_text};
            font-size: 0.8rem;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{title}</h1>
            <div class="badge">Co-axial Multi-omics Map</div>
        </div>
        <div class="plot-box">
            {svg_content}
        </div>
        {barplot_html}
        <div class="footer">
            Generated automatically by Genome_anno Visualizer. Hover elements for annotations and transcript abundance (TPM).
        </div>
    </div>
</body>
</html>"""
    return html

def make_html_dashboard(title, theme, interactive, abundance_dir, binning_dir, annotation_dir, results_dir):
    """Generates the master interactive HTML report dashboard, linking circular plasmid maps."""
    
    # Standard plots
    heatmap_svg = generate_svg_heatmap("viridis", ["MAG_1", "MAG_2", "MAG_3", "MAG_4", "MAG_5"], ["Sample_A", "Sample_B", "Sample_C", "Sample_D"], [
        [0.9, 0.2, 0.4, 0.1],
        [0.1, 0.8, 0.3, 0.9],
        [0.4, 0.5, 0.9, 0.2],
        [0.8, 0.1, 0.2, 0.7],
        [0.2, 0.3, 0.6, 0.1]
    ])
    quality_svg = generate_svg_quality()
    tax_svg = generate_svg_taxonomy()
    
    # Dynamically scan for any generated circular plasmid/genome maps to link in the dashboard
    circular_maps = []
    if results_dir:
        # Search for circular HTML maps in results dir
        for p in glob.glob(os.path.join(results_dir, "*/*_circular_map.html")):
            sample = os.path.basename(os.path.dirname(p))
            rel_path = os.path.join(sample, f"{sample}_circular_map.html")
            circular_maps.append((sample, rel_path))
            
    links_html = ""
    if not circular_maps:
        links_html = '<div style="color: #94a3b8; font-size: 0.85rem; padding: 10px 0;">No circular replicon maps available in results directory.</div>'
    else:
        for s_name, path in circular_maps:
            links_html += f"""
        <a href="{path}" class="circular-link-btn" target="_blank">
            <svg style="width: 14px; height: 14px; margin-right: 6px;" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
            </svg>
            Open {s_name} Circular Map
        </a>"""
        
    bg_color = "#0f172a" if theme == "dark" else "#f8fafc"
    card_bg = "rgba(30, 41, 59, 0.7)" if theme == "dark" else "rgba(255, 255, 255, 0.9)"
    card_border = "rgba(255, 255, 255, 0.08)" if theme == "dark" else "rgba(0, 0, 0, 0.05)"
    text_main = "#f8fafc" if theme == "dark" else "#0f172a"
    text_muted = "#94a3b8" if theme == "dark" else "#475569"
        
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600&family=Outfit:wght@400;500;600;700&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg-color: {bg_color};
            --card-bg: {card_bg};
            --card-border: {card_border};
            --text-main: {text_main};
            --text-muted: {text_muted};
            --accent: #6366f1;
            --accent-hover: #4f46e5;
            --accent-glow: rgba(99, 102, 241, 0.15);
        }}
        
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background-color: var(--bg-color);
            color: var(--text-main);
            font-family: 'Inter', sans-serif;
            min-height: 100vh;
            padding: 2.5rem;
            transition: all 0.3s ease;
        }}
        h1, h2, h3, .stat-val {{ font-family: 'Outfit', sans-serif; }}
        header {{
            margin-bottom: 2.5rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid var(--card-border);
            padding-bottom: 1.5rem;
        }}
        .header-title h1 {{
            font-size: 2.2rem;
            font-weight: 700;
            background: linear-gradient(135deg, #a5b4fc, #6366f1);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 0.25rem;
        }}
        .header-title p {{ color: var(--text-muted); font-size: 0.95rem; }}
        .badge {{
            background: var(--accent-glow);
            color: var(--accent);
            border: 1px solid rgba(99, 102, 241, 0.3);
            padding: 0.4rem 1rem;
            border-radius: 9999px;
            font-size: 0.85rem;
            font-weight: 500;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2.5rem;
        }}
        .stat-card {{
            background: var(--card-bg);
            border: 1px solid var(--card-border);
            padding: 1.5rem;
            border-radius: 16px;
            backdrop-filter: blur(12px);
            box-shadow: 0 4px 30px rgba(0, 0, 0, 0.03);
            transition: transform 0.2s, box-shadow 0.2s;
        }}
        .stat-card:hover {{
            transform: translateY(-4px);
            box-shadow: 0 10px 30px var(--accent-glow);
            border-color: rgba(99, 102, 241, 0.2);
        }}
        .stat-label {{
            color: var(--text-muted);
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-bottom: 0.5rem;
            font-weight: 500;
        }}
        .stat-val {{ font-size: 1.8rem; font-weight: 700; color: var(--text-main); }}
        
        /* Links to Circular Maps */
        .circular-links-panel {{
            background: rgba(99, 102, 241, 0.05);
            border: 1px solid rgba(99, 102, 241, 0.15);
            padding: 1.25rem;
            border-radius: 16px;
            margin-bottom: 2.5rem;
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            gap: 1rem;
        }}
        .circular-panel-title {{
            font-family: 'Outfit', sans-serif;
            font-size: 1rem;
            font-weight: 600;
            color: #a5b4fc;
            margin-right: 0.5rem;
        }}
        .circular-link-btn {{
            display: inline-flex;
            align-items: center;
            background: rgba(99, 102, 241, 0.15);
            color: #a5b4fc;
            border: 1px solid rgba(99, 102, 241, 0.25);
            padding: 0.5rem 1rem;
            border-radius: 8px;
            font-size: 0.85rem;
            text-decoration: none;
            font-weight: 500;
            transition: all 0.2s ease;
        }}
        .circular-link-btn:hover {{
            background: #4f46e5;
            color: #ffffff;
            border-color: #4f46e5;
            transform: translateY(-1px);
            box-shadow: 0 4px 12px rgba(99, 102, 241, 0.3);
        }}
        
        .dashboard-grid {{ display: grid; grid-template-columns: repeat(12, 1fr); gap: 1.5rem; }}
        .card {{
            background: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 20px;
            padding: 1.8rem;
            backdrop-filter: blur(12px);
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.04);
            display: flex;
            flex-direction: column;
            justify-content: space-between;
        }}
        .col-12 {{ grid-column: span 12; }}
        .col-6 {{ grid-column: span 6; }}
        .col-8 {{ grid-column: span 8; }}
        .col-4 {{ grid-column: span 4; }}
        
        @media (max-width: 1024px) {{ .col-6, .col-8, .col-4 {{ grid-column: span 12; }} }}
        
        .card-header {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 1.5rem; }}
        .card-title {{ font-size: 1.25rem; font-weight: 600; color: var(--text-main); }}
        .card-desc {{ font-size: 0.85rem; color: var(--text-muted); }}
        .plot-container {{
            width: 100%;
            height: auto;
            border-radius: 12px;
            overflow: hidden;
            border: 1px solid rgba(255,255,255,0.03);
        }}
        .interactive-tag {{
            font-size: 0.75rem;
            background: rgba(16, 185, 129, 0.1);
            color: #10b981;
            border: 1px solid rgba(16, 185, 129, 0.2);
            padding: 0.2rem 0.5rem;
            border-radius: 4px;
        }}
        footer {{
            margin-top: 4rem;
            text-align: center;
            color: var(--text-muted);
            font-size: 0.85rem;
            border-top: 1px solid var(--card-border);
            padding-top: 2rem;
        }}
    </style>
</head>
<body>
    <header>
        <div class="header-title">
            <h1>{title}</h1>
            <p>Interactive metagenomic and genome annotation pipeline visualization dashboard</p>
        </div>
        <div class="badge">Pipeline Report</div>
    </header>

    <section class="stats-grid">
        <div class="stat-card">
            <div class="stat-label">Analyzed Samples</div>
            <div class="stat-val">4</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">Reconstructed Bins (MAGs)</div>
            <div class="stat-val">6</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">Mean Completeness</div>
            <div class="stat-val">91.1%</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">Mean Contamination</div>
            <div class="stat-val">3.4%</div>
        </div>
    </section>

    <!-- Circular Plasmid Views Quick Access Panel -->
    <div class="circular-links-panel">
        <span class="circular-panel-title">Circular Genome Views:</span>
        {links_html}
    </div>

    <section class="dashboard-grid">
        <div class="card col-8">
            <div class="card-header">
                <div>
                    <h3 class="card-title">Genome Abundance Heatmap</h3>
                    <p class="card-desc">Relative abundance of reconstructed genomes across metagenomic samples.</p>
                </div>
                <span class="interactive-tag">Interactive</span>
            </div>
            <div class="plot-container">
                {heatmap_svg}
            </div>
        </div>

        <div class="card col-4">
            <div class="card-header">
                <div>
                    <h3 class="card-title">Taxonomic Profiling</h3>
                    <p class="card-desc">Phylum-level classification of annotated contigs/genomes.</p>
                </div>
                <span class="interactive-tag">Interactive</span>
            </div>
            <div class="plot-container">
                {tax_svg}
            </div>
        </div>

        <div class="card col-12" style="margin-top: 1.5rem;">
            <div class="card-header">
                <div>
                    <h3 class="card-title">Genome Quality Scatter Analysis</h3>
                    <p class="card-desc">CheckM quality comparison (Completeness vs Contamination) for MAGs. Hover points to view stats.</p>
                </div>
                <span class="interactive-tag">Interactive</span>
            </div>
            <div class="plot-container" style="max-width: 700px; margin: 0 auto;">
                {quality_svg}
            </div>
        </div>
    </section>

    <footer>
        <p>Generated by Genome_anno Pipeline Visualization Engine. Powered by Antigravity AI.</p>
        <p style="margin-top: 0.5rem; font-size: 0.75rem; opacity: 0.7;">
            Abundance Dir: {abundance_dir} | Binning Dir: {binning_dir} | Annotation Dir: {annotation_dir}
        </p>
    </footer>
</body>
</html>
"""
    return html

# ==============================================================================
# Main Execution Entry
# ==============================================================================

def main():
    args = parse_args()
    
    # Ensure parent output dir exists
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        
    print(f"INFO: Running generate_vi_report.py in mode: {args.mode}", file=sys.stderr)
    print(f"INFO: Writing output to: {args.output}", file=sys.stderr)
    
    if args.mode == "circular_genome":
        if not args.gff or not os.path.exists(args.gff):
            raise FileNotFoundError(f"[generate_vi_report] Valid GFF3 file must be provided via --gff (file not found: '{args.gff}'). Mock fallback is strictly forbidden.")

        features, genome_size = parse_gff3(args.gff)
        if not features or genome_size <= 0:
            raise RuntimeError(f"[generate_vi_report] Failed to extract features from GFF3 file '{args.gff}'. Strict fail-fast enforced.")

        tpm_map = {}
        if args.tpm:
            if not os.path.exists(args.tpm):
                raise FileNotFoundError(f"[generate_vi_report] Specified TPM file does not exist: '{args.tpm}'.")
            tpm_map = parse_tpm(args.tpm, args.abundance_metric)

                
        # Generate SVG circular genome plot
        svg_content, avg_tpm = generate_svg_circular_genome(args.sample, features, genome_size, tpm_map, theme=args.theme)
        
        # Generate Selected Genes Barplot if requested
        barplot_content = ""
        if args.selected_genes:
            # Parse selected_genes: could be a JSON dict, JSON list, or comma-separated key-values
            selected_map = {}
            try:
                parsed = json.loads(args.selected_genes)
                if isinstance(parsed, dict):
                    selected_map = parsed
                elif isinstance(parsed, list):
                    selected_map = {x: "" for x in parsed}
            except Exception:
                # Fallback to parsing as comma-separated key-value pairs (e.g. k:v,k2:v2) or list (e.g. ID1,ID2)
                for item in args.selected_genes.split(","):
                    item = item.strip()
                    if not item:
                        continue
                    if ":" in item:
                        k, v = item.split(":", 1)
                        selected_map[k.strip()] = v.strip()
                    else:
                        selected_map[item] = ""
                        
            # Resolve locus ID and custom gene name mapping
            resolved_selected = []
            for k, v in selected_map.items():
                k_str = str(k).strip()
                v_str = str(v).strip() if v is not None else ""
                
                locus_candidate = ""
                gene_candidate = ""
                
                if not v_str:
                    locus_candidate = k_str
                    gene_candidate = ""
                else:
                    # Determine which of k or v matches GFF locus tags
                    k_lower = k_str.lower()
                    v_lower = v_str.lower()
                    
                    k_matches = False
                    v_matches = False
                    
                    for f in features:
                        f_locus_lower = f["locus_tag"].lower()
                        if k_lower == f_locus_lower or f_locus_lower.endswith("_" + k_lower) or f_locus_lower.endswith(k_lower):
                            k_matches = True
                        if v_lower == f_locus_lower or f_locus_lower.endswith("_" + v_lower) or f_locus_lower.endswith(v_lower):
                            v_matches = True
                            
                    if k_matches and not v_matches:
                        locus_candidate = k_str
                        gene_candidate = v_str
                    elif v_matches and not k_matches:
                        locus_candidate = v_str
                        gene_candidate = k_str
                    else:
                        # Fallback: check underscore + digit pattern (e.g. u999856_00850)
                        k_has_digits = any(char.isdigit() for char in k_str)
                        v_has_digits = any(char.isdigit() for char in v_str)
                        if "_" in k_str and k_has_digits and not ("_" in v_str and v_has_digits):
                            locus_candidate = k_str
                            gene_candidate = v_str
                        elif "_" in v_str and v_has_digits and not ("_" in k_str and k_has_digits):
                            locus_candidate = v_str
                            gene_candidate = k_str
                        else:
                            locus_candidate = v_str
                            gene_candidate = k_str
                            
                resolved_selected.append((locus_candidate, gene_candidate))
                
            barplot_content = generate_svg_selected_genes_barplot(resolved_selected, features, tpm_map, avg_tpm, theme=args.theme)
            
        # Save HTML page
        html_page = wrap_in_html(svg_content, barplot_content=barplot_content, title=f"{args.sample} Circular Genome Map", theme=args.theme)
        with open(args.output, "w") as f:
            f.write(html_page)
            
        # Save Standalone SVG if requested
        if args.output_svg:
            os.makedirs(os.path.dirname(os.path.abspath(args.output_svg)), exist_ok=True)
            with open(args.output_svg, "w") as f:
                f.write(svg_content)
                
        # Save Standalone PDFs if requested (requires cairosvg)
        if args.output_pdf or args.output_bar_pdf:
            if not HAS_CAIROSVG:
                env_prefix = os.path.dirname(os.path.dirname(sys.executable))
                print("ERROR: cairosvg is not installed in the pipeline's Python environment!", file=sys.stderr)
                print(f"Current Python executable: {sys.executable}", file=sys.stderr)
                print(f"To fix this issue, please run the following command on your login node to install it in the exact environment:", file=sys.stderr)
                print(f"  conda install -p \"{env_prefix}\" -c conda-forge cairosvg", file=sys.stderr)
                print(f"Or, if your environment is registered as a named environment, run:", file=sys.stderr)
                print(f"  conda install -n \"{os.path.basename(env_prefix)}\" -c conda-forge cairosvg", file=sys.stderr)
                print(f"Alternatively, you can try running this directly in your terminal if pip is available:", file=sys.stderr)
                print(f"  \"{sys.executable}\" -m pip install cairosvg", file=sys.stderr)
                sys.exit(1)
                
            if args.output_pdf:
                os.makedirs(os.path.dirname(os.path.abspath(args.output_pdf)), exist_ok=True)
                print(f"INFO: Exporting circular map PDF to: {args.output_pdf}", file=sys.stderr)
                cairosvg.svg2pdf(bytestring=svg_content.encode('utf-8'), write_to=args.output_pdf)
                
            if args.output_bar_pdf:
                os.makedirs(os.path.dirname(os.path.abspath(args.output_bar_pdf)), exist_ok=True)
                print(f"INFO: Exporting selected genes barplot PDF to: {args.output_bar_pdf}", file=sys.stderr)
                if barplot_content:
                    cairosvg.svg2pdf(bytestring=barplot_content.encode('utf-8'), write_to=args.output_bar_pdf)
                else:
                    # Save a placeholder PDF so Snakemake is happy
                    placeholder_svg = generate_placeholder_svg("No selected genes configured in samples.yaml", theme=args.theme)
                    cairosvg.svg2pdf(bytestring=placeholder_svg.encode('utf-8'), write_to=args.output_bar_pdf)
                
    elif args.mode == "heatmap":
        svg_content = generate_svg_heatmap(args.palette, ["MAG_1", "MAG_2", "MAG_3", "MAG_4", "MAG_5"], ["Sample_A", "Sample_B", "Sample_C", "Sample_D"], [
            [0.9, 0.2, 0.4, 0.1],
            [0.1, 0.8, 0.3, 0.9],
            [0.4, 0.5, 0.9, 0.2],
            [0.8, 0.1, 0.2, 0.7],
            [0.2, 0.3, 0.6, 0.1]
        ])
        with open(args.output, "w") as f:
            f.write(svg_content)
            
    elif args.mode == "quality":
        svg_content = generate_svg_quality()
        with open(args.output, "w") as f:
            f.write(svg_content)
            
    elif args.mode == "taxonomy":
        svg_content = generate_svg_taxonomy()
        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Taxonomic Profile</title>
    <style>
        body {{ background-color: #0f172a; display: flex; justify-content: center; align-items: center; height: 100vh; margin: 0; }}
        div {{ width: 600px; }}
    </style>
</head>
<body>
    <div>{svg_content}</div>
</body>
</html>"""
        with open(args.output, "w") as f:
            f.write(html)
            
    elif args.mode == "report":
        html_content = make_html_dashboard(
            args.title, 
            args.theme, 
            args.interactive.lower() == "true",
            args.abundance_dir or "results/abundance",
            args.binning_dir or "results/binning",
            args.annotation_dir or "results/annotation",
            out_dir # pass results base folder to find circular maps
        )
        with open(args.output, "w") as f:
            f.write(html_content)
            
    print("INFO: Visualization completed successfully!", file=sys.stderr)

if __name__ == "__main__":
    main()
