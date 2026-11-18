import os
import sys
import pandas as pd
from Bio import SeqIO
import argparse

def extract_potential(fasta_in, genomad_tsv, fasta_out, selection_log, min_score):
    gn = pd.read_csv(genomad_tsv, sep='\t')
    gn['norm_id'] = gn['seq_name'].str.split().str[0]
    ai_set = set(gn[gn['plasmid_score'] >= min_score]['norm_id'])
    
    selected_records = []
    source_data = [] 
    
    for record in SeqIO.parse(fasta_in, "fasta"):
        norm_id = record.id.split()[0]
        length = len(record.seq)
        header = record.description
        is_circular = "circular-yes" in header or "circular-possibly" in header
        in_range = 2000 <= length <= 1000000 
        is_ai_hit = norm_id in ai_set
        
        if (is_circular and in_range) or is_ai_hit:
            selected_records.append(record)
            source = "AI+Circular" if (is_circular and is_ai_hit) else ("Circular_Only" if is_circular else "AI_Only")
            source_data.append({"norm_id": norm_id, "Evidence_Source": source})

    with open(fasta_out, "w") as f:
        if not selected_records: 
            f.write(">dummy_placeholder\nNNNNN\n")
        else: 
            SeqIO.write(selected_records, f, "fasta")
            
    pd.DataFrame(source_data).to_csv(selection_log, sep='\t', index=False)

def extract_fastas(fasta_in, report_tsv, out_dir, done_file):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    report = pd.read_csv(report_tsv, sep='\t')
    if report.empty or (len(report)==1 and report.iloc[0]['norm_id']=="-"):
        with open(done_file, 'w') as f: f.write("No plasmids found.")
        return
    # Robust ID matching: simplify IDs from the report to match our target format
    plasmid_ids = set([str(x).split('ctg')[0] for x in report['norm_id'].tolist()])
    count = 0
    for record in SeqIO.parse(fasta_in, "fasta"):
        rec_id = record.id.split()[0]
        # Simplify rec_id for comparison with the report
        rec_id_simple = rec_id.split('ctg')[0]
        
        if rec_id_simple in plasmid_ids:
            # Simplify ID for saving: u987790ctg_len... -> u987790
            simple_id = rec_id_simple
            record.id = simple_id
            record.description = simple_id
            
            out_path = os.path.join(out_dir, f"{simple_id}.fa")
            with open(out_path, "w") as out_f:
                SeqIO.write(record, out_f, "fasta")
            count += 1
    with open(done_file, 'w') as f: f.write(f"Extraction complete. {count} plasmids saved.")

def merge_reports(genomad_tsv, selection_tsv, typer_report, out_tsv):
    import re
    if not os.path.exists(selection_tsv): return
    sel = pd.read_csv(selection_tsv, sep='\t')
    gn = pd.read_csv(genomad_tsv, sep='\t')
    gn['norm_id'] = gn['seq_name'].str.split().str[0]
    gn_cols = ['norm_id', 'plasmid_score', 'n_genes', 'n_hallmarks', 'marker_enrichment', 'conjugation_genes', 'amr_genes']
    gn_data = gn[gn_cols].copy()
    if os.path.exists(typer_report) and os.path.getsize(typer_report) > 0:
        mr = pd.read_csv(typer_report, sep='\t')
        mr_cols = {'norm_id': 'norm_id', 'rep_type(s)': 'rep_type(s)', 'relaxase_type(s)': 'relaxase_type(s)', 'predicted_mobility': 'predicted_mobility','primary_cluster_id': 'primary_cluster_id', 'mash_neighbor_distance': 'mash_dist'}
        valid_cols = [c for c in mr_cols.keys() if c in mr.columns]
        mr_data = mr[valid_cols].copy().rename(columns={k: v for k, v in mr_cols.items() if k in valid_cols})
    else: mr_data = pd.DataFrame(columns=['norm_id', 'primary_cluster_id', 'mash_dist'])
    merged = pd.merge(sel, gn_data, on='norm_id', how='left')
    merged = pd.merge(merged, mr_data, on='norm_id', how='left')
    def parse_header(name):
        if pd.isna(name) or name == "-": return "-", "-", "-"
        circ = re.search(r"circular-([^_ ]+)", name); depth = re.search(r"depth-([^_ ]+)", name); mult = re.search(r"mult=([^_ ]+)", name)
        return (circ.group(1) if circ else "-", depth.group(1) if depth else "-", mult.group(1) if mult else "-")
    meta = merged['norm_id'].apply(parse_header)
    merged['Myloasm_Circular'] = [x[0] for x in meta]; merged['Myloasm_Depth'] = [x[1] for x in meta]; merged['Myloasm_Multiplicity'] = [x[2] for x in meta]
    merged = merged.sort_values(by='plasmid_score', ascending=False) if 'plasmid_score' in merged.columns else merged
    merged = merged.fillna("-")
    final_cols = ['norm_id', 'Evidence_Source', 'plasmid_score', 'n_genes', 'n_hallmarks', 'marker_enrichment', 'primary_cluster_id', 'rep_type(s)', 'relaxase_type(s)', 'predicted_mobility','mash_dist', 'Myloasm_Circular', 'Myloasm_Depth', 'Myloasm_Multiplicity']
    for col in final_cols:
        if col not in merged.columns: merged[col] = "-"
    merged['norm_id'] = merged['norm_id'].apply(lambda x: str(x).split('ctg')[0])
    merged[final_cols].to_csv(out_tsv, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    p1 = subparsers.add_parser("potential"); p1.add_argument("--fasta", required=True); p1.add_argument("--tsv", required=True); p1.add_argument("--out", required=True); p1.add_argument("--log", required=True); p1.add_argument("--score", type=float, default=0.7)
    p2 = subparsers.add_parser("extract"); p2.add_argument("--fasta", required=True); p2.add_argument("--report", required=True); p2.add_argument("--outdir", required=True); p2.add_argument("--done", required=True)
    p3 = subparsers.add_parser("merge"); p3.add_argument("--genomad", required=True); p3.add_argument("--selection", required=True); p3.add_argument("--typer", required=True); p3.add_argument("--out", required=True)
    args = parser.parse_args()
    if args.command == "potential": extract_potential(args.fasta, args.tsv, args.out, args.log, args.score)
    elif args.command == "extract": extract_fastas(args.fasta, args.report, args.outdir, args.done)
    elif args.command == "merge": merge_reports(args.genomad, args.selection, args.typer, args.out)
