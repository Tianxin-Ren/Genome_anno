import os
import pandas as pd
import subprocess
from Bio import SeqIO
import argparse

def run_batch(fasta_in, out_tsv, work_dir, env_path):
    os.makedirs(work_dir, exist_ok=True)
    records = list(SeqIO.parse(fasta_in, "fasta"))
    results = []
    for record in records:
        norm_id = record.id.split()[0]
        if norm_id == "dummy_placeholder": continue
        contig_fa = os.path.join(work_dir, f"{norm_id}.fa")
        contig_out = os.path.join(work_dir, f"{norm_id}_typer.txt")
        with open(contig_fa, "w") as out: SeqIO.write(record, out, "fasta")
        cmd = f"conda run -p {env_path} mob_typer --infile {contig_fa} --out_file {contig_out} --num_threads 1 --plasmid_meta --plasmid_reports"
        subprocess.run(cmd, shell=True, check=False)
        if os.path.exists(contig_out) and os.path.getsize(contig_out) > 0:
            try:
                df = pd.read_csv(contig_out, sep='\t')
                if not df.empty:
                    res = df.iloc[0].to_dict()
                    res['norm_id'] = norm_id
                    results.append(res)
            except: pass
    if results:
        pd.DataFrame(results).to_csv(out_tsv, sep='\t', index=False)
    else:
        pd.DataFrame(columns=['norm_id', 'primary_cluster_id']).to_csv(out_tsv, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True); parser.add_argument("--out", required=True); parser.add_argument("--tmp", required=True); parser.add_argument("--env", required=True)
    args = parser.parse_args()
    run_batch(args.fasta, args.out, args.tmp, args.env)
