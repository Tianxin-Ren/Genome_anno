import os
import sys
from Bio import SeqIO

def fix_ids(gff_in, faa_in, tsv_in, fasta_in, gff_out, faa_out, tsv_out):
    # 1. Get original IDs and simplify them
    original_ids = []
    for record in SeqIO.parse(fasta_in, "fasta"):
        simple_id = record.id.split('ctg')[0]
        original_ids.append(simple_id)
    
    print(f"Simplifying {len(original_ids)} contig IDs for chunk {fasta_in}")

    # 2. Map contig_1, contig_2... to simplified IDs
    mapping = {f"contig_{i+1}": name for i, name in enumerate(original_ids)}

    # 3. Fix GFF3
    if os.path.exists(gff_in):
        with open(gff_in, 'r') as f_in, open(gff_out, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    # Fix ##sequence-region lines
                    if line.startswith('##sequence-region'):
                        parts = line.split()
                        if len(parts) > 1:
                            # If it's contig_i, map it. Otherwise, simplify it.
                            if parts[1] in mapping:
                                parts[1] = mapping[parts[1]]
                            else:
                                parts[1] = parts[1].split('ctg')[0]
                            line = " ".join(parts) + "\n"
                    f_out.write(line)
                    continue
                parts = line.split('\t')
                # Apply simplification to column 1
                if parts[0] in mapping:
                    parts[0] = mapping[parts[0]]
                else:
                    parts[0] = parts[0].split('ctg')[0]
                f_out.write("\t".join(parts))

    # 4. Fix FAA (headers)
    if os.path.exists(faa_in):
        import re
        with open(faa_in, 'r') as f_in, open(faa_out, 'w') as f_out:
            for line in f_in:
                if line.startswith('>'):
                    # Simplify [contig=u124ctg_len...] -> [contig=u124]
                    if '[contig=' in line:
                        line = re.sub(r'\[contig=[^\]]+ctg[^\]]*\]', 
                                     lambda m: m.group(0).split('ctg')[0] + ']', 
                                     line)
                f_out.write(line)

    # 5. Fix TSV
    if os.path.exists(tsv_in):
        with open(tsv_in, 'r') as f_in, open(tsv_out, 'w') as f_out:
            header = f_in.readline()
            f_out.write(header)
            for line in f_in:
                # TSV might have contig column in some versions, but Bakta usually doesn't.
                # If needed, we could add logic here.
                f_out.write(line)

if __name__ == "__main__":
    fix_ids(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
