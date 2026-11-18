#!/usr/bin/env python3
import os, sys, argparse

def get_faa_map(faa):
    """Build a map of sequence -> ID for the original assembly proteins."""
    seq_map = {}
    with open(faa) as f:
        cid = ""
        seq = []
        for line in f:
            if line.startswith(">"):
                if cid: seq_map["".join(seq)] = cid
                cid = line.split()[0][1:]
                seq = []
            else:
                seq.append(line.strip())
        if cid: seq_map["".join(seq)] = cid
    return seq_map

def map_eggnog(in_file, out_file, mapping):
    if in_file == "none" or not os.path.exists(in_file):
        with open(out_file, 'w') as f: f.write("# No EggNOG input provided\n")
        return 0
    
    count = 0
    with open(in_file) as fin, open(out_file, 'w') as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.split("\t")
            old_id = parts[0]
            # Since mapping is refined_id -> assembly_id, we need to reverse it or search it
            # But the script's logic is usually to find which refined gene matches which assembly gene
            # Let's assume the mapping passed in is refined_id -> assembly_id
            for refined_id, assembly_id in mapping.items():
                if assembly_id == old_id:
                    new_parts = [refined_id] + parts[1:]
                    fout.write("\t".join(new_parts))
                    count += 1
    return count

def map_kofam(in_file, out_file, mapping):
    if in_file == "none" or not os.path.exists(in_file):
        with open(out_file, 'w') as f: f.write("# No KOfam input provided\n")
        return 0
    
    count = 0
    with open(in_file) as fin, open(out_file, 'w') as fout:
        for line in fin:
            if line.startswith("#") or line.startswith("-"):
                fout.write(line)
                continue
            # Handle Kofam output which might have leading '*'
            clean_line = line.lstrip("*").strip()
            parts = clean_line.split()
            if not parts: continue
            old_id = parts[0]
            for refined_id, assembly_id in mapping.items():
                if assembly_id == old_id:
                    # Keep original formatting as much as possible
                    fout.write(line.replace(old_id, refined_id, 1))
                    count += 1
    return count

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refined_faa", required=True)
    parser.add_argument("--assembly_faa", required=True)
    parser.add_argument("--assembly_eggnog", default="none")
    parser.add_argument("--assembly_kofam", default="none")
    parser.add_argument("--out_eggnog", required=True)
    parser.add_argument("--out_kofam", required=True)
    args = parser.parse_args()

    print(f"Building hash map from {args.assembly_faa}...")
    assembly_map = get_faa_map(args.assembly_faa)

    print(f"Mapping refined proteins from {args.refined_faa}...")
    final_map = {} # refined_id -> assembly_id
    with open(args.refined_faa) as f:
        cid = ""
        seq = []
        for line in f:
            if line.startswith(">"):
                if cid:
                    s = "".join(seq)
                    if s in assembly_map:
                        final_map[cid] = assembly_map[s]
                cid = line.split()[0][1:]
                seq = []
            else:
                seq.append(line.strip())
        if cid:
            s = "".join(seq)
            if s in assembly_map:
                final_map[cid] = assembly_map[s]

    print(f"Processing EggNOG mapping...")
    e_count = map_eggnog(args.assembly_eggnog, args.out_eggnog, final_map)
    print(f"Mapped {e_count} EggNOG entries.")

    print(f"Processing KofamScan mapping...")
    k_count = map_kofam(args.assembly_kofam, args.out_kofam, final_map)
    print(f"Mapped {k_count} KOfam entries.")

if __name__ == "__main__":
    main()
