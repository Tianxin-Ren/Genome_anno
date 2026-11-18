import sys
import os

def load_ko_list(ko_list_path):
    ko_data = {}
    print(f"Loading metadata from {ko_list_path}...")
    with open(ko_list_path, 'r') as f:
        # Skip header knum threshold ...
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            # Based on user sample: knum=0, threshold=1, definition=11
            if len(parts) >= 2:
                ko_id = parts[0]
                try:
                    threshold = float(parts[1])
                    # Corrected index to 11 for definition
                    definition = parts[11] if len(parts) > 11 else "hypothetical protein"
                    ko_data[ko_id] = {"threshold": threshold, "desc": definition}
                except ValueError:
                    continue
    return ko_data

def patch_hmm_with_bakta_format(input_hmm, output_hmm, ko_data, default_val=25.0):
    print(f"Patching HMM for Bakta compatibility: {input_hmm} -> {output_hmm}")
    current_ko = None
    with open(input_hmm, 'r') as fin, open(output_hmm, 'w') as fout:
        for line in fin:
            if line.startswith("NAME  "):
                current_ko = line.strip().split()[1]
            
            # Skip existing DESC to replace it with Bakta-compliant one
            if line.startswith("DESC  "):
                continue

            # Inject before HMM data starts
            if line.startswith("HMM "):
                data = ko_data.get(current_ko, {"threshold": default_val, "desc": "hypothetical protein"})
                
                # 1. Inject Bakta-compliant DESC: gene~~~product~~~dbxrefs~~~type
                # We use: KO_ID~~~Description~~~KO:KO_ID~~~CDS
                fout.write(f"DESC  {current_ko}~~~{data['desc']}~~~KO:{current_ko}~~~CDS\n")
                
                # 2. Inject thresholds
                thresh_str = f"{data['threshold']:.2f} {data['threshold']:.2f}"
                fout.write(f"GA    {thresh_str}\n")
                fout.write(f"TC    {thresh_str}\n")
                fout.write(f"NC    {thresh_str}\n")
            
            fout.write(line)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python fix_hmm.py <input.hmm> <ko_list> <output.hmm>")
    else:
        input_hmm = sys.argv[1]
        ko_list = sys.argv[2]
        output_hmm = sys.argv[3]
        
        data = load_ko_list(ko_list)
        patch_hmm_with_bakta_format(input_hmm, output_hmm, data)
