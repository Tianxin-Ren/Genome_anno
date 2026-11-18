import sys
import argparse
from Bio import SeqIO

def extract_region(input_gbk, output_gbk, region_str):
    """Extracts a specific coordinate range from a GenBank file."""
    try:
        # Support both 'start-end' and 'start..end'
        if "-" in region_str:
            start, end = map(int, region_str.split('-'))
        elif ".." in region_str:
            start, end = map(int, region_str.split('..'))
        else:
            raise ValueError
    except ValueError:
        print(f"Error: Invalid region format '{region_str}'. Expected 'start-end' (e.g. 1000-5000).")
        sys.exit(1)

    # 1-indexed to 0-indexed conversion
    start_idx = max(0, start - 1)
    end_idx = end

    records = list(SeqIO.parse(input_gbk, "genbank"))
    if not records:
        print(f"Error: No records found in {input_gbk}")
        sys.exit(1)

    # cblaster queries are typically single records. 
    # If multiple records exist, we apply the crop to the first one that fits or the first one total.
    record = records[0]
    
    if start_idx >= len(record):
        print(f"Error: Start coordinate {start} is beyond record length {len(record)}.")
        sys.exit(1)
        
    sub_record = record[start_idx:end_idx]
    
    # Update metadata to reflect the slice
    sub_record.id = f"{record.id}_{start}_{end}"
    sub_record.name = sub_record.id[:16] # GenBank name limit is 16 chars
    sub_record.description = f"{record.description} | Region: {start}-{end}"

    # Write output
    with open(output_gbk, "w") as out_f:
        SeqIO.write(sub_record, out_f, "genbank")
    
    print(f"Successfully extracted region {start}-{end} from {record.id} to {output_gbk}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract a region from a GenBank file.")
    parser.add_argument("input", help="Input GBK file")
    parser.add_argument("output", help="Output GBK file")
    parser.add_argument("region", help="Region to extract (start-end, 1-indexed)")
    args = parser.parse_args()
    
    extract_region(args.input, args.output, args.region)
