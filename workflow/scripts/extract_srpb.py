import json
import os

def main():
    input_json = snakemake.input["json"]
    output_tsv = snakemake.output["tsv"]
    sample_id = snakemake.wildcards.sample
    
    with open(input_json, 'r') as f:
        data = json.load(f)

    srpb_map = {}
    for item in data.get("rearrangements", []):
        item_id = item.get("id")
        if item_id in ("01", "02"):
            srpb = item.get("evidence", {}).get("SRPB")
            srpb_map[item_id] = srpb

    core_val = srpb_map.get("01", "")
    extended_val = srpb_map.get("02", "")

    # Write to tab-separated output file
    # Columns: sample, core, extended
    with open(output_tsv, 'w') as f:
        f.write("sample\tSRPB_core\tSRPB_extended\n")
        f.write(f"{sample_id}\t{core_val}\t{extended_val}\n")
        f.write("# SRPB should be higher than 5 for the core region and 15 for the extended region to signify presence of a DUX4-IGH fusion\n")

if __name__ == "__main__":
    main()

