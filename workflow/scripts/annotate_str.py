#!/usr/bin/env python

import pysam
import logging

log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None
logging.basicConfig(
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
)


def get_indel_length(record):
    """Calculates the length of an indel from Manta."""
    # For symbolic variants (e.g., <DEL>), Manta often uses INFO/SVLEN
    if 'SVLEN' in record.info:
        val = record.info['SVLEN']
        return abs(val[0]) if isinstance(val, (list, tuple)) else abs(val)
    # For explicit sequences (standard indels)
    return max(len(record.ref), len(record.alts[0]))


def calculate_overlap(chrom, start, end, bed_tabix):
    """Finds the total number of overlapping bases in the BED file."""
    overlap_bp = 0
    try:
        for region in bed_tabix.fetch(chrom, start, end):
            fields = region.split('\t')
            b_start = int(fields[1])
            b_end = int(fields[2])
            o_start = max(start, b_start)
            o_end = min(end, b_end)
            if o_start < o_end:
                overlap_bp += (o_end - o_start)
    except ValueError:
        # Chromosome is not present in the BED file
        pass
    return overlap_bp


def main():
    logging.info("Starting STR annotation...")

    # 1. Fetch paths directly from the Snakemake object
    vcf_in_path = snakemake.input.vcf
    vcf_out_path = snakemake.output.vcf
    bed_gz_path = snakemake.input.bed

    logging.info(f"Input VCF: {vcf_in_path}")
    logging.info(f"Input BED (Simple Repeats): {bed_gz_path}")
    logging.info(f"Output VCF: {vcf_out_path}")

    # Open files using pysam
    vcf_in = pysam.VariantFile(vcf_in_path)
    bed = pysam.TabixFile(bed_gz_path)

    # Add a new header for the INFO field
    vcf_in.header.info.add('STR_PERCENT', '1', 'Float', 'Percentage of indel overlapping with STR (1-100)')

    records_processed = 0
    records_annotated = 0

    with pysam.VariantFile(vcf_out_path, 'w', header=vcf_in.header) as vcf_out:
        for record in vcf_in:
            records_processed += 1

            # Manta often starts at POS and ends at POS + length.
            # We convert to 0-based start/end for the interval search.
            v_start = record.start
            v_len = get_indel_length(record)
            v_end = v_start + v_len

            if v_len > 0:
                overlap_bp = calculate_overlap(record.chrom, v_start, v_end, bed)
                percent = round((overlap_bp / v_len) * 100, 2)
                if percent > 1.0:
                    record.info['STR_PERCENT'] = min(percent, 100.0)
                    records_annotated += 1

            vcf_out.write(record)

    logging.info(f"Finished! Processed {records_processed} variants. Annotated {records_annotated} with STR_PERCENT.")


if __name__ == "__main__":
    main()
