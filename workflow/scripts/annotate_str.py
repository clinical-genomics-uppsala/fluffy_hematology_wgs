import pysam
import sys
import os

def get_indel_length(record):
    """Beräknar längden på en indel från Manta."""
    # För symboliska varianter (t.ex. <DEL>) använder Manta ofta INFO/SVLEN
    if 'SVLEN' in record.info:
        val = record.info['SVLEN']
        return abs(val[0]) if isinstance(val, (list, tuple)) else abs(val)
    # För explicita sekvenser (standard indels)
    return max(len(record.ref), len(record.alts[0]))

def calculate_overlap(chrom, start, end, bed_tabix):
    """Hittar totalt antal överlappande baser i BED-filen."""
    overlap_bp = 0
    try:
        for region in bed_tabix.fetch(chrom, start, end):
            # region format: "chrom\tstart\tend"
            b_start = int(region.split('\t')[1])
            b_end = int(region.split('\t')[2])
            
            # Beräkna skärningen mellan [start, end] och [b_start, b_end]
            o_start = max(start, b_start)
            o_end = min(end, b_end)
            
            if o_start < o_end:
                overlap_bp += (o_end - o_start)
    except ValueError:
        # Chromosome finns inte i BED-filen
        pass
    return overlap_bp

def main(vcf_in_path, bed_gz_path, vcf_out_path):
    vcf_in = pysam.VariantFile(vcf_in_path)
    bed = pysam.TabixFile(bed_gz_path)
    
    # Lägg till ny header för INFO-fältet
    vcf_in.header.info.add('STR_PERCENT', '1', 'Float', 'Percentage of indel overlapping with STR (0-100)')
    
    with pysam.VariantFile(vcf_out_path, 'w', header=vcf_in.header) as vcf_out:
        for record in vcf_in:
            # Manta startar ofta vid POS, och slutar vid POS + längd
            # Vi konverterar till 0-baserad start/end för sökning
            v_start = record.start
            v_len = get_indel_length(record)
            v_end = v_start + v_len
            
            if v_len > 0:
                overlap_bp = calculate_overlap(record.chrom, v_start, v_end, bed)
                percent = round((overlap_bp / v_len) * 100, 2)
                
                if percent > 0:
                    record.info['STR_PERCENT'] = min(percent, 100.0)
            
            vcf_out.write(record)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
