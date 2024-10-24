#!/bin/python3

import sys
import csv
import re


# process input as pairs sample, fastq1 from units.tsv
sample_order_duplicates = []
sample_order_index = ["sample", "s_index", "lab_id", "type"]
for sample, type, fastq_path in snakemake.params.filelist:
    fastq = fastq_path.split("/")[-1]
    lab_id = fastq.split("_")[1]
    s_pattern = re.compile("_S([0-9]+)_")
    # In case of missing S-index in fastq1-filename set s_index to 99 (last)
    try:
        s_index = int(s_pattern.search(fastq).group(1))
    except AttributeError:
        s_index = 99
    # If same sample sequenced twice use latest runs s_index
    for old_sample, old_s, old_lab, old_type in sample_order_duplicates:
        if old_sample == sample + "_" + type and old_lab == lab_id and s_index == 99:
            s_index = old_s
    sample_order_duplicates.append([sample + "_" + type, s_index, lab_id, type])

# Remove duplicates and order based on s_index in fastq1 filename
sample_order = [list(x) for x in set(tuple(x) for x in sample_order_duplicates)]
sample_order.sort(key=lambda x: int(x[1]))

with open(snakemake.output.replacement_dna, "w+") as replacement_dna, \
     open(snakemake.output.replacement_rna, "w+") as replacement_rna, \
     open(snakemake.output.order_dna, "w+") as order_dna, \
     open(snakemake.output.order_rna, "w+") as order_rna, \
     open(snakemake.output.dnanumber, "w+") as dna_table, \
     open(snakemake.output.rnanumber, "w+") as rna_table:

    order_dna.write("\t".join(["Sample Order", "Pedegree ID", "DNA number"]) + "\n")
    order_rna.write("\t".join(["Sample Order", "Pedegree ID", "RNA number"]) + "\n")
    dna_table.write("\t".join(["Sample", "ped_id", "dna_number"]) + "\n")
    rna_table.write("\t".join(["Sample", "ped_id", "rna_number"]) + "\n")

    i = 1
    j = 1
    for sample_line in sample_order:
        replacement_line = [sample_line[sample_order_index.index("sample")]]
        order_line = [sample_line[sample_order_index.index("sample")], sample_line[sample_order_index.index("lab_id")]]
        if (
            sample_line[sample_order_index.index("type")].lower() == "t"
            or sample_line[sample_order_index.index("type")].lower() == "n"
        ):
            replacement_dna.write("\t".join(replacement_line + ["sample_" + str(f"{i:03}")]) + "\n")
            order_dna.write("\t".join(["sample_" + str(f"{i:03}")] + order_line) + "\n")
            dna_table.write("\t".join(["sample_" + str(f"{i:03}")] + order_line) + "\n")
            i += 1
        elif sample_line[sample_order_index.index("type")].lower() == "r":
            replacement_rna.write("\t".join(replacement_line + ["sample_" + str(f"{j:03}")]) + "\n")
            order_rna.write("\t".join(["sample_" + str(f"{j:03}")] + order_line) + "\n")
            rna_table.write("\t".join(["sample_" + str(f"{j:03}")] + order_line) + "\n")
            j += 1
