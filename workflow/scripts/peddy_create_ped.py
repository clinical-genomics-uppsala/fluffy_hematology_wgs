#!/usr/bin/env python
# coding: utf-8


import sys


input_file = snakemake.input[0]


with open(input_file, "r") as samplesheet:
    header_line = samplesheet.readline().strip().split("\t")
    for lline in samplesheet:
        line = lline.strip().split("\t")
        if line[header_line.index("sex")] == "M":
            sex = "1"
        elif line[header_line.index("sex")] == "K":
            sex = "2"
        else:
            sex = "0"
        with open("qc/peddy/" + line[header_line.index("sample")] + ".peddy.fam", "w+") as pedfile:
            pedfile.write(
                "\t".join([line[header_line.index("sample")], line[header_line.index("sample")] + "_T", "0", "0", sex, "-9"])
                + "\n"
            )
