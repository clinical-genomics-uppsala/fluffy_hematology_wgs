#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import sys
from pysam import VariantFile
import logging

logging.basicConfig(
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
)

input_file = VariantFile(snakemake.input.vcf)
analysis = snakemake.params.analysis
output_file = snakemake.output.tsv


csq_index = []
for x in input_file.header.records:
    if "CSQ" in str(x):
        csq_index = str(x).split("Format: ")[1].strip().strip('">').split("|")

outlines = []
if len(list(input_file.header.samples)) == 1:
    logging.debug("One sample identified. " + f"{list(input_file.header.samples)=}")
    outlines.append(
        [
            "#SAMPLE",
            "CHROMOSOME",
            "POSITION",
            "REFERENCE",
            "ALTERNATIVE",
            "ALLELEFREQUENCY",
            "DEPTH",
            "GENE",
            "FILTERFLAG",
            "TRANSCRIPT",
            "NUCLEOTIDESEQ",
            "PROTEIN",
            "PROTEINSEQ",
            "CONSEQUENCE",
        ]
    )

elif len(list(input_file.header.samples)) == 2:
    logging.debug("Two sample identified. " + f"{list(input_file.header.samples)=}")
    outlines.append(
        [
            "#SAMPLE",
            "CHROMOSOME",
            "POSITION",
            "REFERENCE",
            "ALTERNATIVE",
            "ALLELEFREQUENCY",
            "NORMAL_AF",
            "DEPTH",
            "NORMAL_DP",
            "GENE",
            "FILTERFLAG",
            "TRANSCRIPT",
            "NUCLEOTIDESEQ",
            "PROTEIN",
            "PROTEINSEQ",
            "CONSEQUENCE",
        ]
    )
else:
    logging.error(snakemake.input.vcf + " contains more than two samples in header.")
    sys.exit()


logging.info("Starting to iterate through calls in " + snakemake.input.vcf)
for row in input_file.fetch():
    csq = row.info["CSQ"][0].split("|")
    gene = csq[csq_index.index("SYMBOL")]

    transcript_name = csq[csq_index.index("HGVSc")].split(":")[0]
    transcript_change = ""
    if len(csq[csq_index.index("HGVSc")].split(":")) > 1:
        transcript_change = csq[csq_index.index("HGVSc")].split(":")[1]

    protein_name = csq[csq_index.index("HGVSp")].split(":")[0]
    protein_change = ""
    if len(csq[csq_index.index("HGVSp")].split(":")) > 1:
        protein_change = csq[csq_index.index("HGVSp")].split(":")[1]

    consequence = csq[csq_index.index("Consequence")]

    tumor_sample = [x for x in row.samples if x.endswith("_T")][0]
    af = row.samples[tumor_sample]["AF"][0]
    dp = row.samples[tumor_sample]["DP"]
    flag = ",".join(row.filter.keys())
    outline = [
        tumor_sample.split("_")[0],
        row.chrom,
        row.pos,
        row.ref,
        row.alts[0],
        af,
        dp,
        gene,
        flag,
        transcript_name,
        transcript_change,
        protein_name,
        protein_change,
        consequence,
    ]

    if len(list(input_file.header.samples)) == 2:
        normal_sample = [x for x in row.samples if x.endswith("_N")][0]
        af_normal = row.samples[normal_sample]["AF"][0]
        dp_normal = row.samples[normal_sample]["DP"]
        outline = [
            tumor_sample.split("_")[0],
            row.chrom,
            row.pos,
            row.ref,
            row.alts[0],
            af,
            af_normal,
            dp,
            dp_normal,
            gene,
            flag,
            transcript_name,
            transcript_change,
            protein_name,
            protein_change,
            consequence,
        ]

    outlines.append(outline)

logging.info("Writing results to " + output_file)
with open(output_file, "wt") as tsv:
    tsv_writer = csv.writer(tsv, delimiter="\t")
    tsv_writer.writerows(outlines)
