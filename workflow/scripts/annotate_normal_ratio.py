#!/bin/python3

from pysam import VariantFile
import logging

logging.basicConfig(
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
)

logging.info("Starting to process " + snakemake.input.vcf)

vcf_in = VariantFile(snakemake.input.vcf)

new_header = vcf_in.header
new_header.info.add("N_ratio", "1", "Float", "Number of times larger normal AF is than tumour sample.")
logging.debug("Creating header for " + snakemake.output.vcf)
vcf_out = VariantFile(snakemake.output.vcf, "w", header=new_header)

sample_tumor = [x for x in list(vcf_in.header.samples) if x.endswith("_T")][0]
sample_normal = [x for x in list(vcf_in.header.samples) if x.endswith("_N")][0]
logging.debug(f"{sample_tumor=}, {sample_normal=}")

logging.info("Starting to iterate through input.")
for record in vcf_in.fetch():
    try:
        n_freq = float(record.samples[sample_normal].get("AF")[0])
    except TypeError:
        n_freq = 0

    try:
        t_freq = float(record.samples[sample_tumor].get("AF")[0])
    except TypeError:
        t_freq = 0

    if n_freq == 0:
        logging.warning("Normal freq 0 or missing, setting ratio to 0")
        ratio = 0
    elif t_freq == 0:
        logging.warning("Tumor freq 0 or missing, setting ratio to 100")
        ratio = 100
    else:
        ratio = float(record.samples[sample_normal].get("AF")[0]) / float(record.samples[sample_tumor].get("AF")[0])

    logging.debug(f"{record.samples[sample_normal].get('AF')[0]=}, {record.samples[sample_tumor].get('AF')[0]=}")
    logging.debug(f"{ratio=}")
    record.info["N_ratio"] = ratio

    vcf_out.write(record)
