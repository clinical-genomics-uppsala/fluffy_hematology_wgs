#!/bin/python3

import xlsxwriter
import datetime
import sys
import logging
import operator

logging.basicConfig(
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
)


def convert_columns_to_letter(nr_columns):
    # Function to convert number of columns to alphabetical coordinates for xlsx-sheets
    if nr_columns < 27:
        letter = chr(nr_columns + 64)
    elif nr_columns < 703:
        i = int((nr_columns - 1) / 26)
        letter = chr(i + 64) + chr(nr_columns - (i * 26) + 64)
    else:
        logging.error(f"Nr columns has to be less than 703, does not support three letter column-index for tables {nr_columns=}")
        sys.exit()
    return letter


""" Prepping input data """
logging.info("Prepping input data")

logging.info("Loading fusioncatcher genelist if exists")
single_genes = set()
gene_pairs = set()
if snakemake.params.fusioncatcher_genelist:
    with open(snakemake.params.fusioncatcher_genelist, "r") as genelist_txt:
        for lline in genelist_txt:
            entry = lline.strip()
            if not entry:
                continue
            if "|" in entry:
                parts = [p.strip().upper() for p in entry.split("|") if p.strip()]
                if len(parts) == 2:
                    gene_pairs.add(frozenset(parts))
                else:
                    logging.error(f"Invalid fusion pair format: {entry}, only pairs allowed")
                    sys.exit(1)
            else:
                single_genes.add(entry.upper())
        logging.debug(f"Loaded {len(single_genes)=} single genes and {len(gene_pairs)=} gene pairs")
        logging.debug("Genes included in Fusioncatcher short table: " + str(single_genes | gene_pairs))


logging.info(f"Loading arriba results: {snakemake.input.arriba=}")
sample_name = snakemake.input.arriba.split("/")[-1].split("_")[0]
arriba_table = {"headers": [], "data": []}
with open(snakemake.input.arriba, "r") as arriba_tsv:
    for lline in arriba_tsv:
        if lline.startswith("#"):
            [arriba_table["headers"].append({"header": column}) for column in lline[1:].strip().split("\t")]
        else:
            arriba_table["data"].append(lline.strip().split("\t"))


# Loading fusioncatcher results in to three different lists/tables
#     - fusioncatcher_table: the whole table
#     - fusioncatcher_table_short: fusioncatcher results with only genes from single_genes and gene_pairs included, inframe fusions in one list and the rest in another to be able to print with in-frame fusions first
#     - fusioncatcher_dux_table: fusions DUX4::IGH or DUX4::ERG
logging.info(f"Loading fusioncatcher results: {snakemake.input.fusioncatcher=}")
fusioncatcher_table = {"headers": [], "data": []}
fusioncatcher_table_short = {"inframe_rows": [], "the_rest_rows": []}
fusioncatcher_dux_table = {"data": []}
dux_pairs = [frozenset(["DUX4", "ERG"]), frozenset(["DUX4", "IGH@"])]
logging.debug(f"{dux_pairs=}")
with open(snakemake.input.fusioncatcher, "r") as fusioncatcher_tsv:
    first_row = True
    for lline in fusioncatcher_tsv:
        line = lline.strip().split("\t")
        if first_row:
            [fusioncatcher_table["headers"].append({"header": column}) for column in line]
            predicted_effect_idx = line.index("Predicted_effect")
            first_row = False
        else:
            fusioncatcher_table["data"].append(line)
            g1 = line[0].strip().upper()
            g2 = line[1].strip().upper()
            predicted_effect = line[predicted_effect_idx].strip().lower()
            if frozenset([g1, g2]) in dux_pairs:
                fusioncatcher_dux_table["data"].append(line)
            elif g1 in single_genes or g2 in single_genes or frozenset([g1, g2]) in gene_pairs:
                if "in-frame" in predicted_effect:
                    fusioncatcher_table_short["inframe_rows"].append(line)
                else:
                    fusioncatcher_table_short["the_rest_rows"].append(line)
logging.debug(f"Number of dux4-pairs found: {len(fusioncatcher_dux_table['data'])=}")
fusioncatcher_table_short["data"] = fusioncatcher_table_short["inframe_rows"] + fusioncatcher_table_short["the_rest_rows"]
logging.debug(f"Number of fusions in Fusioncatcher short table: {len(fusioncatcher_table_short['data'])=}")


logging.info(f"Loading starfusion input: {snakemake.input.star_fusion=}")
starfusion_table = {"headers": [], "data": []}
with open(snakemake.input.star_fusion, "r") as starfusion_tsv:
    for lline in starfusion_tsv:
        if lline.startswith("#"):
            [starfusion_table["headers"].append({"header": column}) for column in lline[1:].strip().split("\t")]
        else:
            starfusion_table["data"].append(lline.strip().split("\t"))


""" Creating xlsx file """
logging.info(f"Creating xlsx-workbbok {snakemake.output.xlsx=}")
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
worksheet_overview = workbook.add_worksheet("Overview")
worksheet_arriba = workbook.add_worksheet("Arriba")
worksheet_fusioncatcher = workbook.add_worksheet("Fusioncatcher")
worksheet_starfusion = workbook.add_worksheet("StarFusion")

format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_bold = workbook.add_format({"bold": True, "text_wrap": True})

# Overview sheet
logging.debug(f"Creating Overview sheet")
worksheet_overview.write(0, 0, sample_name, format_heading)
worksheet_overview.write(1, 0, "Processing date: " + datetime.datetime.now().strftime("%d %B, %Y"))

worksheet_overview.write(4, 0, "Created by: ")
worksheet_overview.write(4, 4, "Valid from: ")
worksheet_overview.write(5, 0, "Signed by: ")
worksheet_overview.write(5, 4, "Document nr: ")

worksheet_overview.write(7, 0, "Sheets:", format_bold)
worksheet_overview.write_url(8, 0, "internal:'Arriba'!A1", string="Arriba fusions")
worksheet_overview.write_url(9, 0, "internal:'Fusioncatcher'!A1", string="Fusioncatcher results")
worksheet_overview.write_url(10, 0, "internal:'StarFusion'!A1", string="StarFusion results")

worksheet_overview.write(12, 0, "DUX4-IGH and DUX4-ERG hits from Fusioncatcher", format_bold)
worksheet_overview.write(13, 0, "Number of hits: " + str(len(fusioncatcher_dux_table["data"])))
worksheet_overview.write(
    15,
    0,
    "Genes included in Fusioncatcher short table: "
    + ", ".join(sorted(list(single_genes) + ["|".join(sorted(pair)) for pair in gene_pairs])),
)

i = 18
logging.debug(f"Creating the DUX4 table")
column_end = convert_columns_to_letter(len(fusioncatcher_table["headers"]))
if len(fusioncatcher_dux_table["data"]) > 0:
    table_area_dux = "A" + str(i) + ":" + column_end + str(len(fusioncatcher_dux_table["data"]) + i)
else:
    table_area_dux = "A" + str(i) + ":" + column_end + str(i + 1)

worksheet_overview.add_table(
    table_area_dux,
    {"columns": fusioncatcher_table["headers"], "data": fusioncatcher_dux_table["data"], "style": "Table Style Light 1"},
)


i = i + len(fusioncatcher_dux_table["data"]) + 3
logging.debug(f"Creating fusioncatcher short table")
if len(fusioncatcher_table_short["data"]) > 0:
    table_area_short = "A" + str(i) + ":" + column_end + str(len(fusioncatcher_table_short["data"]) + i)
else:
    table_area_short = "A" + str(i) + ":" + column_end + str(i + 1)

worksheet_overview.add_table(
    table_area_short,
    {"columns": fusioncatcher_table["headers"], "data": fusioncatcher_table_short["data"], "style": "Table Style Light 1"},
)

# Arriba sheet
logging.debug(f"Creating Arriba sheet")
worksheet_arriba.set_column("E:F", 12)

worksheet_arriba.write("A1", "Fusions detected by Arriba", format_heading)
worksheet_arriba.write("A3", "Sample: " + str(sample_name))

i = 5
column_end = convert_columns_to_letter(len(arriba_table["headers"]))
if len(arriba_table["data"]) > 0:
    table_area = "A" + str(i) + ":" + column_end + str(len(arriba_table["data"]) + i)
else:
    table_area = "A" + str(i) + ":" + column_end + str(i + 1)

worksheet_arriba.add_table(
    table_area, {"columns": arriba_table["headers"], "data": arriba_table["data"], "style": "Table Style Light 1"}
)

# Fusioncatcher sheet
logging.debug(f"Creating Fusioncatcher sheet")
worksheet_fusioncatcher.set_column("C:C", 12)
worksheet_fusioncatcher.set_column("I:J", 12)

worksheet_fusioncatcher.write("A1", "Fusions detected by Fusioncatcher", format_heading)
worksheet_fusioncatcher.write("A3", "Sample: " + str(sample_name))

i = 5
column_end = convert_columns_to_letter(len(fusioncatcher_table["headers"]))
if len(fusioncatcher_table["data"]) > 0:
    table_area = "A" + str(i) + ":" + column_end + str(len(fusioncatcher_table["data"]) + i)
else:
    table_area = "A" + str(i) + ":" + column_end + str(i + 1)

worksheet_fusioncatcher.add_table(
    table_area, {"columns": fusioncatcher_table["headers"], "data": fusioncatcher_table["data"], "style": "Table Style Light 1"}
)

# StarFusion sheet
logging.debug(f"Creating Starfusion sheet")
worksheet_starfusion.set_column("A:A", 18)
worksheet_starfusion.set_column("F:J", 15)

worksheet_starfusion.write("A1", "Fusions detected by StarFusion", format_heading)
worksheet_starfusion.write("A3", "Sample: " + str(sample_name))

i = 5
column_end = convert_columns_to_letter(len(starfusion_table["headers"]))
if len(starfusion_table["data"]) > 0:
    table_area = "A" + str(i) + ":" + column_end + str(len(starfusion_table["data"]) + i)
else:
    table_area = "A" + str(i) + ":" + column_end + str(i + 1)

worksheet_starfusion.add_table(
    table_area, {"columns": starfusion_table["headers"], "data": starfusion_table["data"], "style": "Table Style Light 1"}
)

workbook.close()
logging.info(f"All done!")
