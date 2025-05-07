#!/bin/python3

import xlsxwriter
import datetime
import sys
import logging

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
sample_name = snakemake.input.arriba.split("/")[-1].split("_")[0]
logging.debug(f"Loading arriba results: {snakemake.input.arriba=}")
arriba_table = {"headers": [], "data": []}
with open(snakemake.input.arriba, "r") as arriba_tsv:
    for lline in arriba_tsv:
        if lline.startswith("#"):
            [arriba_table["headers"].append({"header": column}) for column in lline[1:].strip().split("\t")]
        else:
            arriba_table["data"].append(lline.strip().split("\t"))


logging.debug(f"Loading fusioncatcher results: {snakemake.input.fusioncatcher=}")
fusioncatcher_table = {"headers": [], "data": []}
with open(snakemake.input.fusioncatcher, "r") as fusioncatcher_tsv:
    first_row = True
    for lline in fusioncatcher_tsv:
        line = lline.strip().split("\t")
        if first_row:
            [fusioncatcher_table["headers"].append({"header": column}) for column in line]
            first_row = False
        else:
            fusioncatcher_table["data"].append(line)


logging.debug(f"Loading starfusion input: {snakemake.input.star_fusion=}")
starfusion_table = {"headers": [], "data": []}
with open(snakemake.input.star_fusion, "r") as starfusion_tsv:
    for lline in starfusion_tsv:
        if lline.startswith("#"):
            [starfusion_table["headers"].append({"header": column}) for column in lline[1:].strip().split("\t")]
        else:
            starfusion_table["data"].append(lline.strip().split("\t"))


logging.debug(f"Loading dux4 counts: {snakemake.input.dux4_igh_counts=}")
with open(snakemake.input.dux4_igh_counts, "r") as counts_txt:
    lline = counts_txt.readline()
    dux_calls = lline.strip().split("\t")[0]

logging.debug(f"Loading dux4 hits: {snakemake.input.dux4_igh_calls=}")
dux4_table = {"headers": fusioncatcher_table["headers"], "data": []}
with open(snakemake.input.dux4_igh_calls, "r") as calls_txt:
    for lline in calls_txt:
        if lline.strip() != "none":
            line = lline.strip().split("\t")
            dux4_table["data"].append(line)


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
worksheet_overview.write(13, 0, "Number of hits: " + dux_calls)

i = 15
column_end = convert_columns_to_letter(len(dux4_table["headers"]))
if len(dux4_table["data"]) > 0:
    table_area = "A" + str(i) + ":" + column_end + str(len(dux4_table["data"]) + i)
else:
    table_area = "A" + str(i) + ":" + column_end + str(i + 1)
worksheet_overview.add_table(
    table_area, {"columns": dux4_table["headers"], "data": dux4_table["data"], "style": "Table Style Light 1"}
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
