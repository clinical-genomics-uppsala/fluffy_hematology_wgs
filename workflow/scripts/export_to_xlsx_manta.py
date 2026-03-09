#!/bin/python3
# // src/export_to_xlsx_manta.py : v0.8 : 14:55

from export_to_xlsx_create_tables import *
import xlsxwriter
import datetime
from pysam import VariantFile
import sys
import logging
import os

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

def load_target_genes(filepath):
    """Läs in genlista från fil om den finns."""
    genes = []
    if filepath and os.path.exists(filepath):
        try:
            with open(filepath, 'r') as f:
                genes = [line.strip() for line in f if line.strip()]
            logging.info(f"Loaded {len(genes)} target genes from {filepath}")
        except Exception as e:
            logging.error(f"Could not load gene list: {e}")
    return genes

def create_sheet(workbook, sheet_name, title, sample_name, filter_flags, table_data, set_cols=None):
    """Hjälpfunktion för att skapa en standardflik med tabell."""
    if not table_data or "headers" not in table_data:
        return None

    worksheet = workbook.add_worksheet(sheet_name)
    format_heading = workbook.add_format({"bold": True, "font_size": 18})
    
    # Sätt kolumnbredd om angivet
    if set_cols:
        for col_range, width in set_cols.items():
            worksheet.set_column(col_range, width)
    
    worksheet.write("A1", title, format_heading)
    worksheet.write("A3", "Sample: " + str(sample_name))
    worksheet.write("A5", "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))
    
    row_offset = 7
    if "Deletions" in sheet_name:
         worksheet.write("A6", "Calls have to be longer than 100 bp to be included.")
         row_offset = 8

    headers = table_data["headers"]
    data = table_data["data"]
    
    column_end = ":" + convert_columns_to_letter(len(headers))
    end_row = len(data) + row_offset if len(data) > 0 else row_offset + 1
    table_area = f"A{row_offset}{column_end}{end_row}"

    worksheet.add_table(
        table_area,
        {"columns": headers, "data": data, "style": "Table Style Light 1"},
    )
    return worksheet


""" MAIN EXECUTION """

# 1. Prepping data
logging.info(f"Prepping data, such as loading {snakemake.input.vcf}=")
sample_name = snakemake.output.xlsx.split("/")[-1].split(".manta.xlsx")[0]

filter_flags = ["MinQUAL", "MinGQ", "MinSomaticScore", "Ploidy", "MaxDepth", "MaxMQ0Frac", "NoPairSupport", "SampleFT", "HomRef"]

# Ladda in genlista om den finns
target_genes = []
if hasattr(snakemake.input, 'target_genes'):
    target_genes = load_target_genes(snakemake.input.target_genes)

# Hämta tabellerna (nu hanteras In Target Panel internt av den funktionen om target_genes skickas in)
manta_tables_full = create_manta_tables(snakemake.input.vcf, filter_flags, target_genes=target_genes)

# 2. Creating xlsx workbook
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
logging.info(f"Creating xlsx workbook {snakemake.output.xlsx}=")

format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_bold = workbook.add_format({"bold": True, "text_wrap": True})

# Skapa Overview först så den garanterat blir flik nr 1 i Excel
worksheet_overview = workbook.add_worksheet("Overview")

# 3. Create Data Sheets
create_sheet(workbook, "Deletions", "Deletions found by Manta", sample_name, filter_flags, manta_tables_full["del"], {"B:C": 12, "E:E": 12})
create_sheet(workbook, "Insertions", "Insertions found by Manta", sample_name, filter_flags, manta_tables_full["ins"], {"B:B": 12, "F:F": 12})
create_sheet(workbook, "Duplications", "Duplications found by Manta", sample_name, filter_flags, manta_tables_full["dup"], {"B:C": 12, "E:E": 12})
create_sheet(workbook, "Translocations", "Translocations found by Manta", sample_name, filter_flags, manta_tables_full["bnd"], {"B:B": 12, "C:D": 15})

# Translocations (Panels from BED)
for vcf in snakemake.input.vcfs_bed:
    panel = vcf.split(".")[-3]
    logging.debug(f"Creating {panel} sheet")
    
    # Glöm inte att skicka in target_genes även hit om ni vill ha kolumnen i panel-flikarna
    panel_tables = create_manta_tables(vcf, filter_flags, target_genes=target_genes)
    
    sheet_title = "Translocations in " + panel.upper() + " genes"
    sheet_name = "Translocations " + panel.upper()
    create_sheet(workbook, sheet_name, sheet_title, sample_name, filter_flags, panel_tables["bnd"], {"B:B": 12, "C:D": 15})

# 4. Populate Overview Sheet
logging.debug(f"Populating Overview sheet")
worksheet_overview.write(0, 0, sample_name, format_heading)
worksheet_overview.write(1, 0, "Processing date: " + datetime.datetime.now().strftime("%d %B, %Y"))
worksheet_overview.write(4, 0, "Created by: ")
worksheet_overview.write(4, 4, "Valid from: ")
worksheet_overview.write(5, 0, "Signed by: ")
worksheet_overview.write(5, 4, "Document nr: ")
worksheet_overview.write(7, 0, "Sheets:", format_bold)

# Länkar till de statiska flikarna
row_idx = 8
link_map = [
    ("Deletions", "Manta Deletions"),
    ("Insertions", "Manta Insertions"),
    ("Duplications", "Manta Duplications"),
    ("Translocations", "Manta Translocations or breakpoints")
]

for sheet_name, desc in link_map:
    if workbook.get_worksheet_by_name(sheet_name):
        worksheet_overview.write_url(row_idx, 0, f"internal:'{sheet_name}'!A1", string=desc)
        row_idx += 1

# Länkar till de dynamiska panel-flikarna
for vcf in snakemake.input.vcfs_bed:
    panel = vcf.split(".")[-3]
    s_name = "Translocations " + panel.upper()
    worksheet_overview.write_url(row_idx, 0, f"internal:'{s_name}'!A1", string=f"Manta Translocations in {panel.upper()} genes")
    row_idx += 1

# Metadata om körningen
if hasattr(snakemake.input, 'all_bed'):
    worksheet_overview.write(row_idx + 4, 0, "ALL bedfile: " + snakemake.input.all_bed)
if hasattr(snakemake.input, 'aml_bed'):
    worksheet_overview.write(row_idx + 5, 0, "AML bedfile: " + snakemake.input.aml_bed)
if target_genes:
    genes_string = ", ".join(target_genes)
    worksheet_overview.write(row_idx + 6, 0, f"Target Genes filter added: {len(target_genes)} genes loaded - [{genes_string}]")
    
worksheet_overview.write(row_idx + 9, 0, "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))

workbook.close()
logging.info("All done")