#!/usr/bin/env python

from export_to_xlsx_create_tables import *
import xlsxwriter
import datetime
from pysam import VariantFile
import sys
import logging
import os
from collections import Counter

logging.basicConfig(
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
)

logging.info(f"Using xlsxwriter version: {xlsxwriter.__version__} from {xlsxwriter.__file__}")


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
    genes = []
    if filepath and os.path.exists(filepath):
        try:
            with open(filepath, 'r') as f:
                genes = [line.strip() for line in f if line.strip()]
            logging.info(f"Loaded {len(genes)} target genes from {filepath}")
        except Exception as e:
            logging.error(f"Could not load gene list: {e}")
    return genes


def format_manta_table(table, sample_name, format_2dec):
    """
    Pre-processes the table data to:
    1. Prepend 'Sample Name' to all rows.
    2. Convert any frequency columns to exactly two decimal places.
    """
    if not table or "headers" not in table:
        return
    # Check if we already formatted this table (prevents double-adding)
    if len(table["headers"]) > 0 and table["headers"][0].get("header") == "Sample Name":
        return

    table["headers"] = [{"header": "Sample Name"}] + table["headers"]

    # Identify which columns are frequencies
    freq_indices = []
    for i, h in enumerate(table["headers"]):
        header_str = str(h.get("header", "")).lower()
        if "freq" in header_str:
            h["format"] = format_2dec
            freq_indices.append(i)

    # Process data rows
    for row_idx in range(len(table["data"])):
        row = table["data"][row_idx]
        new_row = [sample_name] + row
        for i in freq_indices:
            try:
                # Force to float and round. The Excel format_2dec ensures it displays as 0.XX
                new_row[i] = round(float(new_row[i]), 2)
            except (ValueError, TypeError, IndexError):
                pass  # Keep original if it's "NA" or empty
        table["data"][row_idx] = new_row


def apply_compact_formatting(worksheet, headers, data):
    """
    Scans the data and applies sets column widths based on header name and content.
    """
    for col_idx, h in enumerate(headers):
        header_str = str(h.get("header", ""))
        max_len = len(header_str)

        # Sample the first 100 rows to find max content length
        for row in data[:100]:
            if col_idx < len(row) and row[col_idx] is not None:
                max_len = max(max_len, len(str(row[col_idx])))

        # Apply compact constraints based on column type
        lower_header = header_str.lower()
        if "freq" in lower_header or "af" in lower_header or "chr" in lower_header:
            final_width = 6
        elif "pos" in lower_header or "length" in lower_header or "sample" in lower_header:
            final_width = 12
        elif "mantaid" in lower_header:
            final_width = 16
        elif "gene" in lower_header or "sample" in lower_header:
            final_width = 20
        else:
            final_width = min(max_len + 2, 35)  # Allow slightly longer text cols but cap at 35

        worksheet.set_column(col_idx, col_idx, final_width)


def write_target_summary(worksheet, start_row, title, table_data, sample_name, format_heading):
    """
    Helper function to write a summary table on the Overview sheet
    containing ONLY the variants where 'In Target Panel' is 'Yes' and manta_n_AF <= 0.2.
    """
    if not table_data or "headers" not in table_data:
        return start_row

    headers = table_data["headers"]
    data = table_data["data"]

    # 1. Find the target column
    target_col_idx = -1
    svdb_col_idx = -1
    for idx, h in enumerate(headers):
        if h.get("header") == "In Target Panel":
            target_col_idx = idx
        elif h.get("header") == "manta_N_AF":
            svdb_col_idx = idx

    if target_col_idx == -1:
        return start_row  # Target column not found, skip

    # 2. Extract target data
    target_data = []
    for row in data:
        if row[target_col_idx] == "Yes":
            keep_variant = True

            if svdb_col_idx != -1:
                af_val = row[svdb_col_idx]
                if af_val is not None and str(af_val).strip() != "":
                    try:
                        # Convert string/int/float to float. If > 0.2, discard it.
                        if float(af_val) > 0.2:
                            keep_variant = False
                    except (ValueError, TypeError):
                        pass

            if keep_variant:
                target_data.append(row)

    # Write the title
    worksheet.write(start_row, 0, title, format_heading)
    table_start_idx = start_row + 2  # Skip a line after title

    if not target_data:
        worksheet.write(table_start_idx, 0, "Inga target-varianter hittades för denna typ.")
        return table_start_idx + 3

    # Calculate Excel table coordinates
    column_end = convert_columns_to_letter(len(headers))
    excel_start_row = table_start_idx + 1
    excel_end_row = excel_start_row + len(target_data)
    table_area = f"A{excel_start_row}:{column_end}{excel_end_row}"

    # Draw the table
    worksheet.add_table(
        table_area,
        {"columns": headers, "data": target_data, "style": "Table Style Light 1"},
    )

    # Apply compact formatting to the summary table
    apply_compact_formatting(worksheet, headers, target_data)
    # Return the next available row (with some margin)
    return table_start_idx + len(target_data) + 4


def create_sheet(workbook, sheet_name, title, sample_name, filter_flags, table_data, set_cols=None):
    if not table_data or "headers" not in table_data:
        return None

    worksheet = workbook.add_worksheet(sheet_name)
    format_heading = workbook.add_format({"bold": True, "font_size": 18})

    # Set column widths if specified
    if set_cols:
        for col_range, width in set_cols.items():
            worksheet.set_column(col_range, width)

    worksheet.write("A1", title, format_heading)
    worksheet.write("A3", "Sample: " + str(sample_name))
    worksheet.write("A5", "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))
    worksheet.write("A6", "MaxDepth calls for regions with depth > 3x median are only included if they have PR or SR support >= 5% and manta_N_OCC = 0")
    row_offset = 7
    if "Deletions" in sheet_name:
        worksheet.write("A7", "Calls have to be longer than 100 bp to be included.")
        row_offset = 8

    headers = table_data["headers"]
    data = table_data["data"]

    # 1. Find columns
    svdb_col_idx = -1
    for idx, header_dict in enumerate(headers):
        if header_dict.get("header") == "manta_N_AF":
            svdb_col_idx = idx

    # xlsxwriter's add_table requires 1-based Excel coordinates (e.g., A7:K20)
    column_end = convert_columns_to_letter(len(headers))
    end_row = len(data) + row_offset if len(data) > 0 else row_offset + 1
    table_area = f"A{row_offset}:{column_end}{end_row}"

    # Create the table with default styling. This automatically applies an autofilter.
    worksheet.add_table(
        table_area,
        {"columns": headers, "data": data, "style": "Table Style Light 1"},
    )

    apply_compact_formatting(worksheet, headers, data)

    # 2. Hide rows with High manta_N_AF (but leave target genes visible)
    if svdb_col_idx != -1 and len(data) > 0:
        for i, row_data in enumerate(data):
            excel_row_index = row_offset + i
            try:
                if float(row_data[svdb_col_idx]) > 0.2:
                    worksheet.set_row(excel_row_index, options={'hidden': True})
            except (ValueError, TypeError):
                pass

    return worksheet


def filter_complex_and_junk_bnd(table, max_sites=4):
    if not table or "headers" not in table or not table["data"]:
        return table

    headers = table["headers"]
    data = table["data"]

    chr_idx = -1
    partner_idx = -1
    id_idx = -1

    for i, h in enumerate(headers):
        header_name = str(h.get("header", "")).lower()
        if header_name == "chr":
            chr_idx = i
        elif header_name in {"breakend", "alt"}:
            partner_idx = i
        elif header_name == "mantaid":
            id_idx = i

    if chr_idx == -1 or id_idx == -1:
        logging.warning("Could not filter BND table, missing Chr or MantaID")
        return table

    id_counts = Counter(str(row[id_idx]) for row in data)
    junk_keywords = ["chrun", "random", "gl0", "ki2"]

    filtered_data = []
    for row in data:
        chr_val = str(row[chr_idx]).lower()
        partner_val = str(row[partner_idx]).lower() if partner_idx != -1 else ""
        manta_id = str(row[id_idx])

        is_junk_contig = (
            any(k in chr_val for k in junk_keywords) or
            any(k in partner_val for k in junk_keywords)
        )

        has_too_many_sites = id_counts[manta_id] > max_sites
        # use >= if you want 4-or-more removed

        if not is_junk_contig and not has_too_many_sites:
            filtered_data.append(row)

    table["data"] = filtered_data
    logging.info("Top repeated BND IDs: %s", id_counts.most_common(10))
    return table


def filter_maxdepth_by_support(tables_dict, min_support=0.05):
    """
    Removes variants flagged with 'MaxDepth' UNLESS their PR or SR frequency is >= min_support.
    """
    for sv_type in ["bnd", "del", "dup", "ins"]:
        table_data = tables_dict[sv_type]
        if not table_data["data"]:
            continue

        headers = [h["header"] for h in table_data["headers"]]
        try:
            ann_idx = headers.index("Annotation")
            pr_idx = headers.index("Paired-read freq")
            sr_idx = headers.index("Spanning-read freq")
            occ_idx = headers.index("manta_N_OCC")
        except ValueError:
            continue

        filtered_data = []
        for row in table_data["data"]:
            annotation = str(row[ann_idx])

            # Check if variant is flagged with MaxDepth
            if "MaxDepth" in annotation:
                # Safely extract frequency floats (they might be "" if no PR/SR was found)
                pr_freq = row[pr_idx] if isinstance(row[pr_idx], float) else 0.0
                sr_freq = row[sr_idx] if isinstance(row[sr_idx], float) else 0.0

                # Keep it only if support meets our threshold
                if (pr_freq >= min_support or sr_freq >= min_support) and row[occ_idx] == 0:
                    filtered_data.append(row)
            else:
                # Keep all other non-MaxDepth rows normally
                filtered_data.append(row)

        tables_dict[sv_type]["data"] = filtered_data
    return tables_dict


""" MAIN EXECUTION """

# 1. Prepping data
logging.info(f"Prepping data, such as loading {snakemake.input.manta}=")
sample_name = snakemake.output.xlsx.split("/")[-1].split(".manta.xlsx")[0]

filter_flags = ["MinQUAL", "MinGQ", "MinSomaticScore", "Ploidy", "MaxMQ0Frac", "NoPairSupport", "SampleFT", "HomRef"]

# Load target genes for easy filtering in Excel
target_genes = []
if hasattr(snakemake.input, 'target_genes'):
    target_genes = load_target_genes(snakemake.input.target_genes)

manta_tables_full = create_manta_tables(snakemake.input.manta, filter_flags, target_genes=target_genes)
manta_tables_full = filter_maxdepth_by_support(manta_tables_full, 0.05)
manta_tables_full["bnd"] = filter_complex_and_junk_bnd(manta_tables_full["bnd"], max_sites=4)


# 2. Creating xlsx workbook
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
logging.info(f"Creating xlsx workbook {snakemake.output.xlsx}=")

format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_bold = workbook.add_format({"bold": True, "text_wrap": True})
format_overview_title = workbook.add_format({"bold": True, "font_size": 16})
format_2dec = workbook.add_format({"num_format": "0.00"})

for sv_key in ["del", "ins", "dup", "bnd"]:
    format_manta_table(manta_tables_full[sv_key], sample_name, format_2dec)

panel_tables_dict = {}
for vcf in snakemake.input.vcfs_bed:
    panel = vcf.split(".")[-3]
    panel_tables = create_manta_tables(vcf, filter_flags, target_genes=target_genes)
    format_manta_table(panel_tables["bnd"], sample_name, format_2dec)
    panel_tables_dict[panel] = panel_tables

worksheet_overview = workbook.add_worksheet("Overview")

# 3. Create Data Sheets
create_sheet(
    workbook, "Deletions", "Deletions found by Manta",
    sample_name, filter_flags, manta_tables_full["del"],  {"B:C": 12, "E:E": 12}
)
create_sheet(
    workbook, "Insertions", "Insertions found by Manta",
    sample_name, filter_flags, manta_tables_full["ins"], {"B:B": 12, "F:F": 12}
)
create_sheet(
    workbook, "Duplications", "Duplications found by Manta",
    sample_name, filter_flags, manta_tables_full["dup"], {"B:C": 12, "E:E": 12}
)
create_sheet(
    workbook, "Translocations", "Translocations found by Manta",
    sample_name, filter_flags, manta_tables_full["bnd"], {"B:B": 12, "C:D": 15}
)

# Translocations (Panels from BED)
for vcf in snakemake.input.vcfs_bed:
    panel = vcf.split(".")[-3]
    logging.debug(f"Creating {panel} sheet")
    panel_tables = create_manta_tables(vcf, filter_flags, target_genes=target_genes)
    panel_tables["bnd"] = filter_complex_and_junk_bnd(panel_tables["bnd"], max_sites=4)
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

for vcf in snakemake.input.vcfs_bed:
    panel = vcf.split(".")[-3]
    s_name = "Translocations " + panel.upper()
    worksheet_overview.write_url(row_idx, 0, f"internal:'{s_name}'!A1", string=f"Manta Translocations in {panel.upper()} genes")
    row_idx += 1

if hasattr(snakemake.input, 'all_bed'):
    worksheet_overview.write(row_idx + 4, 0, "ALL bedfile: " + snakemake.input.all_bed)
if hasattr(snakemake.input, 'aml_bed'):
    worksheet_overview.write(row_idx + 5, 0, "AML bedfile: " + snakemake.input.aml_bed)
if target_genes:
    genes_string = ", ".join(target_genes)
    worksheet_overview.write(row_idx + 6, 0,
                             f"Target Genes filter added: {len(target_genes)} genes loaded - [{genes_string}]")

worksheet_overview.write(row_idx + 9, 0,
                         "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))

# -------------------------------------------------------------
# 5. Add Summary Tables for Target Panel == "Yes" on Overview
# -------------------------------------------------------------
row_idx += 12  # Move down below the metadata

worksheet_overview.write(row_idx, 0, "Variants in gene list", format_overview_title)
row_idx += 2

summaries = [
    ("del", "Deletions"),
    ("ins", "Insertions"),
    ("dup", "Duplications"),
    ("bnd", "Translocations")
]

for sv_key, sv_title in summaries:
    row_idx = write_target_summary(
        worksheet_overview,
        row_idx,
        sv_title,
        manta_tables_full[sv_key],
        sample_name,
        format_bold
    )

workbook.set_size(1800, 1200)
workbook.close()
logging.info("All done")
