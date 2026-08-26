#!/usr/bin/env python

from export_to_xlsx_create_tables import *
import xlsxwriter
from xlsxwriter.utility import xl_col_to_name
import datetime
import sys
import logging
import os
from collections import Counter
from manta_maxdepth_rescue import *

MAX_OVERVIEW_NORMAL_AF = 0.2

COLUMN_WIDTHS = {  # This is also used to order columns on Overview
    "Sample Name": 18,
    "Chr": 6,
    "Pos": 10,
    "EndPos": 10,
    "Ref": 12,
    "Alt": 24,
    "SV Length": 11,
    "MantaID": 12,
    "Depth": 6,
    "In Target Panel": 6,
    "manta_N_OCC": 6,
    "manta_T_OCC": 6,
    "manta_N_AF": 8,
    "manta_T_AF": 8,
    "STR %": 4,
    "Paired-read freq": 6,
    "Spanning-read freq": 6,
    "Annotation": 16,
    "Genes": 26,
    "Details": 20,
    "Hom Length": 11,
    "Hom Sequence": 18,
    "BreakEnd": 18,
    "BND Event ID": 30,
    "MaxDepth Rescue": 6,
}

logging.basicConfig(
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
)

logging.info(f"Using xlsxwriter version: {xlsxwriter.__version__} from {xlsxwriter.__file__}")


def add_maxdepth_rescue_row_color(
    worksheet,
    workbook,
    headers,
    data,
    row_offset,
):
    columns = column_indexes({"headers": headers})
    rescue_idx = columns.get("maxdepth rescue")

    if rescue_idx is None or not data:
        return

    rescue_col = xl_col_to_name(rescue_idx)
    last_col = xl_col_to_name(len(headers) - 1)

    first_data_row = row_offset + 1
    last_data_row = row_offset + len(data)
    data_range = f"A{first_data_row}:{last_col}{last_data_row}"

    rescue_format = workbook.add_format({
        "bg_color": "#DDEBF7",
        "font_color": "#1F4E78",
        "bold": True,
    })

    worksheet.conditional_format(data_range, {
        "type": "formula",
        "criteria": (
            f'=${rescue_col}{first_data_row}="Yes"'
        ),
        "format": rescue_format,
        "stop_if_true": True,
    })


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
    if not filepath:
        logging.warning("No target gene list supplied; target-panel summaries will be omitted.")
        return []

    try:
        with open(filepath, encoding="utf-8") as gene_file:
            genes = [line.strip() for line in gene_file if line.strip()]
    except FileNotFoundError as error:
        raise FileNotFoundError(
            f"Target gene list does not exist: {filepath}"
        ) from error
    except OSError as error:
        raise OSError(
            f"Could not read target gene list: {filepath}"
        ) from error

    if not genes:
        raise ValueError(f"Target gene list is empty: {filepath}")

    logging.info("Loaded %d target genes from %s", len(genes), filepath)
    return genes


def format_manta_table(table, sample_name, format_2dec):
    """
    Pre-processes the table data to add 'Sample Name' to all rows and convert freq columns to two decimal places.
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
    """Apply stable, field-specific widths to an exported table."""
    for col_idx, header in enumerate(headers):
        name = header.get("header", "")
        width = COLUMN_WIDTHS.get(name, 14)
        worksheet.set_column(col_idx, col_idx, width)


def add_manta_af_row_colors(worksheet, workbook, headers, data, row_offset):
    """Color and hide rows with high normal-sample allele frequency."""
    columns = column_indexes({"headers": headers})
    af_idx = columns.get("manta_n_af")

    if af_idx is None or not data:
        return

    af_col = xl_col_to_name(af_idx)
    last_col = xl_col_to_name(len(headers) - 1)
    first_data_row = row_offset + 1
    last_data_row = row_offset + len(data)
    data_range = f"A{first_data_row}:{last_col}{last_data_row}"

    format_orange = workbook.add_format({"bg_color": "#ffd280"})

    worksheet.conditional_format(data_range, {
        "type": "formula",
        "criteria": (
            f"=AND(ISNUMBER(${af_col}{first_data_row}),"
            f"${af_col}{first_data_row}>{MAX_OVERVIEW_NORMAL_AF})"
        ),
        "format": format_orange,
    })


def create_sheet(workbook, sheet_name, title, sample_name, filter_flags, table_data, set_cols=None):
    if not table_data or "headers" not in table_data:
        return None

    worksheet = workbook.add_worksheet(sheet_name)
    format_heading = workbook.add_format({"bold": True, "font_size": 18})

    worksheet.write("A1", title, format_heading)
    worksheet.write("A3", "Sample: " + str(sample_name))
    worksheet.write("A5", "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))
    if "Translocations" in sheet_name:
        worksheet.write("A6", "MaxDepth calls passing rescue criteria are included.")

    row_offset = 7
    if "Deletions" in sheet_name:
        worksheet.write("A7", "Deletions have to be longer than 100 bp to be included.")
        row_offset = 8

    headers = table_data["headers"]
    data = table_data["data"]

    # 1. Find columns
    columns = column_indexes(table_data)
    svdb_col_idx = columns.get("manta_n_af")

    # xlsxwriter's add_table requires 1-based Excel coordinates (e.g., A7:K20)
    column_end = convert_columns_to_letter(len(headers))
    end_row = len(data) + row_offset if len(data) > 0 else row_offset + 1
    table_area = f"A{row_offset}:{column_end}{end_row}"

    # Create the table with default styling. This automatically applies an autofilter.
    worksheet.add_table(
        table_area,
        {"columns": excel_headers(headers), "data": data, "style": "Table Style Light 1"},
    )

    apply_compact_formatting(worksheet, headers, data)

    add_maxdepth_rescue_row_color(worksheet, workbook, headers, data, row_offset)

    add_manta_af_row_colors(worksheet, workbook, headers, data, row_offset)

    # Hide rows with high manta_N_AF, except target-gene variants
    if svdb_col_idx is not None and data:
        for i, row_data in enumerate(data):
            excel_row_index = row_offset + i

            if has_high_normal_af(row_data, svdb_col_idx) and not is_in_target_panel(row_data, columns):
                worksheet.set_row(excel_row_index, options={"hidden": True})

    return worksheet


def excel_headers(headers):
    """Create display-only Excel headers."""
    return [
        {
            **header,
            "header": header.get("header", "").replace("manta_", "", 1),
        }
        for header in headers
    ]


def known_fusion_events(table):
    """Return all available rows for BND events containing a KNOWN_FUSION annotation."""
    columns = column_indexes(table)
    details_idx = columns.get("details")

    if details_idx is None:
        return []

    events = {}
    for row in table.get("data", []):
        event_id = bnd_event_key(row, columns)
        events.setdefault(event_id, []).append(row)

    selected_rows = []

    for event_rows in events.values():
        has_known_fusion = any(
            "KNOWN_FUSION" in {
                part.strip()
                for part in str(row[details_idx]).split(",")
            }
            for row in event_rows
        )

        if not has_known_fusion:
            continue

        # Include both breakends of the event.
        selected_rows.extend(event_rows)

    return selected_rows


def has_high_normal_af(row, normal_af_idx):
    return (
        normal_af_idx is not None
        and as_float(row[normal_af_idx], default=0.0) > MAX_OVERVIEW_NORMAL_AF
    )


def is_in_target_panel(row, column_indexes):
    """
    True when the row is flagged as belonging to the target gene panel.
    """
    target_idx = column_indexes.get("in target panel")

    return (
        target_idx is not None
        and str(row[target_idx]).strip().lower() == "yes"
    )


def align_overview_table(table_data, selected_data):
    """Put Overview fields in the fixed COLUMN_WIDTHS order."""
    source_headers = table_data["headers"]
    source_indexes = {
        header.get("header"): index
        for index, header in enumerate(source_headers)
    }

    headers = []
    for name in COLUMN_WIDTHS:
        if name in source_indexes:
            headers.append(
                dict(source_headers[source_indexes[name]])
            )
        else:
            headers.append({"header": name})

    aligned_data = [
        [
            row[source_indexes[name]]
            if name in source_indexes
            else ""
            for name in COLUMN_WIDTHS
        ]
        for row in selected_data
    ]

    return headers, aligned_data


def _filter_overview_rows_by_normal_af(table_data, selected_data, event_atomic=False):
    """
    Remove Overview candidates with normal-panel AF above the configured limit.

    Target-panel rows are exempt, so gene-list calls that stay visible on the
    data sheets are kept here as well. They are still colour-coded by
    manta_N_AF, which flags them as frequent in the normal panel.

    BND events are removed as complete events when event_atomic is True, and
    an event is exempt when any of its breakends is in the target panel.
    """
    columns = column_indexes(table_data)
    normal_af_idx = columns.get("manta_n_af")

    if normal_af_idx is None:
        return list(selected_data)

    if not event_atomic:
        return [
            row
            for row in selected_data
            if not has_high_normal_af(row, normal_af_idx) or is_in_target_panel(row, columns)
        ]

    exempt_event_ids = {
        bnd_event_key(row, columns)
        for row in selected_data
        if is_in_target_panel(row, columns)
    }

    excluded_event_ids = {
        bnd_event_key(row, columns)
        for row in selected_data
        if has_high_normal_af(row, normal_af_idx)
    } - exempt_event_ids

    return [
        row
        for row in selected_data
        if bnd_event_key(row, columns) not in excluded_event_ids
    ]


def write_overview_summary(
    worksheet,
    workbook,
    start_row,
    title,
    table_data,
    selected_data,
    format_heading,
    empty_message=None,
    event_atomic=False,
):
    """
    Apply the common Overview filters and write a summary table.

    BND summaries are filtered as complete events when event_atomic is True.
    When no rows remain, the section is omitted unless empty_message is given.
    """
    if not table_data or "headers" not in table_data:
        return start_row

    selected_data = _filter_overview_rows_by_normal_af(
        table_data,
        selected_data,
        event_atomic=event_atomic,
    )

    if not selected_data and empty_message is None:
        return start_row

    headers, selected_data = align_overview_table(
        table_data,
        selected_data,
    )
    worksheet.write(start_row, 0, title, format_heading)
    table_start_row = start_row + 2

    if not selected_data:
        worksheet.write(
            table_start_row,
            0,
            empty_message,
        )
        return table_start_row + 3

    column_end = convert_columns_to_letter(len(headers))
    excel_start_row = table_start_row + 1
    excel_end_row = excel_start_row + len(selected_data)
    table_area = f"A{excel_start_row}:{column_end}{excel_end_row}"

    worksheet.add_table(
        table_area,
        {
            "columns": excel_headers(headers),
            "data": selected_data,
            "style": "Table Style Light 1",
        },
    )

    apply_compact_formatting(worksheet, headers, selected_data)
    add_maxdepth_rescue_row_color(worksheet, workbook, headers, selected_data, excel_start_row)
    add_manta_af_row_colors(worksheet, workbook, headers, selected_data, excel_start_row)

    return table_start_row + len(selected_data) + 4


def write_target_summary(worksheet, workbook, start_row, title, table_data, format_heading, event_atomic=False):
    """
    Select target-panel variants and write them to the Overview sheet.
    When event_atomic is True, all rows belonging to a selected BND event
    are included if at least one breakend is in the target panel.
    """
    if not table_data or "headers" not in table_data:
        return start_row

    data = table_data.get("data", [])

    columns = column_indexes(table_data)

    if "in target panel" not in columns:
        return start_row

    target_data = [
        row
        for row in data
        if is_in_target_panel(row, columns)
    ]

    if event_atomic and target_data:

        selected_event_ids = {
            bnd_event_key(row, columns)
            for row in target_data
        }

        target_data = [
            row
            for row in data
            if bnd_event_key(row, columns) in selected_event_ids
        ]

    return write_overview_summary(
        worksheet,
        workbook,
        start_row,
        title,
        table_data,
        target_data,
        format_heading,
        empty_message="Inga target-varianter hittades för denna typ.",
        event_atomic=event_atomic,
    )





""" MAIN EXECUTION """

# 1. Prepping data
logging.info(f"Prepping data, such as loading {snakemake.input.manta}=")
sample_name = snakemake.output.xlsx.split("/")[-1].split(".manta_new.xlsx")[0]

filter_flags = ["MinQUAL", "MinGQ", "MinSomaticScore", "Ploidy", "MaxMQ0Frac", "NoPairSupport", "SampleFT", "HomRef", "MaxDepth"]

# Load target genes for easy filtering in Excel
target_genes_path = getattr(snakemake.input, "target_genes", "") or getattr(snakemake.params, "target_genes", "")
target_genes = load_target_genes(target_genes_path)

manta_tables_full = create_manta_tables(snakemake.input.manta, filter_flags, target_genes=target_genes)
manta_tables_all = create_manta_tables(snakemake.input.manta, avoid_filterflags=[], target_genes=target_genes)
blocking_filter_flags = set(filter_flags) - {"MaxDepth"}
manta_tables_maxdepth = create_maxdepth_bnd_rescue_table(
    manta_tables_all,
    blocking_filter_flags,
    min_support=MAXDEPTH_RESCUE_MIN_SUPPORT,
)

# 2. Creating xlsx workbook
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
logging.info(f"Creating xlsx workbook {snakemake.output.xlsx}=")

format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_bold = workbook.add_format({"bold": True, "text_wrap": True})
format_overview_title = workbook.add_format({"bold": True, "font_size": 16})
format_2dec = workbook.add_format({"num_format": "0.00"})

manta_tables_full["bnd"] = merge_manta_tables(manta_tables_full["bnd"], manta_tables_maxdepth)

panel_tables_dict = {}

for vcf in snakemake.input.vcfs_bed:
    panel = vcf.split(".")[-3]

    panel_tables = create_manta_tables(vcf, filter_flags, target_genes=target_genes)
    panel_tables_all = create_manta_tables(vcf, avoid_filterflags=[], target_genes=target_genes)
    panel_rescue = select_rescue_events_for_panel(manta_tables_maxdepth, panel_tables_all["bnd"])

    panel_tables["bnd"] = merge_manta_tables(panel_tables["bnd"], panel_rescue)
    add_maxdepth_rescue_column(panel_tables["bnd"], panel_rescue)
    format_manta_table(panel_tables["bnd"], sample_name, format_2dec)

    panel_tables_dict[panel] = panel_tables

add_maxdepth_rescue_column(manta_tables_full["bnd"], manta_tables_maxdepth)
add_maxdepth_rescue_column(manta_tables_maxdepth, manta_tables_maxdepth)

for sv_key in ["del", "ins", "dup", "bnd"]:
    format_manta_table(manta_tables_full[sv_key], sample_name, format_2dec)

format_manta_table(manta_tables_maxdepth, sample_name, format_2dec)

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
    panel_tables = panel_tables_dict[panel]
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
worksheet_overview.write(row_idx + 10, 0, "MaxDepth calls passing rescue criteria are included.")

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
        workbook,
        row_idx,
        sv_title,
        manta_tables_full[sv_key],
        format_bold,
        event_atomic=(sv_key == "bnd"),
    )

known_fusion_rows = known_fusion_events(manta_tables_full["bnd"])

known_fusion_start = row_idx + 2

new_row_idx = write_overview_summary(
    worksheet_overview,
    workbook,
    known_fusion_start,
    "Translocations (KNOWN_FUSION)",
    manta_tables_full["bnd"],
    known_fusion_rows,
    format_overview_title,
    event_atomic=True,
)

if new_row_idx != known_fusion_start:
    row_idx = new_row_idx

maxdepth_start = row_idx + 2

new_row_idx = write_overview_summary(
    worksheet_overview,
    workbook,
    maxdepth_start,
    "MaxDepth rescue calls",
    manta_tables_maxdepth,
    manta_tables_maxdepth.get("data", []),
    format_overview_title,
    event_atomic=True,
)

if new_row_idx != maxdepth_start:
    row_idx = new_row_idx

workbook.set_size(1800, 1200)
workbook.close()
logging.info("All done")
