#!/usr/bin/env python

from export_to_xlsx_create_tables import *
import xlsxwriter
import datetime
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


def write_target_summary(worksheet, start_row, title, table_data, format_heading, event_atomic=False):
    """
    Write variants annotated as belonging to the target panel.

    When event_atomic is True, all rows belonging to a selected BND event
    are included if at least one breakend is in the target panel.
    """
    if not table_data or "headers" not in table_data:
        return start_row

    headers = table_data["headers"]
    data = table_data.get("data", [])

    target_col_idx = next(
        (
            idx
            for idx, header in enumerate(headers)
            if header.get("header") == "In Target Panel"
        ),
        None,
    )

    if target_col_idx is None:
        return start_row

    target_data = [
        row
        for row in data
        if str(row[target_col_idx]).strip().lower() == "yes"
    ]

    if event_atomic and target_data:
        columns = _column_indexes(table_data)

        selected_event_ids = {
            _bnd_event_key(row, columns)
            for row in target_data
        }

        target_data = [
            row
            for row in data
            if _bnd_event_key(row, columns) in selected_event_ids
        ]

    worksheet.write(start_row, 0, title, format_heading)
    table_start_idx = start_row + 2

    if not target_data:
        worksheet.write(
            table_start_idx,
            0,
            "Inga target-varianter hittades för denna typ.",
        )
        return table_start_idx + 3

    column_end = convert_columns_to_letter(len(headers))
    excel_start_row = table_start_idx + 1
    excel_end_row = excel_start_row + len(target_data)
    table_area = f"A{excel_start_row}:{column_end}{excel_end_row}"

    worksheet.add_table(
        table_area,
        {
            "columns": headers,
            "data": target_data,
            "style": "Table Style Light 1",
        },
    )

    apply_compact_formatting(
        worksheet,
        headers,
        target_data,
    )

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
    row_offset = 7
    if "Deletions" in sheet_name:
        worksheet.write("A7", "Calls have to be longer than 100 bp to be included.")
        row_offset = 8

    headers = table_data["headers"]
    data = table_data["data"]

    # 1. Find columns
    svdb_col_idx = -1
    target_col_idx = -1
    for idx, header_dict in enumerate(headers):
        header_name = header_dict.get("header")

        if header_name == "manta_N_AF":
            svdb_col_idx = idx
        elif header_name == "In Target Panel":
            target_col_idx = idx

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

    # Hide rows with high manta_N_AF, except target-gene variants
    if svdb_col_idx != -1 and data:
        for i, row_data in enumerate(data):
            excel_row_index = row_offset + i

            is_target = (
                target_col_idx != -1
                and str(row_data[target_col_idx]).strip().lower() == "yes"
            )

            try:
                high_normal_af = float(row_data[svdb_col_idx]) > 0.2
            except (ValueError, TypeError):
                high_normal_af = False

            if high_normal_af and not is_target:
                worksheet.set_row(excel_row_index, options={"hidden": True})

    return worksheet


def _column_indexes(table):
    return {
        str(header.get("header", "")).lower(): index
        for index, header in enumerate(table.get("headers", []))
    }


def _as_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _annotation_flags(row, column_indexes):
    return {
        flag.strip()
        for flag in str(row[column_indexes["annotation"]]).split(",")
        if flag.strip()
    }


def _bnd_event_key(row, column_indexes):
    """Return a stable BND event key, falling back to a single-record event."""
    event_idx = column_indexes.get("bnd event id")
    if event_idx is not None and row[event_idx]:
        return str(row[event_idx])

    id_idx = column_indexes.get("mantaid")
    return str(row[id_idx]) if id_idx is not None else id(row)


def _is_junk_bnd(row, column_indexes):
    chr_idx = column_indexes.get("chr")
    partner_idx = column_indexes.get("breakend", column_indexes.get("alt"))

    chr_value = str(row[chr_idx]).lower() if chr_idx is not None else ""
    partner_value = str(row[partner_idx]).lower() if partner_idx is not None else ""

    return any(
        keyword in value
        for keyword in ("chrun", "random", "gl0", "ki2")
        for value in (chr_value, partner_value)
    )


def known_fusion_events(table):
    """Return complete BND events containing a KNOWN_FUSION annotation."""
    columns = _column_indexes(table)
    details_idx = columns.get("details")

    if details_idx is None:
        return []

    events = {}
    for row in table.get("data", []):
        event_id = _bnd_event_key(row, columns)
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


def write_known_fusion_summary(worksheet, start_row, title, headers, selected_data, format_heading):
    """Write selected KNOWN_FUSION BND events to the Overview sheet."""
    if not selected_data:
        return start_row

    worksheet.write(start_row, 0, title, format_heading)
    table_start_row = start_row + 2

    column_end = convert_columns_to_letter(len(headers))
    excel_start_row = table_start_row + 1
    excel_end_row = excel_start_row + len(selected_data)
    table_area = f"A{excel_start_row}:{column_end}{excel_end_row}"

    worksheet.add_table(
        table_area,
        {
            "columns": headers,
            "data": selected_data,
            "style": "Table Style Light 1",
        },
    )

    apply_compact_formatting(worksheet, headers, selected_data)

    return table_start_row + len(selected_data) + 4


def create_maxdepth_bnd_rescue_table(tables_dict, blocking_filter_flags, min_support=0.05):
    """
    Return BND events classified as MaxDepth if they pass our rescue criteria.

    The input tables should contain all BND records so that MATEID-linked
    partners are available. An event is rescued only when:

    - at least one breakpoint is flagged MaxDepth;
    - no event row has another blocking filter;
    - every MaxDepth breakpoint has PR or SR frequency >= min_support;
    - no MaxDepth breakpoint is present in the normal panel;
    - neither breakpoint points to a junk/alternate contig.

    When accepted, every row belonging to the BND event is returned.
    """
    if not 0 <= min_support <= 1:
        raise ValueError("min_support must be between 0 and 1")

    bnd_table = tables_dict.get("bnd", {"headers": [], "data": []})
    rescue_table = {
        "headers": [
            header.copy()
            for header in bnd_table.get("headers", [])
        ],
        "data": [],
    }

    if not bnd_table.get("data"):
        return rescue_table

    columns = _column_indexes(bnd_table)

    required_columns = {
        "bnd event id",
        "chr",
        "annotation",
        "paired-read freq",
        "spanning-read freq",
        "manta_n_occ",
    }

    missing_columns = required_columns - columns.keys()

    if "breakend" not in columns and "alt" not in columns:
        missing_columns.add("breakend/alt")

    if missing_columns:
        raise ValueError(
            "Cannot apply MaxDepth rescue to BND table; "
            f"missing columns: {', '.join(sorted(missing_columns))}"
        )

    events = {}
    for row in bnd_table["data"]:
        event_id = _bnd_event_key(row, columns)
        events.setdefault(event_id, []).append(row)

    rescue_stats = Counter()

    for event_rows in events.values():
        maxdepth_rows = [
            row
            for row in event_rows
            if "MaxDepth" in _annotation_flags(row, columns)
        ]

        if not maxdepth_rows:
            continue

        rescue_stats["candidate_events"] += 1

        has_other_blocking_filter = any(
            _annotation_flags(row, columns) & blocking_filter_flags
            for row in event_rows
        )

        all_maxdepth_rows_have_support = all(
            _as_float(row[columns["paired-read freq"]]) >= min_support
            or _as_float(row[columns["spanning-read freq"]]) >= min_support
            for row in maxdepth_rows
        )

        maxdepth_row_has_normal_panel_hit = any(
            _as_float(row[columns["manta_n_occ"]]) != 0
            for row in maxdepth_rows
        )

        has_junk_contig = any(
            _is_junk_bnd(row, columns)
            for row in event_rows
        )

        if has_other_blocking_filter:
            rescue_stats["other_filter"] += 1
        elif not all_maxdepth_rows_have_support:
            rescue_stats["low_support"] += 1
        elif maxdepth_row_has_normal_panel_hit:
            rescue_stats["normal_panel"] += 1
        elif has_junk_contig:
            rescue_stats["junk_contig"] += 1
        else:
            rescue_table["data"].extend(event_rows)
            rescue_stats["rescued_events"] += 1
            rescue_stats["rescued_records"] += len(event_rows)

    if rescue_stats:
        logging.info(
            "MaxDepth rescue for BND: %s",
            dict(rescue_stats),
        )

    return rescue_table


def write_maxdepth_rescue_summary(worksheet, start_row, title, table_data, format_heading):
    """Write accepted MaxDepth BND rescue calls to the Overview sheet."""
    if not table_data or not table_data.get("data"):
        return start_row

    headers = table_data["headers"]
    data = table_data["data"]

    worksheet.write(start_row, 0, title, format_heading)
    table_start_row = start_row + 2

    column_end = convert_columns_to_letter(len(headers))
    excel_start_row = table_start_row + 1
    excel_end_row = excel_start_row + len(data)
    table_area = f"A{excel_start_row}:{column_end}{excel_end_row}"

    worksheet.add_table(
        table_area,
        {"columns": headers, "data": data, "style": "Table Style Light 1"},
    )
    apply_compact_formatting(worksheet, headers, data)

    return table_start_row + len(data) + 4


""" MAIN EXECUTION """

# 1. Prepping data
logging.info(f"Prepping data, such as loading {snakemake.input.manta}=")
sample_name = snakemake.output.xlsx.split("/")[-1].split(".manta.xlsx")[0]

filter_flags = ["MinQUAL", "MinGQ", "MinSomaticScore", "Ploidy", "MaxMQ0Frac", "NoPairSupport", "SampleFT", "HomRef", "MaxDepth"]

# Load target genes for easy filtering in Excel
target_genes = []
target_genes_path = getattr(snakemake.params, "target_genes", "") or getattr(snakemake.input, "target_genes", "")
if target_genes_path:
    target_genes = load_target_genes(target_genes_path)

manta_tables_full = create_manta_tables(snakemake.input.manta, filter_flags, target_genes=target_genes)
manta_tables_all = create_manta_tables(snakemake.input.manta, avoid_filterflags=[], target_genes=target_genes)
blocking_filter_flags = set(filter_flags) - {"MaxDepth"}
manta_tables_maxdepth = create_maxdepth_bnd_rescue_table(
    manta_tables_all,
    blocking_filter_flags,
    min_support=0.05,
)

# 2. Creating xlsx workbook
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
logging.info(f"Creating xlsx workbook {snakemake.output.xlsx}=")

format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_bold = workbook.add_format({"bold": True, "text_wrap": True})
format_overview_title = workbook.add_format({"bold": True, "font_size": 16})
format_2dec = workbook.add_format({"num_format": "0.00"})

for sv_key in ["del", "ins", "dup", "bnd"]:
    format_manta_table(manta_tables_full[sv_key], sample_name, format_2dec)

format_manta_table(manta_tables_maxdepth, sample_name, format_2dec)

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
        format_bold,
        event_atomic=(sv_key == "bnd"),
    )

known_fusion_rows = known_fusion_events(manta_tables_full["bnd"])

if known_fusion_rows:
    row_idx += 2
    row_idx = write_known_fusion_summary(
        worksheet_overview,
        row_idx,
        "Translocations (KNOWN_FUSION)",
        manta_tables_full["bnd"]["headers"],
        known_fusion_rows,
        format_overview_title,
    )

if manta_tables_maxdepth.get("data"):
    row_idx += 2
    row_idx = write_maxdepth_rescue_summary(
        worksheet_overview,
        row_idx,
        "MaxDepth rescue calls",
        manta_tables_maxdepth,
        format_overview_title,
    )

workbook.set_size(1800, 1200)
workbook.close()
logging.info("All done")
