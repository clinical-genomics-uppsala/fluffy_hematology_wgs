import logging
from collections import Counter

MAXDEPTH_RESCUE_MIN_SUPPORT = 0.05


def column_indexes(table):
    return {
        str(header.get("header", "")).lower(): index
        for index, header in enumerate(table.get("headers", []))
    }


def as_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def bnd_event_key(row, column_indexes):
    """Return a stable BND event key, falling back to a single-record event."""
    event_idx = column_indexes.get("bnd event id")
    if event_idx is not None and row[event_idx]:
        return str(row[event_idx])

    id_idx = column_indexes.get("mantaid")
    return str(row[id_idx]) if id_idx is not None else id(row)


def _annotation_flags(row, column_indexes):
    return {
        flag.strip()
        for flag in str(row[column_indexes["annotation"]]).split(",")
        if flag.strip()
    }


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


def add_maxdepth_rescue_column(table, rescue_table):
    """Mark both breakends belonging to a rescued MaxDepth event."""
    if not table or "headers" not in table:
        return

    columns = column_indexes(table)
    rescue_columns = column_indexes(rescue_table)

    rescue_event_ids = {
        bnd_event_key(row, rescue_columns)
        for row in rescue_table.get("data", [])
    }

    rescue_header_idx = columns.get("maxdepth rescue")

    if rescue_header_idx is None:
        table["headers"].append({"header": "MaxDepth Rescue"})
        rescue_header_idx = len(table["headers"]) - 1

        for row in table["data"]:
            row.append("")

    columns = column_indexes(table)

    for row in table["data"]:
        event_id = bnd_event_key(row, columns)
        row[rescue_header_idx] = (
            "Yes" if event_id in rescue_event_ids else ""
        )


def create_maxdepth_bnd_rescue_table(tables_dict, blocking_filter_flags, min_support=MAXDEPTH_RESCUE_MIN_SUPPORT):
    """
    Return BND events classified as MaxDepth if they pass our rescue criteria.

    Manta flags MaxDepth on breakpoints in regions with depth 3x above the
    chromosome mean. Such regions are enriched for mapping artefacts, but strictly removing events based on
    maxdepth was found to remove some true variants. The calls are still kept off the ordinary sheets and only
    the subset below is surfaced, in its own Overview section, for manual review.

    An event is rescued only when:

    - at least one breakpoint is flagged MaxDepth;
    - no breakpoint has another blocking filter;
    - every MaxDepth breakpoint has PR or SR frequency >= min_support (0.05);
    - neither breakpoint is present in the normal pane or in a junk/alternate contig.

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

    columns = column_indexes(bnd_table)

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
        event_id = bnd_event_key(row, columns)
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
            as_float(row[columns["paired-read freq"]]) >= min_support
            or as_float(row[columns["spanning-read freq"]]) >= min_support
            for row in maxdepth_rows
        )

        event_has_normal_panel_hit = any(
            as_float(row[columns["manta_n_occ"]]) != 0
            for row in event_rows
        )

        has_junk_contig = any(
            _is_junk_bnd(row, columns)
            for row in event_rows
        )

        if has_other_blocking_filter:
            rescue_stats["other_filter"] += 1
        elif not all_maxdepth_rows_have_support:
            rescue_stats["low_support"] += 1
        elif event_has_normal_panel_hit:
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


def merge_manta_tables(base_table, extra_table):
    if base_table["headers"] != extra_table["headers"]:
        raise ValueError("Cannot merge Manta tables with different headers")

    columns = column_indexes(base_table)
    manta_id_idx = columns.get("mantaid")

    if manta_id_idx is None:
        raise ValueError("Cannot merge Manta tables without MantaID")

    existing_ids = {
        str(row[manta_id_idx])
        for row in base_table.get("data", [])
    }

    for row in extra_table.get("data", []):
        manta_id = str(row[manta_id_idx])
        if manta_id not in existing_ids:
            base_table["data"].append(row.copy())
            existing_ids.add(manta_id)

    return base_table


def select_rescue_events_for_panel(rescue_table, panel_table):
    """Select complete rescue events when either breakend occurs in the panel VCF."""
    rescue_columns = column_indexes(rescue_table)
    panel_columns = column_indexes(panel_table)

    panel_event_ids = {
        bnd_event_key(row, panel_columns)
        for row in panel_table.get("data", [])
    }

    return {
        "headers": [header.copy() for header in rescue_table.get("headers", [])],
        "data": [
            row.copy()
            for row in rescue_table.get("data", [])
            if bnd_event_key(row, rescue_columns) in panel_event_ids
        ],
    }


__all__ = [
    "MAXDEPTH_RESCUE_MIN_SUPPORT",
    "add_maxdepth_rescue_column",
    "as_float",
    "bnd_event_key",
    "column_indexes",
    "create_maxdepth_bnd_rescue_table",
    "merge_manta_tables",
    "select_rescue_events_for_panel",
]
