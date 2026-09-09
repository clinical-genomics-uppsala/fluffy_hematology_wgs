#!/usr/bin/env python

import gzip
from pysam import VariantFile
import re

MIN_DEL_SIZE = 100


# VEP fields in list to get index
def index_vep(variantfile):
    csq_index = []
    for x in variantfile.header.records:
        if "CSQ" in str(x):
            csq_index = str(x).split("Format: ")[1].strip().strip('">').split("|")
    return csq_index


# Extract table columns from vcf records
def extract_vcf_values(record, csq_index, sample_tumor, sample_normal=""):
    return_dict = {}
    return_dict["filter_flag"] = ",".join(record.filter.keys())

    # --- PINDEL AF FIX ---
    return_dict["af"] = 0.0
    try:
        return_dict["af"] = float(record.samples[sample_tumor]["AF"][0])
    except (KeyError, TypeError, IndexError):
        ad = record.samples[sample_tumor].get("AD")
        if ad and len(ad) >= 2 and sum(ad) > 0:
            return_dict["af"] = float(ad[1]) / float(sum(ad))

    if sample_normal != "":
        try:
            return_dict["n_af"] = float(record.samples[sample_normal]["AF"][0])
        except (KeyError, TypeError, IndexError):
            ad_n = record.samples[sample_normal].get("AD")
            if ad_n and len(ad_n) >= 2 and sum(ad_n) > 0:
                return_dict["n_af"] = float(ad_n[1]) / float(sum(ad_n))
            else:
                return_dict["n_af"] = ""
    else:
        return_dict["n_af"] = ""

    try:
        return_dict["dp"] = int(record.samples[sample_tumor]["DP"])
    except (KeyError, TypeError):
        ad = record.samples[sample_tumor].get("AD")
        if ad:
            return_dict["dp"] = sum(ad)
        else:
            return_dict["dp"] = 0

    if sample_normal != "":
        try:
            return_dict["n_dp"] = int(record.samples[sample_normal]["DP"])
        except (KeyError, TypeError):
            ad_n = record.samples[sample_normal].get("AD")
            if ad_n:
                return_dict["n_dp"] = sum(ad_n)
            else:
                return_dict["n_dp"] = ""
    else:
        return_dict["n_dp"] = ""

    try:
        return_dict["svlen"] = int(record.info["SVLEN"])
    except KeyError:
        return_dict["svlen"] = ""

    try:
        csq = record.info["CSQ"][0].split("|")
    except KeyError:
        csq = None

    # vep annotation
    if csq:
        return_dict["gene"] = csq[csq_index.index("SYMBOL")]
        return_dict["transcript"] = csq[csq_index.index("HGVSc")].split(":")[0]

        try:
            return_dict["exon_nr"] = csq[csq_index.index("EXON")]
        except KeyError:
            return_dict["exon_nr"] = ""

        if len(csq[csq_index.index("HGVSc")].split(":")) > 1:
            return_dict["coding_name"] = csq[csq_index.index("HGVSc")].split(":")[1]
        else:
            return_dict["coding_name"] = ""
        return_dict["ensp"] = csq[csq_index.index("HGVSp")]
        return_dict["consequence"] = csq[csq_index.index("Consequence")]

        existing = csq[csq_index.index("Existing_variation")].split("&")
        cosmic_list = [cosmic for cosmic in existing if cosmic.startswith("CO")]
        if len(cosmic_list) == 0:
            return_dict["cosmic"] = ""
        else:
            return_dict["cosmic"] = ", ".join(cosmic_list)

        return_dict["clinical"] = csq[csq_index.index("CLIN_SIG")]

        rs_list = [rs for rs in existing if rs.startswith("rs")]
        if len(rs_list) == 0:
            return_dict["rs"] = ""
        else:
            return_dict["rs"] = ", ".join(rs_list)
        return_dict["max_pop_af"] = csq[csq_index.index("MAX_AF")]
        return_dict["max_pops"] = csq[csq_index.index("MAX_AF_POPS")]
    else:
        # --- UPDATE FIX FÖR ATT BEHÅLLA AF ---
        return_dict.update(dict.fromkeys(
            [
                "gene",
                "transcript",
                "exon_nr",
                "coding_name",
                "ensp",
                "consequence",
                "cosmic",
                "clinical",
                "rs",
                "max_pop_af",
                "max_pops",
            ],
            "",
        ))

    return return_dict


#  Helper functions for value extraction
def first_info_value(record, key, default=None):
    """Return the first scalar INFO value, or default when absent."""
    value = record.info.get(key)

    if isinstance(value, (tuple, list)):
        value = value[0] if value else None

    if value in (None, "", "."):
        return default

    return value


def info_as_int(record, key, default=0):
    value = first_info_value(record, key)
    if value is None:
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def info_as_float(record, key, default=0.0, decimals=None):
    value = first_info_value(record, key)
    if value is None:
        return default
    try:
        value = float(value)
    except (TypeError, ValueError):
        return default

    return round(value, decimals) if decimals is not None else value


def support_frequency(sample, field, default=None):
    try:
        values = sample.get(field)
        if not values or len(values) != 2:
            return default
        ref_count, alt_count = values
        if ref_count is None or alt_count is None:
            return default
        total = ref_count + alt_count
        return alt_count / total if total > 0 else default
    except (TypeError, ValueError):
        return default


def index_manta_annotation_fields(variantfile, annotation_id):
    """Read the pipe-delimited field definition from a Manta INFO header."""
    for header_record in variantfile.header.records:
        if header_record.key != "INFO" or header_record.get("ID") != annotation_id:
            continue
        description = header_record.get("Description", "")
        match = re.search(r"'([^']+)'", description)
        if match:
            return [field.strip() for field in match.group(1).split("|")]
        raise ValueError(f"{annotation_id} header has no quoted field definition")
    return None


def _info_values(record, key):
    """Return all INFO values as a tuple."""
    value = record.info.get(key)

    if value in (None, "", "."):
        return ()

    if isinstance(value, (tuple, list)):
        return tuple(value)

    return (value,)


def _extract_manta_annotations(record, ann_index, simple_ann_index):
    """Extract compact gene labels and details, preferring SIMPLE_ANN."""
    simple_annotations = _info_values(record, "SIMPLE_ANN")

    if simple_annotations:
        if not simple_ann_index:
            raise ValueError(
                "Manta VCF contains SIMPLE_ANN records "
                "but no SIMPLE_ANN header definition"
            )

        gene_idx = simple_ann_index.index("GENE(s)")
        transcript_idx = simple_ann_index.index("TRANSCRIPT")
        detail_idx = simple_ann_index.index(
            "DETAIL (exon losses, KNOWN_FUSION, "
            "ON_PRIORITY_LIST, NOT_PRIORITISED)"
        )

        genes = []
        details = []

        for annotation in simple_annotations:
            values = annotation.split("|")

            gene_label = (
                f"{values[gene_idx]}({values[transcript_idx]})"
            )
            if gene_label not in genes:
                genes.append(gene_label)

            detail = values[detail_idx]
            if detail not in details:
                details.append(detail)

        return ", ".join(genes), ", ".join(details)

    annotations = _info_values(record, "ANN")

    if annotations:
        if not ann_index:
            raise ValueError(
                "Manta VCF contains ANN records "
                "but no ANN header definition"
            )

        gene_name_idx = ann_index.index("Gene_Name")
        gene_id_idx = ann_index.index("Gene_ID")

        genes = []

        for annotation in annotations:
            values = annotation.split("|")
            gene_label = (
                f"{values[gene_name_idx]}({values[gene_id_idx]})"
            )
            if gene_label not in genes:
                genes.append(gene_label)

        return ", ".join(genes), ""

    return "NA", "NA"


def extract_manta_vcf_values(record, ann_index, simple_ann_index, sample_tumor, sample_normal=""):
    return_dict = {}
    return_dict["filt_ann"] = ",".join(record.filter.keys())
    return_dict["id"] = record.id or ""
    return_dict["depth"] = first_info_value(record, "BND_DEPTH", default="")
    # helper functions sets svdb vales to 0 if missing:
    return_dict["manta_n_occ"] = info_as_int(record, "manta_N_OCC")
    return_dict["manta_t_occ"] = info_as_int(record, "manta_T_OCC")
    return_dict["manta_n_af"] = info_as_float(record, "manta_N_AF")
    return_dict["manta_t_af"] = info_as_float(record, "manta_T_AF")
    return_dict["str_percent"] = info_as_float(record, "STR_PERCENT", default=0, decimals=2)
    return_dict["genes"], return_dict["detail"] = _extract_manta_annotations(record, ann_index, simple_ann_index)

    # create a common id for BND mates
    mate_id = first_info_value(record, "MATEID")
    if mate_id is not None:
        return_dict["bnd_event_id"] = "|".join(sorted([return_dict["id"], str(mate_id)]))
    else:
        return_dict["bnd_event_id"] = return_dict["id"]

    # extract paired read and spanning read frequncies
    tumor_sample = record.samples[sample_tumor]
    return_dict["pr_freq"] = support_frequency(tumor_sample, "PR")
    return_dict["sr_freq"] = support_frequency(tumor_sample, "SR")

    if sample_normal:
        normal_sample = record.samples[sample_normal]
        return_dict["pr_freq_n"] = support_frequency(normal_sample, "PR")
        return_dict["sr_freq_n"] = support_frequency(normal_sample, "SR")

    return_dict["svlength"] = info_as_int(record, "SVLEN", default=None)
    return_dict["hom_len"] = info_as_int(record, "HOMLEN", default=None)
    return_dict["hom_seq"] = first_info_value(record, "HOMSEQ", default=None)

    return return_dict


def create_snv_table(vcf_input):
    vcf_file = VariantFile(vcf_input)
    sample_tumor = [x for x in list(vcf_file.header.samples) if x.endswith("_T")][0]
    if len(list(vcf_file.header.samples)) > 1:
        sample_normal = [x for x in list(vcf_file.header.samples) if x.endswith("_N")][0]
    else:
        sample_normal = ""
    csq_index = index_vep(vcf_file)
    vep_line = [x.value for x in vcf_file.header.records if x.key == "VEP"][0]

    snv_table = {"data": [], "headers": [], "vep_line": vep_line}
    snv_table["headers"] = [
        {"header": "FilterFlag"},
        {"header": "DNAnr"},
        {"header": "Gene"},
        {"header": "Chr"},
        {"header": "Pos"},
        {"header": "Ref"},
        {"header": "Alt"},
        {"header": "AF"},
        {"header": "Normal AF"},
        {"header": "DP"},
        {"header": "Normal DP"},
        {"header": "Transcript"},
        {"header": "Exon"},
        {"header": "Mutation cds"},
        {"header": "ENSP"},
        {"header": "Consequence"},
        {"header": "COSMIC ids on pos"},
        {"header": "Clinical Significance"},
        {"header": "dbSNP"},
        {"header": "Max Pop AF"},
        {"header": "Max Pop"},
    ]
    for record in vcf_file.fetch():
        record_values = extract_vcf_values(record, csq_index, sample_tumor, sample_normal)
        if record_values["af"] > 0.01:
            outline = [
                record_values["filter_flag"],
                sample_tumor,
                record_values["gene"],
                record.contig,
                int(record.pos),
                record.ref,
                record.alts[0],
                record_values["af"],
                record_values["n_af"],
                record_values["dp"],
                record_values["n_dp"],
                record_values["transcript"],
                record_values["exon_nr"],
                record_values["coding_name"],
                record_values["ensp"],
                record_values["consequence"],
                record_values["cosmic"],
                record_values["clinical"],
                record_values["rs"],
                record_values["max_pop_af"],
                record_values["max_pops"],
            ]
            snv_table["data"].append(outline)
    return snv_table


def create_pindel_table(vcf_input):
    pindel_file = VariantFile(vcf_input)
    sample = list(pindel_file.header.samples)[0]
    csq_index = index_vep(pindel_file)
    vep_line = [x.value for x in pindel_file.header.records if x.key == "VEP"][0]

    pindel_table = {"data": [], "headers": [], "vep_line": vep_line}
    pindel_table["headers"] = [
        {"header": "Filter"},
        {"header": "DNAnr"},
        {"header": "Gene"},
        {"header": "Chr"},
        {"header": "Pos"},
        {"header": "Ref"},
        {"header": "Alt"},
        {"header": "SV length"},
        {"header": "AF"},
        {"header": "DP"},
        {"header": "Transcript"},
        {"header": "Mutation cds"},
        {"header": "ENSP"},
        {"header": "Consequence"},
        {"header": "COSMIC ids on pos"},
        {"header": "Clinical Significance"},
        {"header": "dbSNP"},
        {"header": "Max Pop AF"},
        {"header": "Max Pop"},
    ]
    for record in pindel_file.fetch():
        record_values = extract_vcf_values(record, csq_index, sample)
        if record_values["af"] > 0.01:
            outline = [
                record_values["filter_flag"],
                sample,
                record_values["gene"],
                record.contig,
                int(record.pos),
                record.ref,
                record.alts[0],
                record_values["svlen"],
                record_values["af"],
                record_values["dp"],
                record_values["transcript"],
                record_values["coding_name"],
                record_values["ensp"],
                record_values["consequence"],
                record_values["cosmic"],
                record_values["clinical"],
                record_values["rs"],
                record_values["max_pop_af"],
                record_values["max_pops"],
            ]
            pindel_table["data"].append(outline)
    return pindel_table


def create_manta_tables(
    vcf_input,
    avoid_filterflags=[
        "MinQUAL",
        "MinGQ",
        "MinSomaticScore",
        "Ploidy",
        "MaxDepth",
        "MaxMQ0Frac",
        "NoPairSupport",
        "SampleFT",
        "HomRef",
    ],
    target_genes=None,
):
    vcf_file = VariantFile(vcf_input)
    sample_tumor = [x for x in list(vcf_file.header.samples) if x.endswith("_T")][0]
    if len(list(vcf_file.header.samples)) > 1:
        sample_normal = [x for x in list(vcf_file.header.samples) if x.endswith("_N")][0]
    else:
        sample_normal = None

    ann_index = index_manta_annotation_fields(vcf_file, "ANN")
    simple_ann_index = index_manta_annotation_fields(vcf_file, "SIMPLE_ANN")

    manta_tables = {
        "bnd": {"data": [], "headers": []},
        "del": {"data": [], "headers": []},
        "dup": {"data": [], "headers": []},
        "ins": {"data": [], "headers": []},
    }

    manta_tables["bnd"]["headers"] = [
        {"header": "Chr"},
        {"header": "Pos"},
        {"header": "MantaID"},
        {"header": "BND Event ID"},
        {"header": "BreakEnd"},
        {"header": "Genes"},
        {"header": "Details"},
        {"header": "Depth"},
        {"header": "Annotation"},
        {"header": "manta_N_OCC"},
        {"header": "manta_T_OCC"},
        {"header": "manta_N_AF"},
        {"header": "manta_T_AF"},
        {"header": "STR %"},
        {"header": "Paired-read freq"},
        {"header": "Spanning-read freq"},
    ]
    manta_tables["del"]["headers"] = [
        {"header": "Chr"},
        {"header": "Pos"},
        {"header": "EndPos"},
        {"header": "SV Length"},
        {"header": "MantaID"},
        {"header": "Genes"},
        {"header": "Details"},
        {"header": "Annotation"},
        {"header": "manta_N_OCC"},
        {"header": "manta_T_OCC"},
        {"header": "manta_N_AF"},
        {"header": "manta_T_AF"},
        {"header": "STR %"},
        {"header": "Paired-read freq"},
        {"header": "Spanning-read freq"},
    ]
    manta_tables["dup"]["headers"] = [
        {"header": "Chr"},
        {"header": "Pos"},
        {"header": "EndPos"},
        {"header": "SV Length"},
        {"header": "MantaID"},
        {"header": "Genes"},
        {"header": "Details"},
        {"header": "Hom Length"},
        {"header": "Hom Sequence"},
        {"header": "Annotation"},
        {"header": "manta_N_OCC"},
        {"header": "manta_T_OCC"},
        {"header": "manta_N_AF"},
        {"header": "manta_T_AF"},
        {"header": "STR %"},
        {"header": "Paired-read freq"},
        {"header": "Spanning-read freq"},
    ]
    manta_tables["ins"]["headers"] = [
        {"header": "Chr"},
        {"header": "Pos"},
        {"header": "Ref"},
        {"header": "Alt"},
        {"header": "SV Length"},
        {"header": "MantaID"},
        {"header": "Genes"},
        {"header": "Details"},
        {"header": "Hom Length"},
        {"header": "Hom Sequence"},
        {"header": "Annotation"},
        {"header": "manta_N_OCC"},
        {"header": "manta_T_OCC"},
        {"header": "manta_N_AF"},
        {"header": "manta_T_AF"},
        {"header": "STR %"},
        {"header": "Paired-read freq"},
        {"header": "Spanning-read freq"},
    ]

    if sample_normal:
        for sv_type in ["bnd", "del", "dup", "ins"]:
            manta_tables[sv_type]["headers"] += [
                {"header": "Paired-read Normal Freq"},
                {"header": "Spanning-read Normal Freq"},
            ]

    gene_pattern = None
    if target_genes:
        pattern_str = '|'.join([re.escape(g) for g in target_genes])
        gene_pattern = re.compile(rf'\b({pattern_str})\b', re.IGNORECASE)
        for sv_type in ["bnd", "del", "dup", "ins"]:
            manta_tables[sv_type]["headers"].append({"header": "In Target Panel"})

    for record in vcf_file.fetch():
        record_values = extract_manta_vcf_values(record, ann_index, simple_ann_index, sample_tumor, sample_normal)
        if not any(x in avoid_filterflags for x in record_values["filt_ann"].split(",")):
            in_target = []
            if gene_pattern:
                is_match = "Yes" if gene_pattern.search(record_values["genes"]) else "No"
                in_target = [is_match]

            if "MantaBND" in record_values["id"]:
                outline = [
                    str(record.contig),
                    int(record.pos),
                    record_values["id"],
                    record_values["bnd_event_id"],
                    str(record.alts[0]),
                    record_values["genes"],
                    record_values["detail"],
                    record_values["depth"],
                    record_values["filt_ann"],
                    record_values["manta_n_occ"],
                    record_values["manta_t_occ"],
                    record_values["manta_n_af"],
                    record_values["manta_t_af"],
                    record_values["str_percent"],
                    record_values["pr_freq"],
                    record_values["sr_freq"],
                ]
                if sample_normal:
                    outline = outline + [record_values["pr_freq_n"], record_values["sr_freq_n"]]

                outline = outline + in_target
                manta_tables["bnd"]["data"].append(outline)

            elif (
                "MantaDEL" in record_values["id"]
                and record_values["svlength"] is not None
                and record_values["svlength"] <= -MIN_DEL_SIZE
            ):
                outline = [
                    str(record.contig),
                    int(record.pos),
                    record.stop,
                    record_values["svlength"],
                    record_values["id"],
                    record_values["genes"],
                    record_values["detail"],
                    record_values["filt_ann"],
                    record_values["manta_n_occ"],
                    record_values["manta_t_occ"],
                    record_values["manta_n_af"],
                    record_values["manta_t_af"],
                    record_values["str_percent"],
                    record_values["pr_freq"],
                    record_values["sr_freq"],
                ]
                if sample_normal:
                    outline = outline + [record_values["pr_freq_n"], record_values["sr_freq_n"]]

                outline = outline + in_target
                manta_tables["del"]["data"].append(outline)

            elif "MantaDUP" in record_values["id"]:
                outline = [
                    str(record.contig),
                    int(record.pos),
                    record.stop,
                    record_values["svlength"],
                    record_values["id"],
                    record_values["genes"],
                    record_values["detail"],
                    record_values["hom_len"],
                    record_values["hom_seq"],
                    record_values["filt_ann"],
                    record_values["manta_n_occ"],
                    record_values["manta_t_occ"],
                    record_values["manta_n_af"],
                    record_values["manta_t_af"],
                    record_values["str_percent"],
                    record_values["pr_freq"],
                    record_values["sr_freq"],
                ]
                if sample_normal:
                    outline = outline + [record_values["pr_freq_n"], record_values["sr_freq_n"]]

                outline = outline + in_target
                manta_tables["dup"]["data"].append(outline)

            elif "MantaINS" in record_values["id"]:
                outline = [
                    str(record.contig),
                    int(record.pos),
                    record.ref,
                    record.alts[0],
                    record_values["svlength"],
                    record_values["id"],
                    record_values["genes"],
                    record_values["detail"],
                    record_values["hom_len"],
                    record_values["hom_seq"],
                    record_values["filt_ann"],
                    record_values["manta_n_occ"],
                    record_values["manta_t_occ"],
                    record_values["manta_n_af"],
                    record_values["manta_t_af"],
                    record_values["str_percent"],
                    record_values["pr_freq"],
                    record_values["sr_freq"],
                ]
                if sample_normal:
                    outline = outline + [record_values["pr_freq_n"], record_values["sr_freq_n"]]

                outline = outline + in_target
                manta_tables["ins"]["data"].append(outline)

    return manta_tables
