#!/bin/python3

import gzip
from pysam import VariantFile


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
    try:
        return_dict["af"] = float(record.samples[sample_tumor]["AF"][0])
    except KeyError:
        return_dict["af"] = int(record.samples[sample_tumor].get("AD")[1]) / sum(record.samples[sample_tumor].get("AD"))

    if sample_normal != "":
        return_dict["n_af"] = float(record.samples[sample_normal]["AF"][0])
    else:
        return_dict["n_af"] = ""

    try:
        return_dict["dp"] = int(record.samples[sample_tumor]["DP"])
    except KeyError:
        return_dict["dp"] = sum(record.samples[sample_tumor].get("AD"))

    if sample_normal != "":
        return_dict["n_dp"] = int(record.samples[sample_normal]["DP"])
    else:
        return_dict["n_dp"] = ""

    try:
        return_dict["svlen"] = int(record.info["SVLEN"])
    except KeyError:
        pass

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
        return_dict = dict.fromkeys(
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
        )

    return return_dict


def extract_manta_vcf_values(record, ann_index, simple_ann_index, sample_tumor, sample_normal=""):
    return_dict = {}
    return_dict["filt_ann"] = ",".join(record.filter.keys())

    try:
        return_dict["id"] = ":".join(record.id.split(":")[0:2])
    except KeyError:
        return_dict["id"] = ""

    genes = []
    details = []
    try:
        record.info["SIMPLE_ANN"]
    except KeyError:
        try:
            record.info["ANN"]
        except KeyError:
            return_dict["genes"] = "NA"
            return_dict["detail"] = "NA"
        else:
            for annotation in record.info["ANN"]:
                annotation_values = annotation.split("|")
                gene_name_id = (
                    annotation_values[ann_index.index("Gene_Name")] + "(" + annotation_values[ann_index.index("Gene_ID")] + ")"
                )
                if gene_name_id not in genes:
                    genes.append(gene_name_id)

            return_dict["genes"] = ", ".join(genes)
            return_dict["detail"] = ""
    else:
        for annotation in record.info["SIMPLE_ANN"]:
            annotation_values = annotation.split("|")
            gene_name_id = (
                annotation_values[simple_ann_index.index("GENE(s)")]
                + "("
                + annotation_values[simple_ann_index.index("TRANSCRIPT")]
                + ")"
            )
            if gene_name_id not in genes:
                genes.append(gene_name_id)
            detail = annotation_values[
                simple_ann_index.index("DETAIL (exon losses, KNOWN_FUSION, ON_PRIORITY_LIST, NOT_PRIORITISED)")
            ]
            if detail not in details:
                details.append(detail)

        return_dict["genes"] = ", ".join(genes)
        return_dict["detail"] = ", ".join(details)

    try:
        return_dict["depth"] = record.info["BND_DEPTH"]
    except KeyError:
        return_dict["depth"] = ""

    try:
        pr_values = record.samples[sample_tumor]["PR"]
    except KeyError:
        return_dict["pr_freq"] = ""
    else:
        pr_denominator, pr_numerator = pr_values
        return_dict["pr_freq"] = pr_numerator / (pr_denominator + pr_numerator) if pr_denominator + pr_numerator != 0 else None

    try:
        sr_values = record.samples[sample_tumor]["SR"]
    except KeyError:
        return_dict["sr_freq"] = ""
    else:
        sr_denominator, sr_numerator = sr_values
        return_dict["sr_freq"] = sr_numerator / (sr_denominator + sr_numerator) if sr_denominator + sr_numerator != 0 else None

    if sample_normal:
        try:
            pr_values_n = record.samples[sample_normal]["PR"]
        except KeyError:
            return_dict["pr_freq_n"] = ""
        else:
            pr_denominator, pr_numerator = pr_values_n
            return_dict["pr_freq_n"] = (
                pr_numerator / (pr_denominator + pr_numerator) if pr_denominator + pr_numerator != 0 else None
            )

        try:
            sr_values_n = record.samples[sample_normal]["SR"]
        except KeyError:
            return_dict["sr_freq_n"] = ""
        else:
            sr_denominator, sr_numerator = sr_values_n
            return_dict["sr_freq_n"] = (
                sr_numerator / (sr_denominator + sr_numerator) if sr_denominator + sr_numerator != 0 else None
            )

    try:
        return_dict["svlength"] = record.info["SVLEN"][0]
    except KeyError:
        return_dict["svlength"] = ""

    try:
        return_dict["hom_len"] = record.info["HOMLEN"][0]
    except KeyError:
        return_dict["hom_len"] = ""

    try:
        return_dict["hom_seq"] = record.info["HOMSEQ"][0]
    except KeyError:
        return_dict["hom_seq"] = ""

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
):
    vcf_file = VariantFile(vcf_input)
    sample_tumor = [x for x in list(vcf_file.header.samples) if x.endswith("_T")][0]
    if len(list(vcf_file.header.samples)) > 1:
        sample_normal = [x for x in list(vcf_file.header.samples) if x.endswith("_N")][0]
    else:
        sample_normal = None

    for header_row in vcf_file.header.records:
        if "ID=ANN," in str(header_row):
            ann_index = str(header_row).split("'")[1].strip().split(" | ")
        elif "ID=SIMPLE_ANN," in str(header_row):
            simple_ann_index = str(header_row).split("'")[1].strip().split(" | ")

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
        {"header": "BreakEnd"},
        {"header": "Genes"},
        {"header": "Details"},
        {"header": "Depth"},
        {"header": "Annotation"},
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
        {"header": "Paired-read freq"},
        {"header": "Spanning-read freq"},
    ]
    if sample_normal:
        manta_tables["bnd"]["headers"] = manta_tables["bnd"]["headers"] + [
            {"header": "Paired-read Normal Freq"},
            {"header": "Spanning-read Normal Freq"},
        ]
        manta_tables["del"]["headers"] = manta_tables["del"]["headers"] + [
            {"header": "Paired-read Normal Freq"},
            {"header": "Spanning-read Normal Freq"},
        ]
        manta_tables["dup"]["headers"] = manta_tables["dup"]["headers"] + [
            {"header": "Paired-read Normal Freq"},
            {"header": "Spanning-read Normal Freq"},
        ]
        manta_tables["ins"]["headers"] = manta_tables["ins"]["headers"] + [
            {"header": "Paired-read Normal Freq"},
            {"header": "Spanning-read Normal Freq"},
        ]
    for record in vcf_file.fetch():
        record_values = extract_manta_vcf_values(record, ann_index, simple_ann_index, sample_tumor, sample_normal)
        if not any(x in avoid_filterflags for x in record_values["filt_ann"].split(",")):
            if "MantaBND" in record_values["id"]:
                outline = [
                    str(record.contig),
                    int(record.pos),
                    record_values["id"],
                    str(record.alts[0])[1:-1],
                    record_values["genes"],
                    record_values["detail"],
                    record_values["depth"],
                    record_values["filt_ann"],
                    record_values["pr_freq"],
                    record_values["sr_freq"],
                ]
                if sample_normal:
                    outline = outline + [record_values["pr_freq_n"], record_values["sr_freq_n"]]
                manta_tables["bnd"]["data"].append(outline)
            elif "MantaDEL" in record_values["id"] and record_values["svlength"] <= -100:
                outline = [
                    str(record.contig),
                    int(record.pos),
                    record.stop,
                    record_values["svlength"],
                    record_values["id"],
                    record_values["genes"],
                    record_values["detail"],
                    record_values["filt_ann"],
                    record_values["pr_freq"],
                    record_values["sr_freq"],
                ]
                if sample_normal:
                    outline = outline + [record_values["pr_freq_n"], record_values["sr_freq_n"]]
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
                    record_values["pr_freq"],
                    record_values["sr_freq"],
                ]
                if sample_normal:
                    outline = outline + [record_values["pr_freq_n"], record_values["sr_freq_n"]]
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
                    record_values["pr_freq"],
                    record_values["sr_freq"],
                ]
                if sample_normal:
                    outline = outline + [record_values["pr_freq_n"], record_values["sr_freq_n"]]
                manta_tables["ins"]["data"].append(outline)
    return manta_tables
