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
def extract_vcf_values(record, csq_index, sample_tumor, sample_normal = ""):
    return_dict = {}
    csq = record.info["CSQ"][0].split("|")

    try:
        return_dict["af"] = float(record.samples[sample_tumor]["AF"][0])
    except KeyError:
        return_dict["af"] = int(record.samples[sample_tumor].get("AD")[1]) / sum(record.samples[sample_tumor].get("AD"))
    
    if sample_normal != "":
        return_dict["n_af"] = float(record.samples[sample_normal]["AF"][0])
    else:
        return_dict["n_af"] = ""

    try:
        return_dict["dp"] = int(record.info["DP"])
    except KeyError:
        return_dict["dp"] = sum(record.samples[sample_tumor].get("AD"))

    try:
        return_dict["svlen"] = int(record.info["SVLEN"])
    except KeyError:
        pass

    # try:
    #     return_dict["artifact_callers"] = (
    #         str(record.info["Artifact"]).replace("(", "").replace(")", "").replace(", ", ";").replace("'", "")
    #     )
    #     if type(record.info["ArtifactMedian"]) == str:
    #         return_dict["artifact_median"] = str(round(float(record.info["ArtifactMedian"]), 3))
    #     else:
    #         return_dict["artifact_median"] = ";".join([str(round(float(x), 3)) for x in record.info["ArtifactMedian"]])
    #     return_dict["artifact_nr_sd"] = str(record.info["ArtifactNrSD"]).replace("(", "").replace(")", "").replace(", ", ";")
    # except KeyError:
    #     pass

    # try:
    #     return_dict["background_median"] = record.info["PanelMedian"]
    #     return_dict["background_nr_sd"] = record.info["PositionNrSD"]
    # except KeyError:
    #     return_dict["background_median"] = ""
    #     return_dict["background_nr_sd"] = ""

    # try:
    #     return_dict["callers"] = ";".join(record.info["CALLERS"])
    # except KeyError:
    #     pass

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
    return_dict["filter_flag"] = ",".join(record.filter.keys())

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



