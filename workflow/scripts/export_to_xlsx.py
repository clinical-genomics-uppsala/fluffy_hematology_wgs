#!/bin/python3

import xlsxwriter
import datetime
from pysam import VariantFile


def index_vep(variantfile):
    csq_index = []
    for x in variantfile.header.records:
        if "CSQ" in str(x):
            csq_index = str(x).split("Format: ")[1].strip().strip('">').split("|")
    return csq_index


""" Prepping input data """
bedfiles = {}
bedfiles["all"] = snakemake.input.all_bed
bedfiles["aml"] = snakemake.input.aml_bed
bedfiles["tm"] = snakemake.input.tm_bed

vcfs = {}
vcfs["all"] = [x for x in snakemake.input.vcfs if "include.all" in x][0]
vcfs["aml"] = [x for x in snakemake.input.vcfs if "include.aml" in x][0]
vcfs["tm"] = [x for x in snakemake.input.vcfs if "include.tm" in x][0]
subsections = ["all", "aml", "tm"]
vcf_tables = {}
for subsection in subsections:
    vcf_tables[subsection] = []
    vcf = VariantFile(vcfs[subsection])
    csq_index = index_vep(vcf)
    sample_tumor = [x for x in list(vcf.header.samples) if x.endswith("_T")][0]
    if len(list(vcf.header.samples)) > 1:
        sample_normal = [x for x in list(vcf.header.samples) if x.endswith("_N")][0]
    else:
        af_normal = ""
        dp_normal = ""

    for record in vcf.fetch():
        filterflag = ",".join(record.filter.keys())
        af = float(record.samples[sample_tumor].get("AF")[0])
        if af > 0.01:
            csq = record.info["CSQ"][0].split("|")
            gene = csq[csq_index.index("SYMBOL")]
            transcript = csq[csq_index.index("HGVSc")].split(":")[0]
            exon = csq[csq_index.index("EXON")]

            if len(csq[csq_index.index("HGVSc")].split(":")) > 1:
                coding_name = csq[csq_index.index("HGVSc")].split(":")[1]
            else:
                coding_name = ""
            ensp = csq[csq_index.index("HGVSp")]
            consequence = csq[csq_index.index("Consequence")]
            existing = csq[csq_index.index("Existing_variation")].split("&")

            cosmic_list = [cosmic for cosmic in existing if cosmic.startswith("CO")]
            if len(cosmic_list) == 0:
                cosmic = ""
            else:
                cosmic = ", ".join(cosmic_list)

            clinical = csq[csq_index.index("CLIN_SIG")]

            rs_list = [rs for rs in existing if rs.startswith("rs")]
            if len(rs_list) == 0:
                rs = ""
            else:
                rs = ", ".join(rs_list)
            max_pop_af = csq[csq_index.index("MAX_AF")]
            max_pops = csq[csq_index.index("MAX_AF_POPS")]

            if len(list(vcf.header.samples)) > 1:
                af_normal = float(record.samples[sample_normal].get("AF")[0])
                dp_normal = int(record.samples[sample_normal].get("DP"))

            outline = [
                filterflag,
                sample_tumor,
                gene,
                record.contig,
                int(record.pos),
                record.ref,
                record.alts[0],
                af,
                af_normal,
                int(record.info["DP"]),
                dp_normal,
                transcript,
                exon,
                coding_name,
                ensp,
                consequence,
                cosmic,
                clinical,
                rs,
                max_pop_af,
                max_pops,
            ]
            vcf_tables[subsection].append(outline)

# Adding bedfiles
bed_tables = {}
for subsection in subsections:
    bed_table = []
    with open(bedfiles[subsection], "r") as file:
        for line in file:
            bed_table.append(line.strip().split("\t"))
    bed_tables[subsection] = bed_table

""" Creating xlsx file """
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
worksheet_overview = workbook.add_worksheet("Overview")

heading_format = workbook.add_format({"bold": True, "font_size": 18})
bold_format = workbook.add_format({"bold": True, "text_wrap": True})

# Overview sheet
worksheet_overview.write(0, 0, sample_tumor, heading_format)
worksheet_overview.write(1, 0, "Processing date: " + datetime.datetime.now().strftime("%d %B, %Y"))

worksheet_overview.write(4, 0, "Created by: ")
worksheet_overview.write(4, 4, "Valid from: ")
worksheet_overview.write(5, 0, "Signed by: ")
worksheet_overview.write(5, 4, "Document nr: ")


worksheet_overview.write(7, 0, "Sheets:", bold_format)
worksheet_overview.write_url(8, 0, "internal:'ALL'!A1", string="Variants in ALL genes")
worksheet_overview.write_url(9, 0, "internal:'AML'!A1", string="Variants in AML genes")
worksheet_overview.write_url(10, 0, "internal:'TM'!A1", string="Variants in TM exons")
worksheet_overview.write_url(11, 0, "internal:'ALL bedfile'!A1", string="Gene regions included in ALL bedfile")
worksheet_overview.write_url(12, 0, "internal:'AML bedfile'!A1", string="Gene regions included in AML bedfile")
worksheet_overview.write_url(13, 0, "internal:'TM bedfile'!A1", string="Gene regions included in TM bedfile")

worksheet_overview.write(16, 0, "ALL bedfile: " + bedfiles["all"])
worksheet_overview.write(17, 0, "AML bedfile: " + bedfiles["aml"])
worksheet_overview.write(18, 0, "TM exons bedfile: " + bedfiles["tm"])

# Add snv variants sheets
for sheet in subsections:
    vcf_data = vcf_tables[sheet]
    worksheet = workbook.add_worksheet(sheet.upper())
    worksheet.set_column("B:B", 12)
    worksheet.set_column("E:E", 10)
    worksheet.set_column("L:L", 15)
    worksheet.write("A1", "Variants in " + sheet.upper() + " regions", heading_format)
    worksheet.write("A3", "Sample: " + str(sample_tumor))

    tableheading = [
        "FilterFlag",
        "DNAnr",
        "Gene",
        "Chr",
        "Pos",
        "Ref",
        "Alt",
        "AF",
        "Normal AF",
        "DP",
        "Normal DP",
        "Transcript",
        "Exon",
        "Mutation cds",
        "ENSP",
        "Consequence",
        "COSMIC ids on position",
        "Clinical significance",
        "dbSNP",
        "Max popAF",
        "Max Pop",
    ]
    heading_list = [{"header": a} for a in tableheading]

    worksheet.add_table(
        "A5:U" + str(len(vcf_data) + 6), {"data": vcf_data, "columns": heading_list, "style": "Table Style Light 1"}
    )
# Add bedfile sheets
for sheet in subsections:
    bed_data = bed_tables[sheet]
    worksheet = workbook.add_worksheet(sheet.upper() + " bedfile")
    worksheet.set_column("B:C", 10)
    worksheet.write("A1", "Bedfile used to filter out " + sheet.upper() + " regions", heading_format)
    tableheading = ["Chr", "Start", "End", "Region"]
    heading_list = [{"header": a} for a in tableheading]
    worksheet.add_table(
        "A4:D" + str(len(bed_data) + 5), {"data": bed_data, "columns": heading_list, "style": "Table Style Light 1"}
    )

workbook.close()
