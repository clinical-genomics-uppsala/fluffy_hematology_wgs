#!/bin/python3.6

import sys
import xlsxwriter


def float_or_na(value):
    return float(value) if value != '' else None


# import bedfile for CNA genes from bedfile.
bedtable = []
bed_header = ['chr', 'start', 'end', 'annotation']
with open(snakemake.input.gene_interest) as bedfile:
    for line in bedfile:
        bedtable.append(line.strip().split('\t'))

# Import genomic pos to cytoCoord translation to list
chr_bands = []
with open(snakemake.input.cyto, 'r') as chr_band_file:
    for line in chr_band_file:
        chr_bands.append(line.split("\t"))

chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
relevant_cnvs = {i: [] for i in chromosomes}
relevant_cnvs_header = [
    'ChromosomeBands',
    'Chromosome',
    'Start',
    'End',
    'Log2',
    'BAF',
    "CopyNumber",
    'CopyAllel1',
    'CopyAllel2',
    'Depth',
    'Probes',
    'Weight'
]
# import pdb; pdb.set_trace()
with open(snakemake.input.cns, 'r+') as cnsfile:
    cns_header = next(cnsfile).rstrip().split("\t")
    for cnv_line in cnsfile:
        cnv = cnv_line.strip().split("\t")
        if not (((int(snakemake.params.ploidy) % 2) == 0 and int(snakemake.params.ploidy) == int(cnv[cns_header.index('cn')])
                 and cnv[cns_header.index('cn1')] == cnv[cns_header.index('cn2')])
                or ((int(snakemake.params.ploidy) % 2) != 0 and int(snakemake.params.ploidy) == int(cnv[cns_header.index('cn')])
                    and cnv[cns_header.index('cn1')]+1 == cnv[cns_header.index('cn2')])):
            cnv_chr = cnv[cns_header.index('chromosome')]
            cnv_start = int(cnv[cns_header.index('start')])
            cnv_end = int(cnv[cns_header.index('end')])
            cnv_baf = float_or_na(cnv[cns_header.index('baf')])

            # Translate genomic coordinate to cytogen coordinates
            cyto_coord = ['', '']
            for chr_band in chr_bands:
                if chr_band[0] == cnv_chr:
                    if (cnv_start >= int(chr_band[1]) and cnv_start <= int(chr_band[2])):
                        cyto_coord[0] = chr_band[3]
                    if (cnv_end >= int(chr_band[1]) and cnv_end <= int(chr_band[2])):
                        cyto_coord[1] = chr_band[3]
            if cyto_coord[0] == cyto_coord[1]:
                cyto_coordinates = cnv_chr[3:]+cyto_coord[0]
            else:
                cyto_coordinates = cnv_chr[3:]+cyto_coord[0]+'-'+cyto_coord[1]

            if (cnv_end - cnv_start) >= 100000:
                outline = [cyto_coordinates, cnv_chr, cnv_start, cnv_end, float(cnv[cns_header.index('log2')]), cnv_baf,
                           cnv[cns_header.index('cn')], cnv[cns_header.index('cn1')], cnv[cns_header.index('cn2')],
                           cnv[cns_header.index('depth')], cnv[cns_header.index('probes')], cnv[cns_header.index('weight')]]
                relevant_cnvs[cnv_chr].append(outline)
                continue
            else:
                for bedline in bedtable:
                    if (cnv_chr == bedline[bed_header.index('chr')]):
                        bed_start = int(bedline[bed_header.index('start')])
                        bed_end = int(bedline[bed_header.index('end')])
                        if ((cnv_start >= bed_start and cnv_end <= bed_end) or
                           (cnv_start >= bed_start and cnv_start <= bed_end) or
                           (cnv_end >= bed_start and cnv_end <= bed_end) or
                           (cnv_start <= bed_start and cnv_end >= bed_end)):
                            outline = [cyto_coordinates, cnv_chr, cnv_start, cnv_end, float(cnv[cns_header.index('log2')]),
                                       cnv_baf, cnv[cns_header.index('cn')], cnv[cns_header.index('cn1')],
                                       cnv[cns_header.index('cn2')], cnv[cns_header.index('depth')],
                                       cnv[cns_header.index('probes')], cnv[cns_header.index('weight')]]
                            relevant_cnvs[cnv_chr].append(outline)
                            break

''' Creating xlsx file '''
low_log = float(snakemake.params.log.split(',')[0])
high_log = float(snakemake.params.log.split(',')[1])
sample = str(snakemake.input.cns).split("/")[-1].split("_")[0]
workbook = xlsxwriter.Workbook(snakemake.output[0])
heading_format = workbook.add_format({'bold': True, 'font_size': 18})
tablehead_format = workbook.add_format({'bold': True, 'text_wrap': True})
red_format = workbook.add_format({'font_color': 'red', 'italic': True})
worksheet_whole = workbook.add_worksheet("Whole genome")
worksheet_whole.set_column('B:E', 10)
worksheet_whole.write('A1', 'CNVkit calls', heading_format)
worksheet_whole.write('A3', 'Sample: ' + str(sample))
worksheet_whole.write('A4', 'Tumor content: ' + str(snakemake.params.tc))
worksheet_whole.write('A6', 'Calls larger then 100kb or in CNA bedfile included')
worksheet_whole.write('A7', 'CNV calls with log2 values inbetween ' + str(low_log) + ' and ' + str(high_log) +
                      ' are to be considered uncertain.', red_format)
image_path = snakemake.params.cnvkit_scatter_folder + '/' + sample + '_T'+'.png'
worksheet_whole.insert_image('A9', image_path)
worksheet_whole.write('A31', 'CNA bedfile: '+snakemake.input.gene_interest)
# create sheets for each chromosome
for chromosome in chromosomes:
    worksheet = workbook.add_worksheet(chromosome)
    worksheet.set_column('C:F', 10)
    worksheet.write('A1', 'CNVkit calls for '+chromosome, heading_format)
    worksheet.write('A3', 'Sample: '+str(sample))
    worksheet.write('A5', 'Calls larger then 100kb or in CNA bedfile included')
    worksheet.write('A6', 'CNV calls with log2 values inbetween ' + str(low_log) + ' and ' + str(high_log) +
                    ' are to be considered uncertain.', red_format)
    image_path = snakemake.params.cnvkit_scatter_folder + '/' + sample + '_T_'+chromosome+'.png'
    worksheet.insert_image('A8', image_path)
    worksheet.write_row('A30', relevant_cnvs_header, tablehead_format)
    row = 30
    col = 0
    for line in relevant_cnvs[chromosome]:
        if (low_log < float(line[4]) < high_log):
            worksheet.write_row(row, col, line, red_format)
        else:
            worksheet.write_row(row, col, line)
        row += 1

workbook.set_size(1800, 1200)
workbook.close()
