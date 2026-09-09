# 🐶 🐶 🐶 🎶 WGS Leukemia Fluffy 🐍

Snakemake workflow to analyse hematological malignancies in whole genome data

![snakefmt](https://github.com/clinical-genomics-uppsala/fluffy_hematology_wgs/actions/workflows/snakefmt.yaml/badge.svg?branch=develop)
![snakemake dry run](https://github.com/clinical-genomics-uppsala/fluffy_hematology_wgs/actions/workflows/snakemake-dry-run.yaml/badge.svg?branch=develop)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## 💬 Introduction

This snakemake workflow uses modules from [hydra-genetics](https://github.com/hydra-genetics/) to process `.fastq` files and call SNVs, indels, CNVs and SVs from whole genome DNA, plus fusions from RNA. Alongside diagnosis-filtered `.vcf` files, the workflow produces produces Excel reports for manual review, a MultiQC report `.html` file and CNV plots. One of the modules contains the **commercial**
[parabricks toolkit](https://docs.nvidia.com/clara/parabricks/3.7.0/index.html) which can be replaced by sentieon or opensource GATK tools if required.

The workflow runs in two modes depending on what `units.tsv` provides for a sample:

- **Tumor/normal (`tn`)** — matched tumor and normal DNA. Somatic calling subtracts the patient's own germline, and
  Manta's somatic scoring is available.
- **Tumor-only (`t`)** — DNA only, no matched normal. Somatic/germline separation relies on population allele
  frequency and the normal-panel annotations described below.

The following hydra-genetics modules are part of this pipeline: `alignment`, `annotation`, `cnv_sv`, `compression`,
`filtering`, `fusions`, `misc`, `parabricks` (or `sentieon`), `prealignment`, `qc`, `reports`.

## ❗️ Dependencies

In order to use this module, the following dependencies are required:

[![hydra-genetics](https://img.shields.io/badge/hydragenetics-3.1.1-blue)](https://github.com/hydra-genetics/)
[![snakemake](https://img.shields.io/badge/snakemake-7.32-blue)](https://snakemake.readthedocs.io/en/stable/)
[![python](https://img.shields.io/badge/python-3.9-blue)](https://www.python.org/)

## 🎒 Preparations

### Sample and unit data

Input data should be added to
[`samples.tsv`](https://github.com/clinical-genomics-uppsala/fluffy_hematology_wgs/blob/develop/config/samples.tsv)
and
[`units.tsv`](https://github.com/clinical-genomics-uppsala/fluffy_hematology_wgs/blob/develop/config/units.tsv).
The following information need to be added to these files:

| Column Id                 | Description                                                                                      |
| ------------------------- | ------------------------------------------------------------------------------------------------ |
| **`samples.tsv`** |                                                                                                  |
| sample                    | unique sample/patient id, one per row                                                            |
| tumor_content             | ratio of tumor cells to total cells                                                              |
| **`units.tsv`**   |                                                                                                  |
| sample                    | same sample/patient id as in`samples.tsv`                                                      |
| type                      | data type identifier (one letter), can be one of**T**umor, **N**ormal, **R**NA |
| platform                  | type of sequencing platform, e.g.`NovaSeq`                                                     |
| machine                   | specific machine id, e.g. NovaSeq instruments have`@Axxxxx`                                    |
| flowcell                  | identifer of flowcell used                                                                       |
| lane                      | flowcell lane number                                                                             |
| barcode                   | sequence library barcode/index, connect forward and reverse indices by`+`, e.g. `ATGC+ATGC`  |
| fastq1/2                  | absolute path to forward and reverse reads                                                       |
| adapter                   | adapter sequences to be trimmed, separated by comma                                              |

A sample with only a `T` unit is run tumor-only; a sample with both `T` and `N` units is run tumor/normal.

### Reference data

Reference files should be specified in
[`config.yaml`](https://github.com/clinical-genomics-uppsala/fluffy_hematology_wgs/blob/develop/config/config.yaml)

1. A `.fasta` reference file of the human genome is required as well as an `.fai` file and an bwa index of this
   file.
2. A `.vcf` file containing known indel sites. For GRCh38, this file is available as part of the Broad GATK
   resource bundle at
   [google cloud](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz).
3. An `.interval_list` file containing all whole genome calling regions. The GRCh38 version is also available at
   [google cloud](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list).
4. The `trimmer_software` should be specified by indicating a rule which should be used for trimming. This
   pipeline uses `fastp_pe`.
5. `.bed` files defining regions of interest for different diagnoses. This pipeline is assuming `ALL` and `AML`
   and different gene lists for SNVs and SVs.
6. For pindel, a `.bed` file containing the region that the analysis should be limited to.
7. [simple_sv_annotation](https://github.com/AstraZeneca-NGS/simple_sv_annotation) comes with panel and a fusion
   pair list which should also be included in the `config.yaml`.
8. Annotation with [SnpEff](http://pcingola.github.io/SnpEff/) a database is needed which can be downloaded through
   the cli.
9. For [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), a cache resource should be downloaded prior
   to running the workflow.
10. A bgzipped and tabix-indexed `simple_repeats` `.bed` file (with a `.tbi` next to it), used to flag
    Manta calls overlapping simple repeats.
11. Optionally, `target_genes`: a text file with one gene symbol per line. Manta calls hitting these genes are flagged `In Target Panel` in the SV report.
12. `svdb_query.db_manta`: one or more SVDB databases (comma-separated) of Manta SVs from a normal-sample panel, used to flag recurrent/artefactual calls. Without this, SV calls are not checked against a normal panel.

## 🚀 Usage

To run the workflow,
[`resources.yaml`](https://github.com/clinical-genomics-uppsala/fluffy_hematology_wgs/blob/develop/config/resources.yaml)
is needed which defines different resources as default and for different rules. For parabricks, the `gres`
stanza is needed and should specify the number of GPUs available. You also need a [`config.yaml`](https://github.com/clinical-genomics-uppsala/fluffy_hematology_wgs/blob/develop/config/config.yaml) where all run-variables are defined.

```bash
snakemake --profile my-profile --configfile config/config.yaml
```

To run the integration test you only need to add lines in the `tests/integration/config.yaml` that differs from the original `config.yaml`. As of now it is only a dryrun test, no small dataset is available.

```bash
cd .tests/integration/
snakemake --snakefile ../../workflow/Snakefile --configfiles ../../config/config.yaml config.yaml -n
```

### Output files

.fastq files are archived as compressed file pair as .spring: `Archive/{project}/{sample}_{flowcell}_{lane}_{barcode}_{type}.spring`

The MultiQC html report can be found here: `Results/MultiQC_TN.html`

All results (as described in table below) are located in: `Results/{project}/{sample}/`

| File                                                                 | Description                                                                                                                                                         |
| -------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Cram/{sample}_{type}.crumble.cram`                                | crumbled`.cram` file                                                                                                                                              |
| `Cram/{sample}_{type}.crumble.cram.crai`                           | index for crumbled`.cram` file                                                                                                                                    |
| `SNV_indels/{sample}_T.vep.vcf.gz`                                 | `.vcf` output for SNV and small indels annotated with VEP for tumor_only                                                                                          |
| `SNV_indels/{sample}_T.vep.vcf.gz.tbi`                             | index for`.vcf` output from VEP for tumor_only                                                                                                                    |
| `SNV_indels/{sample}_TN.vep.vcf.gz`                                | `.vcf` output from VEP for tumor/normal                                                                                                                           |
| `SNV_indels/{sample}_TN.vep.vcf.gz.tbi`                            | index for`.vcf` output from VEP for tumor/normal                                                                                                                  |
| `SNV_indels/{sample}_TN.vep.all.vcf.gz`                            | `.vcf` output from VEP for tumor/normal, hard-filtered for ALL genes                                                                                              |
| `SNV_indels/{sample}_TN.vep.all.vcf.gz.tbi`                        | index for`.vcf` output from VEP for tumor/normal, hard-filtered for ALL genes                                                                                     |
| `SNV_indels/{sample}_TN.vep.aml.vcf.gz`                            | `.vcf` output from VEP for tumor/normal, hard-filtered for AML genes                                                                                              |
| `SNV_indels/{sample}_TN.vep.aml.vcf.gz.tbi`                        | index for`.vcf` output from VEP for tumor/normal, hard-filtered for AML genes                                                                                     |
| `SNV_indels/{sample}_mutectcaller_TN.all.tsv`                      | `.tsv` file for excel containing SNVs and small indels from mutect2 for ALL                                                                                       |
| `SNV_indels/{sample}_mutectcaller_TN.aml.tsv`                      | `.tsv` file for excel containing SNVs and small indels from mutect2 for AML                                                                                       |
| `SNV_indels/{sample}.pindel.vcf.gz`                                | `.vcf` output from pindel                                                                                                                                         |
| `SNV_indels/{sample}.pindel.vcf.gz.tbi`                            | index for`.vcf` output from pindel                                                                                                                                |
| `CNV/{sample}_T.vcf.gz`                                            | `.vcf` output from cnvkit                                                                                                                                         |
| `CNV/{sample}_T.vcf.gz.tbi`                                        | index for`.vcf` output from cnvkit                                                                                                                                |
| `CNV/{sample}_{type}.CNV.xlsx`                                     | Excel file containing overview of CNVkit results                                                                                                                    |
| `CNV/{sample}_T.png`                                               | scatter plot from cnvkit for entire genome                                                                                                                          |
| `CNV/{sample}_T_chr{chr}.png`                                      | scatter plot per chromosome from cnvkit                                                                                                                             |
| `SV/{sample}_manta_T.ssa.svdb_query.str_annotated.vcf.gz`          | `.vcf` output from Manta, tumor-only, annotated with the normal-panel SVDB frequency (`manta_N_AF`/`manta_N_OCC`) and simple-repeat overlap (`STR_PERCENT`) |
| `SV/{sample}_manta_T.ssa.svdb_query.str_annotated.vcf.gz.tbi`      | index for the file above                                                                                                                                            |
| `SV/{sample}_manta_T.ssa.svdb_query.str_annotated.all.vcf.gz`      | as above, filtered for ALL genes                                                                                                                                    |
| `SV/{sample}_manta_T.ssa.svdb_query.str_annotated.all.vcf.gz.tbi`  | index for the file above                                                                                                                                            |
| `SV/{sample}_manta_T.ssa.svdb_query.str_annotated.aml.vcf.gz`      | as above, filtered for AML genes                                                                                                                                    |
| `SV/{sample}_manta_T.ssa.svdb_query.str_annotated.aml.vcf.gz.tbi`  | index for the file above                                                                                                                                            |
| `SV/{sample}_manta_TN.ssa.svdb_query.str_annotated.vcf.gz`         | as above, tumor/normal                                                                                                                                              |
| `SV/{sample}_manta_TN.ssa.svdb_query.str_annotated.vcf.gz.tbi`     | index for the file above                                                                                                                                            |
| `SV/{sample}_manta_TN.ssa.svdb_query.str_annotated.all.vcf.gz`     | as above, filtered for ALL genes                                                                                                                                    |
| `SV/{sample}_manta_TN.ssa.svdb_query.str_annotated.all.vcf.gz.tbi` | index for the file above                                                                                                                                            |
| `SV/{sample}_manta_TN.ssa.svdb_query.str_annotated.aml.vcf.gz`     | as above, filtered for AML genes                                                                                                                                    |
| `SV/{sample}_manta_TN.ssa.svdb_query.str_annotated.aml.vcf.gz.tbi` | index for the file above                                                                                                                                            |
| `SV/{sample}_T.manta.xlsx`                                         | Excel report of structural variants from Manta, tumor-only                                                                                                          |
| `SV/{sample}_TN.manta.xlsx`                                        | Excel report of structural variants from Manta, tumor/normal                                                                                                        |
| `DNA_fusions/{sample}_pelops.txt`                                  | DUX4-IGH rearrangement evidence from Pelops                                                                                                                         |

#### The Manta Excel report

`SV/{sample}_{T,TN}.manta_new.xlsx` has one sheet per SV type (Deletions, Insertions, Duplications, Translocations)
plus panel-restricted translocation sheets, and an Overview sheet summarising target-gene hits and known fusions.
A few things to know before reading it:

- **Rows with a high normal-panel allele frequency are hidden, not removed.** Any call with `manta_N_AF` above the
  configured threshold is hidden by default (visible via "unhide" in Excel) unless it also matches a
  `target_genes` entry, in which case it stays visible and highlighted. On real samples this commonly hides close
  to half the rows on the Deletions and Insertions sheets
- **MaxDepth calls are rescued selectively, not dropped or kept wholesale.** Manta flags breakpoints in very
  high-depth regions as `MaxDepth`, which are enriched for mapping artefacts. Rather than filtering all of them
  out (losing real variants) or keeping all of them (drowning the report in noise), an event is rescued and shown
  — marked `MaxDepth Rescue` — only if it also has read support, no other blocking filter, no normal-panel hit, and
  isn't on a junk/alternate contig. See `workflow/scripts/manta_maxdepth_rescue.py` for the exact criteria.
- Coloring reflects `manta_N_AF`: orange means "seen in the normal panel above the visible-row threshold".

#### General Statistics - DNA

The general statistics table are ordered based on the "s-index" in fastq-filename. This is done by renaming the samples in two steps using the script `sample_order_multiqc.py`. To toggle between "Sample Order" and "Sample Name" use the buttons just above General Stats header.

<br />

| Column Name                         | Origin                                                               | Comment                                                                                                                                                                |
| ----------------------------------- | -------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| M Reads                             | [Picard](https://broadinstitute.github.io/picard/) HSMetrics          | Total number of reads in inputfile (`alignment/samtools_merge_bam/{sample}_{type}.bam`)                                                                              |
| % Mapped                            | [Samtools stats](http://www.htslib.org/doc/samtools-stats.html)       | Only reads on target (`config[reference][design_bed]`)                                                                                                               |
| % Proper pairs                      | [Samtools stats](http://www.htslib.org/doc/samtools-stats.html)       | Only reads on target (`config[reference][design_bed]`)                                                                                                               |
| Average Quality                     | [Samtools stats](http://www.htslib.org/doc/samtools-stats.html)       | Ratio between sum of base quality over total length. Only reads on target (`config[reference][design_bed]`)                                                          |
| Median                              | [Mosdepth](https://github.com/brentp/mosdepth)                        | Median Coverage over reference                                                                                                                                         |
| >= 10X                              | [Mosdepth](https://github.com/brentp/mosdepth)                        | Fraction of reference with coverage over 10x                                                                                                                           |
| >= 30X                              | [Mosdepth](https://github.com/brentp/mosdepth)                        | Fraction of reference with coverage over 30x                                                                                                                           |
| >=50X                               | [Mosdepth](https://github.com/brentp/mosdepth)                        | Fraction of reference with coverage over 50x                                                                                                                           |
| Error sex check                     | [Peddy](https://github.com/brentp/peddy)                              | Result of sex check based on sex in`units.tsv`                                                                                                                       |
| Predicted sex sex check             | [Peddy](https://github.com/brentp/peddy)                              |                                                                                                                                                                        |
| Bases on Target                     | [Picard](https://broadinstitute.github.io/picard/) HSMetrics          | Bases inside the capture design (`config[reference][design_intervals]`)                                                                                              |
| Fold80                              | [Picard](https://broadinstitute.github.io/picard/) HSMetrics          | The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets (`config[reference][design_intervals]`) |
| % Dups                              | [Picard](https://broadinstitute.github.io/picard/) DuplicationMetrics |                                                                                                                                                                        |
| Mean Insert Size                    | [Picard](https://broadinstitute.github.io/picard/) InsertSizeMetrics  |                                                                                                                                                                        |
| Target Bases with zero coverage [%] | [Picard](https://broadinstitute.github.io/picard/) HSMetrics          | Percent target (`config[reference][design_intervals]`) bases with 0 coverage                                                                                         |
| % Adapter                           | [fastp](https://github.com/OpenGene/fastp)                            |                                                                                                                                                                        |

### Program versions

default container: `docker://hydragenetics/common:0.1.9`

| Program                 | Version    | Container                                                                     |
| ----------------------- | ---------- | ----------------------------------------------------------------------------- |
| Arriba                  | 2.3.0      | `docker://hydragenetics/arriba:2.3.0`                                       |
| bbduk                   | 38.98      | `docker://kincekara/bbduk:38.98`                                            |
| CNVkit                  | 0.9.9      | `docker://hydragenetics/cnvkit:0.9.9` `docker://python:3.9.9-slim-buster` |
| Crumble                 | 0.8.3      | `docker://hydragenetics/crumble:0.8.3`                                      |
| fastp                   | 0.20.1     | `docker://hydragenetics/fastp:0.20.1`                                       |
| FastQC                  | 0.11.9     | `docker://hydragenetics/fastqc:0.11.9`                                      |
| FusionCatcher           | 1.33       | `docker://blcdsdockerregistry/fusioncatcher:1.33`                           |
| GATK                    | 4.2.2.0    | `docker://hydragenetics/gatk4:4.2.2.0`                                      |
| Manta                   | 1.6.0      | `docker://hydragenetics/manta:1.6.0`                                        |
| Mosdepth                | 0.3.2      | `docker://hydragenetics/mosdepth:0.3.2`                                     |
| MultiQC                 | 1.21       | `docker://hydragenetics/multiqc:1.11`                                       |
| Parabricks              | 4.5.1-1    | `docker://nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1`                    |
| Pelops                  | 0.8.0      | `docker://hydragenetics/pelops:0.8.0`                                       |
| Peddy                   | 0.4.8      | `docker://hydragenetics/peddy:0.4.8`                                        |
| Picard                  | 2.25.0     | `docker://hydragenetics/picard:2.25.0`                                      |
| Pindel                  | 0.2.5b9    | `docker://hydragenetics/pindel:0.2.5b9`                                     |
| RSeQC                   | 4.0.0      | `docker://hydragenetics/rseqc:4.0.0`                                        |
| simple_sv_annotation.py | 2019.02.18 | `docker://hydragenetics/simple_sv_annotation:2019.02.18`                    |
| snpEff                  | 5.0        | `docker://hydragenetics/snpeff:5.0`                                         |
| SortMeRNA               | 4.3.4      | `docker://hydragenetics/sortmerna:4.3.4`                                    |
| SPRING                  | 1.0.1      | `docker://hydragenetics/spring:1.0.1`                                       |
| STAR                    | 2.7.10a    | `docker://hydragenetics/star:2.7.10a`                                       |
| STAR-Fusion             | 1.10.1     | `docker://trinityctat/starfusion:1.10.1`                                    |
| svdb                    | 2.6.0      | `docker://hydragenetics/svdb:2.6.0`                                         |
| VEP                     | 111        | `docker://hydragenetics/vep:111`                                            |
| vt                      | 2015.11.10 | `docker://hydragenetics/vt:2015.11.10`                                      |

## :judge: Rule Graph Parabricks version

![rule_graph](images/rulegraph.svg)

## :judge: Rule Graph Sentieon version
