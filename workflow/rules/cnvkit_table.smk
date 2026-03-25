__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_table:
    input:
        cns="cnv_sv/cnvkit_call/{sample}_{type}.pathology.loh.cns",
        gene_interest=config["cnvkit_table"]["bedfile"],
        cnv_scatter=expand(
            "cnv_sv/cnvkit_scatter/{{sample}}_{{type}}_{locus}.png",
            locus=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"],
        ),
        cnvkit_scatter_whole="cnv_sv/cnvkit_scatter/{sample}_{type}.png",
        cyto=config["cnvkit_table"]["cyto_coordinates"],
    output:
        temp("cnv_sv/cnvkit_table/{sample}_{type}.CNV.xlsx"),
    params:
        cnvkit_scatter_folder=lambda wildcards, input: os.path.dirname(input.cnvkit_scatter_whole),
        log=config.get("cnvkit_table", {}).get("log_thresholds", "-0.25,0.2"),
        ploidy=config.get("cnvkit_table", {}).get("ploidy", "2"),
        tc=lambda wildcards: get_sample(samples, wildcards)["tumor_content"],
        extra=config.get("cnvkit_table", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_table/{sample}_{type}.CNV.xlsx.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_table/{sample}_{type}.CNV.xlsx.benchmark.tsv",
            config.get("cnvkit_table", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_table", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_table", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_table", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_table", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_table", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_table", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_table", {}).get("container", config["default_container"])
    message:
        "{rule}: Create output table from {input.cns} and {input.gene_interest}."
    script:
        "../scripts/cnvkit_to_table.py"
