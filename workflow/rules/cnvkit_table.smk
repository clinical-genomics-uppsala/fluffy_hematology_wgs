__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_table:
    input:
        cns="cnv_sv/cnvkit_call/{sample}_{type}.loh.cns",
        gene_interest=config["cnvkit_table"]["bedfile"],
    output:
        temp("cnv_sv/cnvkit_table/{sample}_{type}.CNV.xlsx"),
    params:
        extra=config.get("cnvkit_table", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_table/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_table/{sample}_{type}.output.benchmark.tsv",
            config.get("cnvkit_table", {}).get("benchmark_repeats", 1)
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
    conda:
        "../envs/cnvkit_table.yaml"
    message:
        "{rule}: Create output table from "
    script:
        "../scripts/cnvkit_to_table.py"