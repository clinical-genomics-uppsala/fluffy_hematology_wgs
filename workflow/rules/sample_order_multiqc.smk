__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule sample_order_multiqc:
    output:
        replacement_rna="qc/multiqc/sample_replacement_RNA.tsv",
        order_rna="qc/multiqc/sample_order_RNA.tsv",
        replacement_dna="qc/multiqc/sample_replacement_DNA.tsv",
        order_dna="qc/multiqc/sample_order_DNA.tsv",
        dnanumber="qc/multiqc/DNA_number.table.tsv",
        rnanumber="qc/multiqc/RNA_number.table.tsv",
    params:
        filelist=[(u.sample, u.type, u.fastq1) for u in units.itertuples()],
    log:
        "qc/multiqc/sample_order.tsv.log",
    benchmark:
        repeat("qc/multiqc/sample_order.tsv.benchmark.tsv", config.get("sample_order_multiqc", {}).get("benchmark_repeats", 1))
    threads: config.get("sample_order_multiqc", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sample_order_multiqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sample_order_multiqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sample_order_multiqc", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sample_order_multiqc", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sample_order_multiqc", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sample_order_multiqc", {}).get("container", config["default_container"])
    message:
        "{rule}: Create a sample order tsv based on units.tsv."
    script:
        "../scripts/sample_order_multiqc.py"
