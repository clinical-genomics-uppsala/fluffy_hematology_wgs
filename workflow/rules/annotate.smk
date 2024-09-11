__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule annotate_normal_ratio:
    input:
        vcf="parabricks/pbrun_mutectcaller_tn/{sample}.normalized.vep.vcf",
    output:
        vcf="parabricks/pbrun_mutectcaller_tn/{sample}.normalized.vep.ratio.vcf",
    log:
        "parabricks/pbrun_mutectcaller_tn/{sample}.normalized.vep.ratio.vcf.gz.log",
    benchmark:
        repeat(
            "parabricks/pbrun_mutectcaller_tn/{sample}.normalized.vep.ratio.vcf.gz.benchmark.tsv",
            config.get("annotate_normal_ratio", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("annotate_normal_ratio", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("annotate_normal_ratio", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("annotate_normal_ratio", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("annotate_normal_ratio", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("annotate_normal_ratio", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("annotate_normal_ratio", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("annotate_normal_ratio", {}).get("container", config["default_container"])
    message:
        "{rule}: Add normal ratio annotation to {input.vcf}"
    script:
        "../scripts/annotate_normal_ratio.py"
