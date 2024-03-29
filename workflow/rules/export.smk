__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2023, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule export_to_xlsx:
    input:
        all="parabricks/pbrun_mutectcaller_{analysis}/{sample_type}.vep.include.all.vcf.gz",
        all_bed=config["bcftools_SNV"]["all"],
        aml="parabricks/pbrun_mutectcaller_{analysis}/{sample_type}.vep.include.aml.vcf.gz",
        aml_bed=config["bcftools_SNV"]["aml"],
        tm="parabricks/pbrun_mutectcaller_{analysis}/{sample_type}.vep.include.tm.vcf.gz",
        tm_bed=config["bcftools_SNV"]["tm"],
    output:
        xlsx=temp("export_to_xlsx/{analysis}/{sample_type}.snvs.xlsx"),
    params:
        extra=config.get("export_to_xlsx", {}).get("extra", ""),
    log:
        "export_to_xlsx/{analysis}/{sample_type}.snvs.xslx.log",
    benchmark:
        repeat(
            "export_to_xlsx/{analysis}/{sample_type}.snvs.xslx.benchmark.tsv",
            config.get("export_to_xlsx", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("export_to_xlsx", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("export_to_xlsx", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("export_to_xlsx", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("export_to_xlsx", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("export_to_xlsx", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("export_to_xlsx", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("export_to_xlsx", {}).get("container", config["default_container"])
    message:
        "{rule}: merge {input.all}, {input.aml}, {input.tm} into {output.xlsx}"
    script:
        "../scripts/export_to_xlsx.py"
