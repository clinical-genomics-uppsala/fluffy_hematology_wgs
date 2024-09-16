__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2023, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


def check_if_tn(wildcards):
    if wildcards.analysis == "tn":
        vcfs = expand(
            "parabricks/pbrun_mutectcaller_{{analysis}}/{{sample_type}}.normalized.vep.ratio.filter.somatic.include.{bed}.vcf.gz",
            bed=["all", "aml", "tm"],
        )+ ["cnv_sv/pindel_vcf/{sample_type}_T.no_tc.vep_annotated.vcf"]
    else:
        vcfs = expand(
            "parabricks/pbrun_mutectcaller_{{analysis}}/{{sample_type}}.normalized.vep.filter.somatic.include.{bed}.vcf.gz",
            bed=["all", "aml", "tm"],
        ) + ["cnv_sv/pindel_vcf/{sample_type}.no_tc.vep_annotated.vcf"]
    return vcfs


rule export_to_xlsx:
    input:
        vcfs=check_if_tn,
        all_bed=config["bcftools_SNV"]["all"],
        aml_bed=config["bcftools_SNV"]["aml"],
        tm_bed=config["bcftools_SNV"]["tm"],
        pindel_bed=config["pindel_call"]["include_bed"],
    output:
        xlsx=temp("export_to_xlsx/{analysis}/{sample_type}.snvs.xlsx"),
    params:
        filterfile = config["filter_vcf"]["somatic"],
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
        "{rule}: merge {input.vcfs} into {output.xlsx}"
    script:
        "../scripts/export_to_xlsx.py"
