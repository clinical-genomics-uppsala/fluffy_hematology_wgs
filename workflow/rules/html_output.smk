__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2023, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


rule merge_cnv_json_chr:
    input:
        json=get_json_for_merge_cnv_json_chr,
        fai=config.get("reference", {}).get("fai", ""),
        annotation_bed=list(config.get("annotate_cnv", {}).values()),
        germline_vcf="parabricks/pbrun_mutectcaller_t/{sample}_{type}.normalized.vep.filter.germline.fix_af.vcf",
        cnv_vcfs=get_unfiltered_cnv_vcfs_for_merge_json,
        filtered_cnv_vcfs=get_filtered_cnv_vcfs_for_merge_json,
        filtered_cnv_vcfs_tbi=get_filtered_cnv_vcfs_tbi_for_merge_json,
        cnv_vcfs_tbi=get_unfiltered_cnv_vcfs_tbi_for_merge_json,
    output:
        json=temp("reports/cnv_html_report/{sample}_{type}.{tc_method}.{locus}.merged.json"),
    params:
        skip_chromosomes=lambda wildcards: [
            chromosome
            for chromosome in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
            if chromosome != wildcards.locus
        ],
    log:
        "reports/cnv_html_report/{sample}_{type}.{tc_method}.{locus}.merged.json.log",
    benchmark:
        repeat(
            "reports/cnv_html_report/{sample}_{type}.{tc_method}.{locus}.merged.json.benchmark.tsv",
            config.get("merge_cnv_json", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("merge_cnv_json", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("merge_cnv_json", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge_cnv_json", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge_cnv_json", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge_cnv_json", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge_cnv_json", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge_cnv_json", {}).get("container", config["default_container"])
    message:
        "{rule}: Merge CNV JSON data for {wildcards.sample}_{wildcards.type} {wildcards.locus}"
    script:
        "../scripts/merge_cnv_json.py"
