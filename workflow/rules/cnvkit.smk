__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2021, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_call_no_custom_purity:
    input:
        segment="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        vcf=lambda wildcards: get_vcf_input(wildcards),
    output:
        segment=temp("cnv_sv/cnvkit_call/{sample}_{type}.{tc_method}.loh.cns"),
    log:
        "cnv_sv/cnvkit_call/{sample}_{type}.{tc_method}.loh.cns.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_call/{sample}_{type}.{tc_method}.loh.cns.benchmark.tsv",
            config.get("cnvkit_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_call", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_call", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_call", {}).get("container", config["default_container"])
    message:
        "{rule}: Call cnvs with loh info into cnv_sv/cnvkit_call/{wildcards.sample}_{wildcards.type}.loh.cns"
    shell:
        "(cnvkit.py call {input.segment} "
        "-v {input.vcf} "
        "-o {output.segment}) &> {log}"


rule cnvkit_scatter:
    input:
        segments="cnv_sv/cnvkit_call/{sample}_{type}.pathology.loh.cns",
        segment_regions="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr",
        vcf=lambda wildcards: get_vcf_input(wildcards),
    output:
        plot=temp("cnv_sv/cnvkit_scatter/{sample}_{type}_{locus}.png"),
    params:
        gene=lambda wildcards: config.get("cnvkit_scatter", {}).get("cna", {}).get(wildcards.locus, ""),
        extra=config.get("cnvkit_scatter", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_scatter/{sample}_{type}_{locus}.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_scatter/{sample}_{type}_{locus}.benchmark.tsv",
            config.get("cnvkit_scatter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_scatter", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_scatter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_scatter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_scatter", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_scatter", {}).get("container", config["default_container"])
    message:
        "{rule}: Plot cnvs into cnv_sv/cnvkit_scatter/{wildcards.sample}_{wildcards.type}_{wildcards.locus}.png"
    shell:
        """
        if [ -z {params.gene} ]
        then
        cnvkit.py scatter -s {input.segments} {input.segment_regions}  -v {input.vcf}  -o {output.plot} {params.extra}
        else
        cnvkit.py scatter -s {input.segments} {input.segment_regions}  -v {input.vcf}  -o {output.plot} -g {params.gene} {params.extra}
        fi
        """
