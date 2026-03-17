rule fix_af:
    input:
        vcf="parabricks/pbrun_mutectcaller_t/{sample}_{type}.normalized.vep.filter.germline.vcf",
    output:
        vcf=temp("parabricks/pbrun_mutectcaller_t/{sample}_{type}.normalized.vep.filter.germline.fix_af.vcf"),
    log:
        "parabricks/pbrun_mutectcaller_t/{sample}_{type}.normalized.vep.filter.germline.fix_af.log",
    benchmark:
        repeat(
            "parabricks/pbrun_mutectcaller_t/{sample}_{type}.normalized.vep.filter.germline.fix_af.benchmark.tsv",
            config.get("fix_af", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fix_af", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fix_af", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fix_af", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fix_af", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fix_af", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fix_af", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fix_af", {}).get("container", config["default_container"])
    message:
        "{rule}: fix missing allele frequency field in format column in {input.vcf}"
    script:
        "../scripts/fix_af.py"
