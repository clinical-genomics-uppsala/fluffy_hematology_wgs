__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2026, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


rule pelops:
    input:
        bam="parabricks/pbrun_fq2bam_recal/{sample}_T.bam",
    output:
        json=temp("fusions/pelops/{sample}/{sample}.json"),
        sam=temp("fusions/pelops/{sample}/01_CoreDUX4-IGH.sam"),
    log:
        "fusions/pelops/{sample}/{sample}.log",
    benchmark:
        repeat(
            "fusions/pelops/{sample}/{sample}.benchmark.tsv",
            config.get("pelops", {}).get("benchmark_repeats", 1),
        )
    params:
        extra=config.get("pelops", {}).get("extra", ""),
        dir="fusions/pelops/{sample}/",
    threads: config.get("pelops", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pelops", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pelops", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pelops", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pelops", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pelops", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: calculate SRPB for DUX4 rearrangements with pelops in {input.bam}"
    shell:
        "(pelops dux4r "
        "--export {params.dir} "
        "--json {output.json} "
        "{params.extra} {input.bam}) &> {log}"


rule extract_srpb:
    input:
        json="fusions/pelops/{sample}/{sample}.json",
    output:
        tsv=temp("fusions/pelops/{sample}/{sample}.tsv"),
    log:
        "fusions/pelops/{sample}/{sample}.tsv.log",
    benchmark:
        repeat(
            "fusions/pelops/{sample}/{sample}.tsv.benchmark.tsv",
            config.get("extract_srpb", {}).get("benchmark_repeats", 1),
        )
    params:
        extra=config.get("extract_srpb", {}).get("extra", ""),
    threads: config.get("extract_srpb", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("extract_srpb", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("extract_srpb", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("extract_srpb", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("extract_srpb", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("extract_srpb", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("extract_srpb", {}).get("container", config["default_container"])
    message:
        "{rule}: extract SRPB for DUX4 - IGH rearrangements {input.json} and save to a tsv file"
    script:
        "../scripts/extract_srpb.py"

