__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2023, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


rule peddy_create_ped:
    input:
        config["samples"],
    output:
        [f"qc/peddy/{sample}.peddy.fam" for sample in get_samples(samples)],
    log:
        "qc/peddy/peddy_create_ped.fam.log",
    benchmark:
        repeat(
            "qc/peddy/peddy_create_ped.fam.benchmark.tsv",
            config.get("peddy_create_ped", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("peddy_create_ped", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("peddy_create_ped", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("peddy_create_ped", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("peddy_create_ped", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("peddy_create_ped", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("peddy_create_ped", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("peddy_create_ped", {}).get("container", config["default_container"])
    message:
        "{rule}: Create fam file for all samples for peddy input"
    script:
        "../scripts/peddy_create_ped.py"


rule combine_peddy_results:
    input:
        sex=lambda wildcards: [
            "qc/peddy/%s/peddy.sex_check.csv" % (s)
            for s in set([row["sample"] for index, row in units.iterrows() if row["type"] == "T"])
        ],
        fam=lambda wildcards: [
            "qc/peddy/%s.peddy.fam" % (s) for s in set([row["sample"] for index, row in units.iterrows() if row["type"] == "T"])
        ],
    output:
        peddy_sex_check="qc/peddy/peddy.sex_check.csv",
        ped="qc/peddy/all.ped",
    log:
        "qc/peddy/all.ped.log",
    benchmark:
        repeat(
            "qc/peddy/all.ped.benchmark.tsv",
            config.get("combine_peddy_results", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("combine_peddy_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("combine_peddy_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("combine_peddy_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("combine_peddy_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("combine_peddy_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("combine_peddy_results", {}).get("container", config["default_container"])
    message:
        "{rule}: creates combined all.ped and qc/peddy/peddy.sex_check.csv for T samples for sex check"
    shell:
        """
        cat {input.sex} | awk 'NR==1; !/^sample_id/' > {output.peddy_sex_check}
        cat {input.fam} > {output.ped}
        """


rule create_peddy_mqc_tsv:
    input:
        peddy_sex_check="qc/peddy/peddy.sex_check.csv",
        ped="qc/peddy/all.ped",
        yaml=config["peddy"]["config"],
    output:
        sex_check_mqc=temp("qc/peddy/peddy_sex_check_mqc.tsv"),
    log:
        "qc/peddy/peddy.log",
    benchmark:
        repeat(
            "qc/peddy/create_peddy_mqc_tsv.benchmark.tsv",
            config.get("create_peddy_mqc_tsv", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("create_peddy_mqc_tsv", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_peddy_mqc_tsv", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_peddy_mqc_tsv", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_peddy_mqc_tsv", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_peddy_mqc_tsv", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_peddy_mqc_tsv", {}).get("container", config["default_container"])
    message:
        "{rule}: Create multiqc custom content embedded config tsv files from peddy sex_check and ped_check files"
    script:
        "../scripts/create_peddy_mqc_config.py"
