__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2026, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


rule bbduk:
    input:
        fq1="prealignment/merged/{sample}_{type}_fastq1.fastq.gz",
        fq2="prealignment/merged/{sample}_{type}_fastq2.fastq.gz",
        ref=config.get("bbduk", {}).get("fasta", ""),
    output:
        out=temp("prealignment/bbduk/{sample}_{type}.non-rrna.fq.gz"),
        stats=temp("prealignment/bbduk/{sample}_{type}.rrna.log"),
    params:
        extra=config.get("bbduk", {}).get("extra", ""),
        ref=get_bbduk_refs,
        kmer=config.get("bbduk", {}).get("kmer", "31"),
    log:
        "prealignment/bbduk/{sample}_{type}.non-rrna.fq.gz.log",
    benchmark:
        repeat(
            "prealignment/bbduk/{sample}_{type}.rrna.fq.gz.benchmark.tsv",
            config.get("bbduk", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("bbduk", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bbduk", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bbduk", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bbduk", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bbduk", {}).get("time", config["default_resources"]["time"]),
    threads: config.get("bbduk", {}).get("threads", config["default_resources"]["threads"])
    container:
        config.get("bbduk", {}).get("container", config["default_container"])
    message:
        "{rule}: identify contamination from {input.ref} in {input.fq1} and {input.fq2}"
    shell:
        "bbduk.sh "
        "in={input.fq1} "
        "in2={input.fq2} "
        "ref={params.ref} "
        "out={output.out} "
        "stats={output.stats} "
        "threads={threads} "
        "k={params.kmer} "
        "{params.extra} &> {log}"
