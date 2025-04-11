__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2021, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


if aligner == "bwa_gpu": 

    rule somalier_extract:
        input:
            sites=config.get("somalier_extract").get("sites"),
            fasta=config.get("reference").get("fasta"),
            sample="parabricks/pbrun_fq2bam/{sample}_{type}.bam",
        output:
            "qc/somalier/cohort/{sample}_{type}.somalier",
        log:
            "qc/somalier/log/{sample}_{type}.cohort.log",
        params:
            extract_folder="qc/somalier/cohort/",
        benchmark:
            repeat(
                "qc/somalier/log/{sample}_{type}.cohort.benchmark.tsv",
                config.get("somalier_extract", {}).get("benchmark_repeats", 1),
            )
        threads:
            config.get("somalier_extract", {}).get("threads", config["default_resources"]["threads"]),
        resources:
            threads=config.get("somalier_extract", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("somalier_extract", {}).get("time", config["default_resources"]["time"]),
            mem_mb=config.get("somalier_extract", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("somalier_extract", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            partition=config.get("somalier_extract", {}).get("partition", config["default_resources"]["partition"]),
        container:
            config.get("somalier_extract", {}).get("container", config["default_container"])
        message:
            "{rule}: extracts sites for somalier in sample {wildcards.sample}_{wildcards.type}.bam"
        shell:
            "somalier extract -s {input.sites} -f {input.fasta} -d {params.extract_folder} {input.sample}"


rule somalier_relate:
    input:
        samples=[
                "qc/somalier/cohort/%s_%s.somalier" % (sample,t)
                for sample in get_samples(samples)
                for t in get_unit_types(units, sample)
		if t in ["N", "T"]                
		],
    output:
        pairs="qc/somalier/somalier_relate.pairs.tsv",
        groups="qc/somalier/somalier_relate.groups.tsv",
        samples="qc/somalier/somalier_relate.samples.tsv",
        html="qc/somalier/somalier_relate.html",
    wildcard_constraints:
        type="T|N",
    log:
        "qc/somalier/somalier_relate.log",
    params:
        extra=config.get("somalier_relate", {}).get("extra", ""),
        outname="qc/somalier/somalier_relate",
    benchmark:
        repeat(
            "qc/somalier/somalier_relate.becnhmark.tsv",
            config.get("somalier_relate", {}).get("benchmark_repeats", 1),
        )
    threads:
        config.get("somalier_relate", {}).get("threads", config["default_resources"]["threads"]),
    resources:
        threads=config.get("somalier_relate", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("somalier_relate", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("somalier_relate", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("somalier_relate", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("somalier_relate", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("somalier_relate", {}).get("container", config["default_container"])
    message:
        "{rule}: Running somalier relate for inferring sex and checking T/N"
    shell:
        "somalier relate {params.extra} -o {params.outname} {input.samples}"


rule somalier_tn_test: 
    input:
        pairs="qc/somalier/somalier_relate.pairs.tsv",
    output:
        tncheck="qc/somalier/TMmismatch.txt",
    log:
        "qc/somalier/TMmismatch.log",
    params:
        extra=config.get("somalier_tn_test", {}).get("extra", ""),
    benchmark:
        repeat(
            "qc/somalier/somalier_tn_test.benchmark.tsv",
            config.get("somalier_tn_test", {}).get("benchmark_repeats", 1),
        )
    threads:
        config.get("somalier_tn_test", {}).get("threads", config["default_resources"]["threads"]),
    resources:
        threads=config.get("somalier_tn_test", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("somalier_tn_test", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("somalier_tn_test", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("somalier_tn_test", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("somalier_tn_test", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("somalier_tn_test", {}).get("container", config["default_container"])
    message:
        "{rule}: using awk to extract matched T/N samples from somalier that are not from the same individual"
    shell:
        """
        awk -F"[\t_]" '$1==$3 && $5<=0.2 {{print $1}}' {input.pairs} >{output.tncheck} &> {log}
        """
