__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2021, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"

import os


rule somalier_create_ped_T:
    input:
        config.get("samples", ""),
    output:
        [f"qc/somalier/{sample}_T.fam" for sample in get_samples(samples)],
    log:
        "qc/somalier/somalier_create_ped_T.fam.log",
    benchmark:
        repeat(
            "qc/somalier/somalier_create_ped_T.fam.benchmark.tsv",
            config.get("somalier_create_ped_T", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("somalier_create_ped_T", {}).get("threads", config.get("default_resources", {}).get("threads"))
    params:
        sample_type="T",
    resources:
        mem_mb=config.get("somalier_create_ped_T", {}).get("mem_mb", config.get("default_resources", {}).get("mem_mb")),
        mem_per_cpu=config.get("somalier_create_ped_T", {}).get("mem_per_cpu", config.get("default_resources", {}).get("mem_per_cpu")),
        partition=config.get("somalier_create_ped_T", {}).get("partition", config.get("default_resources", {}).get("partition")),
        threads=config.get("somalier_create_ped_T", {}).get("threads", config.get("default_resources", {}).get("threads")),
        time=config.get("somalier_create_ped_T", {}).get("time", config.get("default_resources", {}).get("time")),
    container:
        config.get("somalier_create_ped_T", {}).get("container", config.get("default_container", ""))
    message:
        "{rule}: Create fam file for all T samples for somalier input"
    script:
        "../scripts/somalier_create_ped.py"


rule somalier_create_ped_N:
    input:
        config.get("samples", ""),
    output:
        [f"qc/somalier/{sample}_N.fam" for sample in get_samples(samples)],
    log:
        "qc/somalier/somalier_create_ped_N.fam.log",
    benchmark:
        repeat(
            "qc/somalier/somalier_create_ped_N.fam.benchmark.tsv",
            config.get("somalier_create_ped_N", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("somalier_create_ped_N", {}).get("threads", config.get("default_resources", {}).get("threads"))
    params:
        sample_type="N",
    resources:
        mem_mb=config.get("somalier_create_ped_N", {}).get("mem_mb", config.get("default_resources", {}).get("mem_mb")),
        mem_per_cpu=config.get("somalier_create_ped_N", {}).get("mem_per_cpu", config.get("default_resources", {}).get("mem_per_cpu")),
        partition=config.get("somalier_create_ped_N", {}).get("partition", config.get("default_resources", {}).get("partition")),
        threads=config.get("somalier_create_ped_N", {}).get("threads", config.get("default_resources", {}).get("threads")),
        time=config.get("somalier_create_ped_N", {}).get("time", config.get("default_resources", {}).get("time")),
    container:
        config.get("somalier_create_ped_N", {}).get("container", config.get("default_container", ""))
    message:
        "{rule}: Create fam file for all T samples for somalier input"
    script:
        "../scripts/somalier_create_ped.py"


rule somalier_combine_fam:
    input:
        fam=[
            "qc/somalier/%s_%s.fam" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
            if t in ["N", "T"]
        ],
    output:
        ped="qc/somalier/somalier_all.ped",
    log:
        "qc/peddy/somalier_all.ped.log",
    benchmark:
        repeat(
            "qc/somalier/somalier_all.ped.benchmark.tsv",
            config.get("somalier_combine_fam", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("somalier_combine_fam", {}).get("threads", config.get("default_resources", {}).get("threads"))
    resources:
        mem_mb=config.get("somalier_combine_fam", {}).get("mem_mb", config.get("default_resources", {}).get("mem_mb")),
        mem_per_cpu=config.get("somalier_combine_fam", {}).get("mem_per_cpu", config.get("default_resources", {}).get("mem_per_cpu")),
        partition=config.get("somalier_combine_fam", {}).get("partition", config.get("default_resources", {}).get("partition")),
        threads=config.get("somalier_combine_fam", {}).get("threads", config.get("default_resources", {}).get("threads")),
        time=config.get("somalier_combine_fam", {}).get("time", config.get("default_resources", {}).get("time")),
    container:
        config.get("somalier_combine_fam", {}).get("container", config.get("default_container", ""))
    message:
        "{rule}: creates combined somalier_all.ped for sex check"
    shell:
        """
        cat {input.fam} > {output.ped}
        """


rule somalier_create_groupfile:
    input:
        samples=config.get("samples", ""),
        units=config.get("units", ""),
    output:
        "qc/somalier/somalier.groups",
    log:
        "qc/somalier/somalier.groups.log",
    benchmark:
        repeat(
            "qc/somalier/somalier.groups.benchmark.tsv",
            config.get("somalier_create_groupfile", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("somalier_create_groupfile", {}).get("threads", config.get("default_resources", {}).get("threads"))
    resources:
        mem_mb=config.get("somalier_create_groupfile", {}).get("mem_mb", config.get("default_resources", {}).get("mem_mb")),
        mem_per_cpu=config.get("somalier_create_groupfile", {}).get("mem_per_cpu", config.get("default_resources", {}).get("mem_per_cpu")),
        partition=config.get("somalier_create_groupfile", {}).get("partition", config.get("default_resources", {}).get("partition")),
        threads=config.get("somalier_create_groupfile", {}).get("threads", config.get("default_resources", {}).get("threads")),
        time=config.get("somalier_create_groupfile", {}).get("time", config.get("default_resources", {}).get("time")),
    container:
        config.get("somalier_create_groupfile", {}).get("container", config.get("default_container", ""))
    message:
        "{rule}: Create group file for somalier input"
    shell:
        """
        for i in $( cut -f1 {input.samples} | tail -n+2 )
        do
        var=$(grep $i {input.units} | cut -f2 | uniq | tr "\\n" "," | sed "s/,$/\\n/")
        if [ $var == "N,T" ] || [ $var == "T,N" ]
        then echo ${{i}}_N,${{i}}_T
        fi
        done > {output}
        """


rule somalier_custom_multiqc:
    input:
        conf=config.get("somalier_mqc", {}).get("config", ""),
        samples="qc/somalier/somalier_relate.samples.tsv",
    output:
        "qc/somalier/somalier_samples_mqc.tsv",
    log:
        "qc/somalier/somalier_custom_multiqc.log",
    params:
        extra=config.get("somalier_custom_multiqc", {}).get("extra", ""),
    benchmark:
        repeat(
            "qc/somalier/somalier_custom_multiqc.becnhmark.tsv",
            config.get("somalier_custom_multiqc", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("somalier_custom_multiqc", {}).get("threads", config.get("default_resources", {}).get("threads"))
    resources:
        threads=config.get("somalier_custom_multiqc", {}).get("threads", config.get("default_resources", {}).get("threads")),
        time=config.get("somalier_custom_multiqc", {}).get("time", config.get("default_resources", {}).get("time")),
        mem_mb=config.get("somalier_custom_multiqc", {}).get("mem_mb", config.get("default_resources", {}).get("mem_mb")),
        mem_per_cpu=config.get("somalier_custom_multiqc", {}).get("mem_per_cpu", config.get("default_resources", {}).get("mem_per_cpu")),
        partition=config.get("somalier_custom_multiqc", {}).get("partition", config.get("default_resources", {}).get("partition")),
    container:
        config.get("somalier_custom_multiqc", {}).get("container", config.get("default_container", ""))
    message:
        "{rule}: creating custom input for somalier to MultiQC general stats"
    script:
        "../scripts/somalier_mqc_config.py"


if aligner == "bwa_gpu":

    rule somalier_extract:
        input:
            sites=config.get("somalier_extract").get("sites"),
            fasta=config.get("reference").get("fasta"),
            sample="parabricks/pbrun_fq2bam_recal/{sample}_{type}.bam",
        output:
            "qc/somalier/cohort/{sample}_{type}.somalier",
        log:
            "qc/somalier/log/{sample}_{type}.cohort.log",
        params:
            extract_folder=lambda wildcards, output: os.path.dirname(output[0]),
        benchmark:
            repeat(
                "qc/somalier/log/{sample}_{type}.cohort.benchmark.tsv",
                config.get("somalier_extract", {}).get("benchmark_repeats", 1),
            )
        threads: config.get("somalier_extract", {}).get("threads", config.get("default_resources", {}).get("threads"))
        resources:
            threads=config.get("somalier_extract", {}).get("threads", config.get("default_resources", {}).get("threads")),
            time=config.get("somalier_extract", {}).get("time", config.get("default_resources", {}).get("time")),
            mem_mb=config.get("somalier_extract", {}).get("mem_mb", config.get("default_resources", {}).get("mem_mb")),
            mem_per_cpu=config.get("somalier_extract", {}).get("mem_per_cpu", config.get("default_resources", {}).get("mem_per_cpu")),
            partition=config.get("somalier_extract", {}).get("partition", config.get("default_resources", {}).get("partition")),
        container:
            config.get("somalier_extract", {}).get("container", config.get("default_container", ""))
        message:
            "{rule}: extracts sites for somalier in sample {wildcards.sample}_{wildcards.type}.bam"
        shell:
            "somalier extract -s {input.sites} -f {input.fasta} -d {params.extract_folder} {input.sample}"


rule somalier_relate:
    input:
        samples=[
            "qc/somalier/cohort/%s_%s.somalier" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
            if t in ["N", "T"]
        ],
        ped="qc/somalier/somalier_all.ped",
        group="qc/somalier/somalier.groups",
    output:
        pairs="qc/somalier/somalier_relate.pairs.tsv",
        samples="qc/somalier/somalier_relate.samples.tsv",
        html="qc/somalier/somalier_relate.html",
    log:
        "qc/somalier/somalier_relate.log",
    params:
        extra=config.get("somalier_relate", {}).get("extra", ""),
        outname=lambda wildcards, output: output.html.replace(".html", ""),
    benchmark:
        repeat(
            "qc/somalier/somalier_relate.benchmark.tsv",
            config.get("somalier_relate", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("somalier_relate", {}).get("threads", config.get("default_resources", {}).get("threads"))
    resources:
        threads=config.get("somalier_relate", {}).get("threads", config.get("default_resources", {}).get("threads")),
        time=config.get("somalier_relate", {}).get("time", config.get("default_resources", {}).get("time")),
        mem_mb=config.get("somalier_relate", {}).get("mem_mb", config.get("default_resources", {}).get("mem_mb")),
        mem_per_cpu=config.get("somalier_relate", {}).get("mem_per_cpu", config.get("default_resources", {}).get("mem_per_cpu")),
        partition=config.get("somalier_relate", {}).get("partition", config.get("default_resources", {}).get("partition")),
    container:
        config.get("somalier_relate", {}).get("container", config.get("default_container", ""))
    message:
        "{rule}: Running somalier relate for inferring sex and checking T/N"
    shell:
        "somalier relate {params.extra} --ped {input.ped} -g {input.group} -o {params.outname} {input.samples}"


rule somalier_tn_test:
    input:
        pairs="qc/somalier/somalier_relate.pairs.tsv",
    output:
        tncheck="qc/somalier/TNmismatch.txt",
    log:
        "qc/somalier/TNmismatch.log",
    benchmark:
        repeat(
            "qc/somalier/somalier_tn_test.benchmark.tsv",
            config.get("somalier_tn_test", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("somalier_tn_test", {}).get("threads", config.get("default_resources", {}).get("threads"))
    resources:
        threads=config.get("somalier_tn_test", {}).get("threads", config.get("default_resources", {}).get("threads")),
        time=config.get("somalier_tn_test", {}).get("time", config.get("default_resources", {}).get("time")),
        mem_mb=config.get("somalier_tn_test", {}).get("mem_mb", config.get("default_resources", {}).get("mem_mb")),
        mem_per_cpu=config.get("somalier_tn_test", {}).get("mem_per_cpu", config.get("default_resources", {}).get("mem_per_cpu")),
        partition=config.get("somalier_tn_test", {}).get("partition", config.get("default_resources", {}).get("partition")),
    container:
        config.get("somalier_tn_test", {}).get("container", config.get("default_container", ""))
    message:
        "{rule}: using awk to extract matched T/N samples from somalier that are not from the same individual"
    shell:
        """
        awk -F"[\t_]" '$1==$3 && $5<=0.2 {{print $1}}' {input.pairs} >{output.tncheck} &> {log}
        """
