__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2023, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


aligner = config.get("aligner", None)


if aligner == "bwa_gpu":

    rule dux4_igh_read_count:
        input:
            bam="parabricks/pbrun_fq2bam/{sample}_{type}.bam",
        output:
            cnt=temp("reports/dux_read_counts/{sample}_{type}.dux4_igh.txt"),
        log:
            "reports/dux_read_counts/{sample}_{type}.dux4_igh.txt.log",
        benchmark:
            repeat(
                "reports/dux_read_counts/{sample}_{type}.dux4_igh.txt.benchmark.tsv",
                config.get("fix_af", {}).get("benchmark_repeats", 1),
            )
        params:
            extra=config.get("dux4_igh", {}).get("extra", ""),
            region=config.get("dux4_igh", {}).get("region", ""),
        threads: config.get("dux4_igh", {}).get("threads", config["default_resources"]["threads"])
        resources:
            mem_mb=config.get("fdux4_igh", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("dux4_igh", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            partition=config.get("dux4_igh", {}).get("partition", config["default_resources"]["partition"]),
            threads=config.get("dux4_igh", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("dux4_igh", {}).get("time", config["default_resources"]["time"]),
        container:
            config.get("dux4_igh", {}).get("container", config["default_container"])
        message:
            "{rule}: get read count of paired reads mapping to DUX4 and IGH {input.bam}"
        shell:
            """
            samtools view -F 1024 -c -e \
            '(rnext == "4" && pnext > 190173000 && pnext < 190176000) || ([SA] =~ "4,19017[345][0-9]{{3}}")' \
            {input.bam} {params.region} > {output.cnt}
            """


if aligner == "bwa_gpu":

    rule dux4_erg_read_count:
        input:
            bam="parabricks/pbrun_fq2bam/{sample}_{type}.bam",
        output:
            cnt=temp("reports/dux_read_counts/{sample}_{type}.dux4_erg.txt"),
        log:
            "reports/dux_read_counts/{sample}_{type}.dux4_erg.txt.log",
        benchmark:
            repeat(
                "reports/dux_read_counts/{sample}_{type}.dux4_erg.txt.benchmark.tsv",
                config.get("fix_af", {}).get("benchmark_repeats", 1),
            )
        params:
            extra=config.get("dux4_erg", {}).get("extra", ""),
            region=config.get("dux4_erg", {}).get("region", ""),
        threads: config.get("dux4_erg", {}).get("threads", config["default_resources"]["threads"])
        resources:
            mem_mb=config.get("fdux4_erg", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("dux4_erg", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            partition=config.get("dux4_erg", {}).get("partition", config["default_resources"]["partition"]),
            threads=config.get("dux4_erg", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("dux4_erg", {}).get("time", config["default_resources"]["time"]),
        container:
            config.get("dux4_erg", {}).get("container", config["default_container"])
        message:
            "{rule}: get read count of paired reads mapping to DUX4 and ERG {input.bam}"
        shell:
            """
            samtools view -F 1024 -c -e \
            '(rnext == "4" && pnext > 190173000 && pnext < 190176000) || ([SA] =~ "4,19017[345][0-9]{{3}}")' \
            {input.bam} {params.region} > {output.cnt}
            """


rule dux4_igh_fusioncatcher_count:
    input:
        fq="fusions/fusioncatcher/{sample}_{type}/final-list_candidate-fusion-genes.txt",
    output:
        cnts=temp("fusions/fusioncatcher/{sample}_{type}/dux4_counts.txt"),
        hits=temp("fusions/fusioncatcher/{sample}_{type}/dux4_hits.txt"),
    log:
        "fusions/fusioncatcher/{sample}_{type}/dux4-igh.log",
    benchmark:
        repeat(
            "fusions/fusioncatcher/{sample}_{type}/dux4-igh.benchmark.tsv",
            config.get("dux4-igh_fusioncatcher", {}).get("benchmark_repeats", 1),
        )
    params:
        extra=config.get("dux4-igh_fusioncatcher", {}).get("extra", ""),
        region=config.get("dux4-igh_fusioncatcher", {}).get("region", ""),
    threads: config.get("dux4-igh_fusioncatcher", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("dux4-igh_fusioncatcher", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("dux4-igh_fusioncatcher", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("dux4-igh_fusioncatcher", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("dux4-igh_fusioncatcher", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("dux4-igh_fusioncatcher", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("dux4-igh_fusioncatcher", {}).get("container", config["default_container"])
    message:
        "{rule}: get count and results of fusioncatcher translocations in DUX4 with IGH in {input.fq}"
    shell:
        """
        if grep -E "IGH@.*DUX4|DUX4.*IGH@|ERG.*DUX4|DUX4.*ERG" {input.fq}; then
        var=$(grep -E "IGH@.*DUX4|DUX4.*IGH@|ERG.*DUX4|DUX4.*ERG" {input.fq})
        echo "$var" > {output.hits}
        wc -l {output.hits} > {output.cnts}
        else
        echo "none" > {output.hits}
        echo "0" > {output.cnts}
        fi
        """
