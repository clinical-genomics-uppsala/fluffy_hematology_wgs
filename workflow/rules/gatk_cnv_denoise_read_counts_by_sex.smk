__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2023, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


aligner = config.get("aligner", None)


if aligner == "bwa_gpu":

    rule gatk_cnv_collect_allelic_counts:
        input:
            bam="parabricks/pbrun_fq2bam/{sample}_T.bam",
            bai="parabricks/pbrun_fq2bam/{sample}_T.bam.bai",
            interval=config.get("gatk_cnv_collect_allelic_counts", {}).get("SNP_interval", ""),
            ref=config["reference"]["fasta"],
        output:
            temp("cnv_sv/gatk_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv"),
        params:
            extra=config.get("gatk_cnv_collect_allelic_counts", {}).get("extra", ""),
        log:
            "cnv_sv/gatk_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv.log",
        benchmark:
            repeat(
                "cnv_sv/gatk_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv.benchmark.tsv",
                config.get("gatk_cnv_collect_allelic_counts", {}).get("benchmark_repeats", 1),
            )
        threads: config.get("gatk_cnv_collect_allelic_counts", {}).get("threads", config["default_resources"]["threads"])
        resources:
            threads=config.get("gatk_cnv_collect_allelic_counts", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("gatk_cnv_collect_allelic_counts", {}).get("time", config["default_resources"]["time"]),
            mem_mb=config.get("gatk_cnv_collect_allelic_counts", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("gatk_cnv_collect_allelic_counts", {}).get(
                "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
            ),
            partition=config.get("gatk_cnv_collect_allelic_counts", {}).get("partition", config["default_resources"]["partition"]),
        container:
            config.get("gatk_cnv_collect_allelic_counts", {}).get("container", config["default_container"])
        message:
            "{rule}: Use gatk_cnv to obtain cnv_sv/gatk_collect_allelic_counts/{wildcards.sample}_{wildcards.type}.clean.allelicCounts.tsv"
        shell:
            "(gatk --java-options '-Xmx32g' CollectAllelicCounts "
            "-L {input.interval} "
            "-I {input.bam} "
            "-R {input.ref} "
            "-O {output}"
            "{params.extra}) &> {log}"


rule gatk_cnv_denoise_read_counts_by_sex:
    input:
        hdf5PoN_f=config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get("pon_female", ""),
        hdf5PoN_m=config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get("pon_male", ""),
        hdf5Tumor="cnv_sv/gatk_collect_read_counts/{sample}_{type}.counts.hdf5",
        sex="qc/somalier/somalier_relate.samples.tsv",
    output:
        denoisedCopyRatio=temp("cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv"),
        stdCopyRatio=temp("cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.standardizedCR.tsv"),
    params:
        extra=config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get("extra", ""),
    log:
        "cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv.benchmark.tsv",
            config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("gatk_cnv_denoise_read_counts_by_sex", {}).get("container", config["default_container"])
    message:
        "{rule}: Use gatk_cnv to obtain cnv_sv/gatk_cnv_denoise_read_counts/{wildcards.sample}_{wildcards.type}.clean.denoisedCR.tsv"
    shell:
        """
        if [[ $(awk -v st="{{wildcards.sample}}_{{wildcards.type}}" '($2 == st) {{print $7}}' {input.sex}) == "male" ]]
        then 
        (gatk --java-options '-Xmx7g' DenoiseReadCounts \
        -I {input.hdf5Tumor} \
        --count-panel-of-normals {input.hdf5PoN_m} \
        --standardized-copy-ratios {output.stdCopyRatio} \
        --denoised-copy-ratios {output.denoisedCopyRatio} \
        {params.extra}) &> {log}
        elif [[ $(awk -v st="{{wildcards.sample}}_{{wildcards.type}}" '($2 == st) {{print $7}}' {input.sex}) == "female" ]]
        then
        (gatk --java-options '-Xmx7g' DenoiseReadCounts \
        -I {input.hdf5Tumor} \
        --count-panel-of-normals {input.hdf5PoN_f} \
        --standardized-copy-ratios {output.stdCopyRatio} \
        --denoised-copy-ratios {output.denoisedCopyRatio} \
        {params.extra}) &> {log}
        else
        if [[ $(awk -v st="{{wildcards.sample}}_{{wildcards.type}}" '($2 == st) {{print $5}}' {input.sex}) == "1" ]]
        then 
        (gatk --java-options '-Xmx7g' DenoiseReadCounts \
        -I {input.hdf5Tumor} \
        --count-panel-of-normals {input.hdf5PoN_m} \
        --standardized-copy-ratios {output.stdCopyRatio} \
        --denoised-copy-ratios {output.denoisedCopyRatio} \
        {params.extra}) &> {log}
        else 
        (gatk --java-options '-Xmx7g' DenoiseReadCounts \
        -I {input.hdf5Tumor} \
        --count-panel-of-normals {input.hdf5PoN_f} \
        --standardized-copy-ratios {output.stdCopyRatio} \
        --denoised-copy-ratios {output.denoisedCopyRatio} \
        {params.extra}) &> {log}
        fi
        fi
        """


rule gatk_model_segments:
    input:
        denoisedCopyRatio="cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv",
        allelicCounts="cnv_sv/gatk_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv",
    output:
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.seg"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.cr.seg"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.af.igv.seg"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.cr.igv.seg"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.hets.tsv"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelBegin.cr.param"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelBegin.af.param"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelBegin.seg"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.af.param"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.cr.param"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        outprefix="{sample}_{type}.clean",
        extra=config.get("gatk_model_segments", {}).get("extra", ""),
    log:
        "cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.seg.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.seg.benchmark.tsv",
            config.get("gatk_model_segments", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_model_segments", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("gatk_model_segments", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_model_segments", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("gatk_model_segments", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_model_segments", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_model_segments", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("gatk_model_segments", {}).get("container", config["default_container"])
    message:
        "{rule}: Use gatk_cnv to obtain cnv_sv/gatk_model_segments/{wildcards.sample}_{wildcards.type}.clean.modelFinal.seg"
    shell:
        "(gatk --java-options '-Xmx16g' ModelSegments "
        "--denoised-copy-ratios {input.denoisedCopyRatio} "
        "--allelic-counts {input.allelicCounts} "
        "--output {params.outdir} "
        "--output-prefix {params.outprefix}"
        "{params.extra}) &> {log}"
