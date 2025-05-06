__author__ = "Martin Rippin, Arielle R. Munters, Nina Hollfelder"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "arielle.munters@scilifelab.uu.se, nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


import itertools
import numpy as np
import pandas as pd
import pathlib
import re
from snakemake.utils import validate
from snakemake.utils import min_version
import yaml
from datetime import datetime

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics import min_version as hydra_min_version

from hydra_genetics.utils.misc import replace_dict_variables
from hydra_genetics.utils.misc import export_config_as_file
from hydra_genetics.utils.software_versions import add_version_files_to_multiqc
from hydra_genetics.utils.software_versions import add_software_version_to_config
from hydra_genetics.utils.software_versions import export_pipeline_version_as_file
from hydra_genetics.utils.software_versions import export_software_version_as_file
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import use_container
from hydra_genetics.utils.software_versions import touch_software_version_file
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name


hydra_min_version("3.0.0")
min_version("7.32.0")


### Set and validate config file
config = replace_dict_variables(config)

try:
    validate(config, schema="../schemas/config.schema.yaml")
except WorkflowError as we:
    # Probably a validation error, but the original exception in lost in
    # snakemake. Pull out the most relevant information instead of a potentially
    # *very* long error message.
    if not we.args[0].lower().startswith("error validating config file"):
        raise
    error_msg = "\n".join(we.args[0].splitlines()[:2])
    parent_rule_ = we.args[0].splitlines()[3].split()[-1]
    if parent_rule_ == "schema:":
        sys.exit(error_msg)
    else:
        schema_hiearachy = parent_rule_.split()[-1]
        schema_section = ".".join(re.findall(r"\['([^']+)'\]", schema_hiearachy)[1::2])
        sys.exit(f"{error_msg} in {schema_section}")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


## get version information on pipeline, containers and software
pipeline_name = "fluffy_hematology_wgs"
pipeline_version = get_pipeline_version(workflow, pipeline_name=pipeline_name)
version_files = touch_pipeline_version_file_name(
    pipeline_version, date_string=pipeline_name, directory="Results/versions/software"
)
if use_container(workflow):
    version_files.append(touch_software_version_file(config, date_string=pipeline_name, directory="Results/versions/software"))
add_version_files_to_multiqc(config, version_files)


onstart:
    export_pipeline_version_as_file(pipeline_version, date_string=pipeline_name, directory="Results/versions/software")
    if use_container(workflow):
        update_config, software_info = add_software_version_to_config(config, workflow, False)
        export_software_version_as_file(software_info, date_string=pipeline_name, directory="Results/versions/software")
    date_string = datetime.now().strftime("%Y%m%d")
    export_config_as_file(update_config, date_string=date_string, directory="Results/versions")


### Read and validate samples file
samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file
units = (
    pd.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


### Read output_files for cp rules
with open(config["output"]) as output:
    output_spec = yaml.safe_load(output.read())


### Set wildcard constraints
wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",
    bed="aml|all|tm",


def type_generator(types):
    if "N" in types and "T" in types:
        types.add("TN")
        return types
    else:
        return types


def get_bam_input(wildcards, t_n=None, use_sample_wildcard=True):
    if use_sample_wildcard is True and t_n is None:
        sample_str = "{}_{}".format(wildcards.sample, wildcards.type)
    elif use_sample_wildcard is True and t_n:
        sample_str = "{}_{}".format(wildcards.sample, t_n)
    else:
        sample_str = wildcards.file

    aligner = config.get("aligner", None)
    if aligner is None:
        sys.exit("aligner missing from config, valid options: bwa_gpu or bwa_sentieon")
    elif aligner == "bwa_gpu":
        bam_input = "parabricks/pbrun_fq2bam/{}.bam".format(sample_str)
    elif aligner == "bwa_sentieon":
        bam_input = "sentieon/realign/{}_REALIGNED.bam".format(sample_str)
    else:
        sys.exit("valid options for aligner are: bwa_gpu or bwa_sentieon")

    bai_input = "{}.bai".format(bam_input)

    return (bam_input, bai_input)


def get_num_gpus(rule, wildcards):
    gres = config.get(rule, {"gres": "--gres=gres:gpu:1"}).get("gres", "--gres=gres:gpu:1")[len("--gres=") :]
    gres_dict = dict()
    for item in gres.split(","):
        items = item.split(":")
        gres_dict[items[1]] = items[2]
    return gres_dict["gpu"]


def get_vcf_input(wildcards):
    aligner = config.get("aligner", None)
    if aligner is None:
        sys.exit("aligner missing from config, valid options: bwa_gpu or bwa_sentieon")
    elif aligner == "bwa_gpu":
        vcf_input = "parabricks/pbrun_mutectcaller_t/{}_{}.normalized.vep.filter.germline.vcf".format(
            wildcards.sample, wildcards.type
        )
    elif aligner == "bwa_sentieon":
        vcf_input = "sentieon/tnscope/{}_TNscope_tn_ML.vcf".format(wildcards.sample)
    else:
        sys.exit("valid options for aligner are: bwa_gpu or bwa_sentieon")

    return vcf_input


def get_cnv_callers(tc_method):
    for tcm in config.get("svdb_merge", {}).get("tc_method", []):
        if tcm["name"] == tc_method:
            return tcm["cnv_caller"]
    raise ValueError(f"no cnv caller config available for tc_method {tc_method}")


def get_json_for_merge_cnv_json_chr(wildcards):
    callers = get_cnv_callers(wildcards.tc_method)
    return [
        "reports/cnv_html_report/{sample}_{type}.{caller}.{tc_method}.{locus}.json".format(caller=c, **wildcards) for c in callers
    ]


def get_tc(wildcards):
    if wildcards.tc_method == "pathology":
        try:
            return get_sample(samples, wildcards)["tumor_content"]
        except KeyError:
            return None


def get_unfiltered_cnv_vcfs_for_merge_json(wildcards):
    cnv_vcfs = []
    tags = config.get("cnv_html_report", {}).get("cnv_vcf", [])
    for t in tags:
        cnv_vcfs.append(
            f"cnv_sv/svdb_query/{wildcards.sample}_{wildcards.type}.{wildcards.tc_method}.svdb_query."
            f"annotate_cnv.{t['annotation']}.vcf.gz"
        )
    return sorted(cnv_vcfs)


def get_filtered_cnv_vcfs_for_merge_json(wildcards):
    cnv_vcfs = []
    tags = config.get("cnv_html_report", {}).get("cnv_vcf", [])
    for t in tags:
        cnv_vcfs.append(
            f"cnv_sv/svdb_query/{wildcards.sample}_{wildcards.type}.{wildcards.tc_method}.svdb_query."
            f"annotate_cnv.{t['annotation']}.filter.{t['filter']}.vcf.gz"
        )
    return sorted(cnv_vcfs)


def get_unfiltered_cnv_vcfs_tbi_for_merge_json(wildcards):
    cnv_vcfs = []
    tags = config.get("cnv_html_report", {}).get("cnv_vcf", [])
    for t in tags:
        cnv_vcfs.append(
            f"cnv_sv/svdb_query/{wildcards.sample}_{wildcards.type}.{wildcards.tc_method}.svdb_query."
            f"annotate_cnv.{t['annotation']}.vcf.gz.tbi"
        )
    return sorted(cnv_vcfs)


def get_filtered_cnv_vcfs_tbi_for_merge_json(wildcards):
    cnv_vcfs = []
    tags = config.get("cnv_html_report", {}).get("cnv_vcf", [])
    for t in tags:
        cnv_vcfs.append(
            f"cnv_sv/svdb_query/{wildcards.sample}_{wildcards.type}.{wildcards.tc_method}.svdb_query."
            f"annotate_cnv.{t['annotation']}.filter.{t['filter']}.vcf.gz.tbi"
        )
    return sorted(cnv_vcfs)


def get_json_for_merge_cnv_json(wildcards):
    callers = get_cnv_callers(wildcards.tc_method)
    return ["reports/cnv_html_report/{sample}_{type}.{caller}.{tc_method}.json".format(caller=c, **wildcards) for c in callers]


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec["directory"])
    output_files = []
    output_fullpath = []

    chromosome_numbers = ["X", "Y"]
    chromosome_numbers.extend(range(1, 23))
    for filedef in output_spec["files"]:
        # add all output that is not TN
        output_files += set(
            [
                filedef["output"].format(
                    sample=sample,
                    type=unit_type,
                    chr=chromosome_number,
                    flowcell=flowcell,
                    barcode=barcode,
                    lane=lane,
                )
                for chromosome_number in chromosome_numbers
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
                if unit_type in set(filedef["types"])
                for flowcell in set([u.flowcell for u in units.loc[(sample, unit_type)].dropna().itertuples()])
                for barcode in set([u.barcode for u in units.loc[(sample, unit_type)].dropna().itertuples()])
                for lane in set([u.lane for u in units.loc[(sample, unit_type)].dropna().itertuples()])
            ]
        )

    # Iterate all files again and add all TN files for samples that have both T and N in units
    for filedef in output_spec["files"]:
        output_files += set(
            [
                filedef["output"].format(
                    sample=sample,
                    type=unit_type,
                    chr=chromosome_number,
                )
                for chromosome_number in chromosome_numbers
                for sample in get_samples(samples)
                for unit_type in type_generator(get_unit_types(units, sample))
                if unit_type in set(filedef["types"]) and unit_type == "TN"
            ]
        )
    # Add directory to beginning of each outputfile
    for op in output_files:
        output_fullpath.append(outdir / Path(op))

    return list(set(output_fullpath))


def generate_copy_rules(output_spec):
    output_directory = pathlib.Path(output_spec["directory"])
    rulestrings = []

    for f in output_spec["files"]:
        if f["input"] is None:
            continue

        rule_name = "copy_{}".format("_".join(re.sub(r"[\"'-.,]", "", f["name"].strip().lower()).split()))
        input_file = pathlib.Path(f["input"])
        output_file = output_directory / pathlib.Path(f["output"])

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        rule_code = "\n".join(
            [
                f'@workflow.rule(name="{rule_name}")',
                f'@workflow.input("{input_file}")',
                f'@workflow.output("{output_file}")',
                f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
                f'@workflow.container("{copy_container}")',
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb="{mem_mb}", '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
                '@workflow.shellcmd("cp --preserve=timestamps -r {input} {output}")',
                "@workflow.run\n",
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, "
                "log, version, rule, conda_env, container_img, singularity_args, use_singularity, "
                "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
                "shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):",
                '\tshell("(cp --preserve=timestamps -r {input[0]} {output[0]}) &> {log}", bench_record=bench_record, '
                "bench_iteration=bench_iteration)\n\n",
            ]
        )
        rulestrings.append(rule_code)

    exec(compile("\n".join(rulestrings), "copy_result_files", "exec"), workflow.globals)


generate_copy_rules(output_spec)
