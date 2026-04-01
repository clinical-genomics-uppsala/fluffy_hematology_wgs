#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import yaml
import numpy as np
import sys


def comment_the_config_keys(config_dict):

    commented_config = '\n'.join(
        ['# ' + line for line in yaml.dump(config_dict).rstrip('\n').split('\n')])

    return commented_config


def get_sex_check_df(sex_check_file_path):
    sex_check_df = pd.read_csv(sex_check_file_path)
    sex_check_df["sex_check_test"] = np.where(~sex_check_df.error,
                                              'Pass', 'Fail')  # report peddy  error check as pass/fail for simplicity
    sex_check_df.rename(columns={'sample_id': 'Sample'}, inplace=True)

    return sex_check_df


def write_peddy_mqc(peddy_df, peddy_config, outfile):

    with open(outfile, 'w') as outfile:
        print(comment_the_config_keys(peddy_config), file=outfile)

        peddy_df.to_csv(outfile, sep='\t', mode='a', index=False)


def main():

    try:

        config = snakemake.config["peddy"]["config"]

        with open(config, 'r') as report_configs:
            peddy_mqc_configs = yaml.load(report_configs, Loader=yaml.FullLoader)

            sex_check_df = get_sex_check_df(snakemake.input.peddy_sex_check)
            peddy_sex_config = peddy_mqc_configs.get('peddy_sex_check')
            write_peddy_mqc(sex_check_df, peddy_sex_config, snakemake.output.sex_check_mqc)

    except FileNotFoundError:
        sys.exit('Path to peddy config file not found/specified in config.yaml')


if __name__ == '__main__':
    main()
