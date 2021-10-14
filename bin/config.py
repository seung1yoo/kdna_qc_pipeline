

import os
import sys
import yaml
import json
import subprocess

import pandas as pd
from collections import defaultdict


class Utils:
    def __init__(self, config):
        self.config = config

    def get_targets(self, targets):
        ls = list()
        for sample, lane_dic in self.config['samples'].items():
            for lane, info_dic in lane_dic.items():
                ls.append(f"analysis/bwa/{sample}/{sample}_{lane}.bam")
                ls.append(f"analysis/bwa/{sample}/{sample}_{lane}.sort.bam")
                ls.append(f"analysis/bwa/{sample}/{sample}_{lane}.dedup.bam")
                ls.append(f"analysis/bwa/{sample}/{sample}_{lane}.dedup.bai")
                ls.append(f"analysis/bwa/{sample}/{sample}_{lane}.dedup.txt")
            ls.append(f"analysis/bwa/{sample}/{sample}.merge.bam")
            ls.append(f"analysis/haplotypecaller/{sample}/{sample}.vcf.gz")
        return ls

    def get_ref(self, config):
        infh = open(config['ref'])
        ref_dic = yaml.safe_load(infh)
        for key, value in ref_dic.items():
            key = f'ref_{key}'
            config[key] = value
        return config

    def get_paths(self, config):
        conda_root = subprocess.check_output('conda info --root', shell=True).decode('utf-8').strip()
        if not "conda_bin_path" in config or not config["conda_bin_path"]:
            config["conda_bin_path"] = os.path.join(conda_root, 'envs', 'kdna', 'bin')
        return config

    def add_config(self, config):
        config["od_smp_s"] = [sample for sample in config["samples"]]
        return config

    def print_config(self, config):
        print("="*100)
        for key, value in sorted(config.items()):
            if not value:
                pass
                #continue
            if key in ['samples']:
                for subkey, subvalue in sorted(value.items()):
                    for subsubkey, subsubvalue in sorted(subvalue.items()):
                        for subsubsubkey, subsubsubvalue in sorted(subsubvalue.items()):
                            print(f"{key} : {subkey} : {subsubkey} : {subsubsubkey} : {subsubsubvalue}")
                #print(f"{key} : {value}")
            else:
                print(f"{key} : {value}")
        print("="*100)
        return config

