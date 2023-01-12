#!/usr/bin/env python
import os
import re
import subprocess
import argparse
def cli_interface():
    parser = argparse.ArgumentParser(description='Run som.py comparison on two folders with same filenames')
    parser.add_argument('--truth', dest='truth_dir', type=str, help='truth directory', required=True)
    parser.add_argument('--query', dest='query_dir', type=str, help='query directory', required=True)
    args = parser.parse_args()
    return args
def main():
    args = cli_interface()
    truedir = args.truth_dir
    querydir = args.query_dir
    for truename in os.listdir(truedir):
        if truename.endswith("-annotated.vcf.gz"):
            if '-PB-' not in truename:
                truth = re.sub("-annotated.vcf.gz", "", truename)
                if truename in os.listdir(querydir):
                    print("bcftools isec -e -n=2  %s/%s %s/%s -p %s\n" % (truedir,truename,querydir,truename,truth))
                    nsame = subprocess.Popen("bcftools isec -e -n=2  %s/%s %s/%s -p %s" % (truedir,truename,querydir,truename,truth), shell=True)
                    nsame.communicate()
if __name__ == "__main__":
    main()
