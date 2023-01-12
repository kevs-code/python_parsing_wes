#!/usr/bin/env python
import os
import re
import subprocess
import argparse
def cli_interface():
    parser = argparse.ArgumentParser(description='File for vcf2maf')
    parser.add_argument('--input', dest='input_filename', nargs='+', type=str, help='input VCF files uncompressed', required=True)
    parser.add_argument('--mafdir', dest='mafdirect', type=str, help='input vcf2maf.pl directory path default is ~/mskcc-vcf2maf-747a1bb/', required=False)
    args = parser.parse_args()
    return args
def main():
    args = cli_interface()
    files = args.input_filename
    try:
        mafdir = args.mafdirect
    except AttributeError:
        None
    if mafdir == None:
        mafdir = ("~/mskcc-vcf2maf-747a1bb")
    print (mafdir)
    for filename in files:
        if filename.endswith("vcf"):
            outname=re.sub("vcf","vep.maf",filename)
            with open(filename, 'r') as f:
                for line in f:
                    if re.search("normal_sample", line):
                        normal = re.search(r'sample=(.*)\s*$', line)
                    if re.search("tumor_sample", line):
                        tumor = re.search(r'sample=(.*)\s*$', line)
            print (normal.group(1),tumor.group(1))
            print ("perl %s/vcf2maf.pl --input-vcf %s --output-maf %s --tumor-id %s --normal-id %s" % (mafdir, filename, outname, tumor.group(1), normal.group(1)))
            nsame = subprocess.Popen("perl %s/vcf2maf.pl --input-vcf %s --output-maf %s --tumor-id %s --normal-id %s" % (mafdir, filename, outname, tumor.group(1), normal.group(1)), shell=True)
            nsame.communicate()
if __name__ == "__main__":
    main()
