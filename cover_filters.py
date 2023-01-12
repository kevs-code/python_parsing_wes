#!/usr/bin/env python
import argparse
import re
import subprocess
def cli_interface():
    parser = argparse.ArgumentParser(description='VCF Filter.')
    parser.add_argument('--input', dest='input_filename', nargs='+', type=str, help='input VCF files', required=True)
    args = parser.parse_args()
    return args
def main():
    args = cli_interface()
    files = args.input_filename
    for filename in files:
        if filename.endswith("-annotated.vcf.gz"):
            budge=filename
            budge=re.sub("-annotated.vcf.gz","-annotated-Covered.vcf.gz",filename)
            print("/usr/local/bcbio/anaconda/bin/gatk-launch SelectVariants -V %s -L ~/Agilent_SureselectXT_v6_r2_Covered.bed -O %s\n" % (filename,budge))
	    nsame = subprocess.Popen("/usr/local/bcbio/anaconda/bin/gatk-launch SelectVariants -V %s -L ~/Agilent_SureselectXT_v6_r2_Covered.bed -O %s" % (filename,budge), shell=True)
	    nsame.communicate()
if __name__ == "__main__":
    main()
#/usr/local/bcbio/anaconda/bin/gatk-launch SelectVariants -V MRD-16-028-ensemble-annotated.vcf.gz -L ~/Agilent_SureselectXT_v6_r2_Covered.bed -O whho.vcf.gz
