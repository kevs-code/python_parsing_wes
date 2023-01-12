#!/usr/bin/env python
import argparse
import vcf
import re
import numpy as np
def cli_interface():
    parser = argparse.ArgumentParser(description='VCF Filter.')
    parser.add_argument('--input', dest='input_filename', nargs='+', type=str, help='input VCF files', required=True)
    parser.add_argument('--af', dest='af_filter', type=float, help='AF Filter >=', required=True)
    parser.add_argument('--dp', dest='dp_filter', type=int, help='DP Filter >=', required=True)
    args = parser.parse_args()
    return args
def main():
    args = cli_interface()
    files = args.input_filename
    aft = args.af_filter
    dpt = args.dp_filter
    for filename in files:
        if filename.endswith(".vcf.gz"):
            budge=filename
            nudge=("-filter%s.vcf" % aft)
            budge=re.sub(".vcf.gz",nudge,filename)
            vcf_reader = vcf.Reader(open(filename,'rb'))
            vcf_writer = vcf.Writer(open(budge, 'w'), vcf_reader)
            a=b=0
            for r in vcf_reader:
                for rs in r.samples:
                    #if rs.sample.endswith('PB')==False:#foraP2
                    if rs.sample.endswith('aPC') or rs.sample.endswith('aP2'):
                        if type(r.genotype(rs.sample)['AD'])==list:
                            VD=r.genotype(rs.sample)['AD'][1];DP=np.sum(r.genotype(rs.sample)['AD']);AF=VD/DP
                        else:
                            VD=r.genotype(rs.sample)['AD'];DP=r.genotype(rs.sample)['DP'];AF=VD/DP
                        if AF>0:
                            a+=1
                        value = int(round(DP * aft))
                        if DP>=dpt and AF>=aft or VD>value:
                            try:
                                vcf_writer.write_record(r)
                            except AttributeError:
                                None
if __name__ == "__main__":
    main()
