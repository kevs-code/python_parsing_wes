#!/usr/bin/env python
import argparse
import vcf
import re
import os
import numpy as np
import pandas as pd

os.chdir("/home/kevin/mrd_wes/results/2018-05-09/set/other_way")
directory = "/home/kevin/mrd_wes/results/2018-05-09/set/other_way"
#names=set()
for file in os.listdir(directory):
    if os.path.isdir(file) and file.startswith('MRD'):
        caller=re.sub('$','.vcf.gz',file)
        frame=re.sub('$','-df.csv',file)
        #word=re.sub(r'-\w+$','',file)
        #names.add(word)
        file_path=("%s/%s" % (file,caller))
        vcf_reader = vcf.Reader(open(file_path,'rb'))
        #vcf_writer = vcf.Writer(open(budge, 'w'), vcf_reader)
        a=b=c=0
        results=[]
        for r in vcf_reader:
            a+=1
            for rs in r.samples:
                if rs.sample.endswith('aPC'):# and r.ID.startswith("rs"):
                    if 'strelka2' not in file_path:
                        if type(r.genotype(rs.sample)['AD'])==list:
                            VD=r.genotype(rs.sample)['AD'][1];DP=np.sum(r.genotype(rs.sample)['AD']);AF=VD/DP
                            RD=DP-VD
                        else:
                            VD=r.genotype(rs.sample)['AD'];DP=r.genotype(rs.sample)['DP'];AF=VD/DP;RD=DP-VD
                    else:
                        DP=r.genotype(rs.sample)['DP'];AF=r.genotype(rs.sample)['AF'];VD=AF*DP;RD=DP-VD
                    if AF>0:
                        b+=1
                    if AF>0.1:
                        c+=1
                tmp_record = [r.CHROM, r.POS, r.REF, r.ALT, r.ID, AF, VD, RD, DP]
                results.append(tmp_record)
        print(caller,a,b,c)
        df=pd.DataFrame(results, columns=['chrom', 'pos', 'ref', 'alt', 'rsid', 'af', 'vd', 'rd','dp'])
        df.to_csv(frame)