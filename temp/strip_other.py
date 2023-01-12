#!/usr/bin/env python
#import argparse
import vcf
import os
import re
import subprocess
os.chdir("/home/kevin/mrd_wes/results/2018-05-09/set/other_way")
directory = "/home/kevin/mrd_wes/results/2018-05-09/set/other_way"
e=[];
a=b=0
for file in os.listdir(directory):
    if os.path.isdir(file) and file.startswith('MRD'):
        caller=re.sub('$','-annotated.vcf',file)
        file_path=("%s/%s" % (file,caller))
        out_path=re.sub('-annotated','',file_path)
        temp = open(out_path,'w')
        print (file_path, out_path)
        with open(file_path, 'r') as f:
            for line in f:
                if re.search("ANN=", line):
                    simples=re.split("\t",line)
                    if ',' not in simples[15]:
                        sim='\t'.join(simples[11:21])
                        temp.write("%s\n" % sim)
                else:
                    temp.write(line)
        temp.close()
        nsame = subprocess.Popen("gzip %s" % out_path, shell=True)
        nsame.communicate()