#!/usr/bin/env python
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
import vcf
import os, sys
import re
os.chdir("/home/kevin/mrd_wes/results/2018-05-09/make_a_table")
directory =("/home/kevin/mrd_wes/results/2018-05-09/make_a_table")
e=[];
a=b=0
for filename in os.listdir(directory):
    if filename.endswith("-ensemble-annotated.txt"):
        with open(filename, 'r') as f:
            for line in f:
                if re.search("ANN=", line):
                    starthead="CHROM;POS;ID;REF;ALT;QUAL;FILTER"
                    simples=re.split("\t",line)
                    head=';'.join(simples[0:7])
                    if simples[2].startswith('rs'):
                        print(starthead)
                        print(head)
                        #sub checkformat
                        if simples[8]=='GT:AD:ADJAF:AF:ALD:BIAS:DP:HIAF:MQ:NM:ODDRATIO:PMEAN:PSTD:QSTD:QUAL:RD:SBF:SN:VD':
                            vardicth=re.split(":",simples[8])
                            vardict1=re.split(":",simples[9])
                            vardict2=re.split(":",simples[10])
                            ADi=1;DPi=6
                            DPsum=np.sum(vardict1[6],vardict2[6])
                        if simples[8]=='GT:AD:DP:DP4:FREQ:RD':
                            varscanh=re.split(":",simples[8])
                            varscan1=re.split(":",simples[9])
                            varscan2=re.split(":",simples[10])
                        if re.search("GT:AD:AF:DP:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS", simples[8]):
                            mutecth=re.split(":",simples[8])
                            mutect1=re.split(":",simples[9])
                            mutect2=re.split(":",simples[10])
#######works waste of time!!!!
#######################################################################################
                        if type(r.genotype(rs.sample)['AD'])==list:
                            VD=r.genotype(rs.sample)['AD'][1];DP=np.sum(r.genotype(rs.sample)['AD']);AF=VD/DP
                            ND=r.genotype(rs.sample)['AD'][0]
                            DPY=r.genotype(rs.sample)['DP']
                        else:
                            VD=r.genotype(rs.sample)['AD'];DP=r.genotype(rs.sample)['DP'];AF=VD/DP
                            ND=DP-VD
                            DPY=r.genotype(rs.sample)['DP']
                        if rs.sample.endswith('PB'):
                            normalDP.append(DP)
                            normalDPY.append(DPY)
                        if rs.sample.endswith('aPC'):
                            tumorDP.append(DP)
                            tumorDPY.append(DPY)
#######################################################################################
                        #if simples[8] not in e:#gives info headers and field length
                        #    e.append(simples[8])


print "GT:AD:AF:DP:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB $d\n";#MUTECT2
print "GT:AD:ADJAF:AF:ALD:BIAS:DP:HIAF:MQ:NM:ODDRATIO:PMEAN:PSTD:QSTD:QUAL:RD:SBF:SN:VD $e\n";#VARDICT
print "GT:AD:DP:DP4:FREQ:RD $f\n";#VARSCAN
print "GT:AD:AF:DP:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:PGT:PID:SA_MAP_AF:SA_POST_PRO $g\n";#MUTECT2

            first = re.search(r'(^.*)(CALLERS=.*)(ANN=.*\s+)(GT.*)$', line)
            base = re.sub("\s+", ",", first.group(1))#start
            annot = re.sub("^ANN=","", first.group(3))#ANNls
##############################################################################
            format, allelle1, allelle2 = re.split("\s+", first.group(4))
            mow=format
            mow = re.sub(":PGT:PID:",":",mow)
            joe=mow
            toe=mow#maybe keep format as style dependent
            joe = re.sub(":", ":t1_", joe)
            joe = re.sub("^GT:","t1_GT:",joe)
            toe = re.sub("^GT:","t2_GT:",toe)
            toe = re.sub(":",":t2_",toe)
            toe = re.sub(":t2_SA_MAP_AF:t2_SA_POST_PROB","",toe)
            #####################
            if re.search(r'GT:AD:AF:DP:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:PGT:PID:SA_MAP_AF:SA_POST_PROB',format):
                allelle2 = re.sub(":.?:.*?$","",allelle2)
                allelle1 = re.sub(":\d+\|\d+:.*?:",":",allelle1)
            #####################
            seq = (allelle1, allelle2)
            head = (joe, toe)
            s = ":"
            tupee = s.join(seq)
            oukey = s.join(head)#format header
            infohead=re.sub("=.*?;",";",first.group(2))#info header
            infohead=re.sub(";$","",infohead)
            if infohead not in e:#gives info headers and field length
                if 'DECOMP' not in infohead:#removes decomposed
                    if not re.search(";RAW_MQ;ReadPosRankSum;TLOD",infohead):
                        e.append(infohead)# not decompose = 4 now 3
            if oukey not in t:
                t.append(oukey)
                new=re.split(":",oukey)
                q.append(new)
            if format not in g:
                simples=re.split(":",format)
                print(len(simples),simples)
                g.append(format)#4 formats logged in g length vardict=80, mutect=70 +62, varscan=20
#######################still messy really
            a.append(base)#start bit content "," sep
            b.append(first.group(2))#info content ";" sep
            d.append(tupee)#format content ":" sep
            setstart=re.sub(";$","",base)
            setstart=re.sub(";",":",setstart)
            setstart=re.split(",",setstart);####start with or without AF from info
            setinfo=re.sub(";$","",first.group(2))
            parse = open('parse.csv','a')
            if re.search(r'CALLERS=vardict.*;DECOMPOSED;LEN=\d+',setinfo):
                setinfo=re.sub("$",";"*10,setinfo)
                parse.write("%s\n" % setinfo)
            if re.search(r'CALLERS=mutect2.*;DECOMPOSED;LEN=\d+',setinfo):
                setinfo=re.sub("$",";"*12,setinfo)
                parse.write("%s\n" % setinfo)
            if re.search(r'RAW_MQ=[+-]?((\d+(\.\d*)?)|\.\d+)([eE][+-]?[0-9]+)?;ReadPosRankSum',setinfo):
                setinfo=re.sub("ReadPosRankSum",";;ReadPosRankSum",setinfo)
                parse.write("%s\n" % setinfo)
            if re.search(r'ReadPosRankSum=[+-]?((\d+(\.\d*)?)|\.\d+)([eE][+-]?[0-9]+)?;TLOD=',setinfo):
                setinfo=re.sub("TLOD",";TLOD",setinfo)
            setinfo=re.split(";",setinfo);####info 4 styles 3 formats
            setformat=re.split(":",tupee);####format 4 styles 3 formats
###########################################################
            #parse = open('parse.csv','a')
            annot = re.sub("LOF=\(.*\)\s+$","",annot)#LOF for now removed save pfft trouble below
                #pfft = re.findall(r';(LOF=\(.*\))\s+$',annot)
                #parse.write("%s\n" % pfft)
                #pfft = re.search(r';(LOF=\(.*\))\s+$',annot)
                #face = re.sub("\|", ":", pfft.group(1))
                #annot = re.sub("LOF=\(.*$",face,annot)
                #pfft   #;LOF=(C2|ENSG00000166278|17|0.35),(CFB|ENSG00000243649|15|0.07),(CFB|ENSG00000244255|2|1.00)
            setannot=re.split(",",annot)

            for l in range(len(setannot)):
                setannot[l] = re.sub("\|", ",", setannot[l])
                c.append(setannot[l])
                orange=setstart+setinfo
                h.append(orange)
                m.append(setformat)
            i += 1
parse.close()