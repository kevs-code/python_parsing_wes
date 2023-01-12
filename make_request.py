#!/usr/bin/env python
import argparse
import vcf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os

def get_AF_fields(ADY,DPY):
    if type(ADY)==list:
        VD=ADY[1];DP=np.sum(ADY);AF=VD/DP;RD=DP-VD#or ADY[0]
        if len(ADY)>2:
            print ("error at index %s GT:AD>2" % c)
    else:
        VD=ADY;DP=DPY;AF=VD/DP;RD=DP-VD
    fields = [AF,VD,DP,RD]
    return fields

def iquart(name, series, sumseries, kind):
    las=sumseries-1
    q1=series[int(round(sumseries/4))];q2=series[int(round(sumseries/2))];q3=series[int(round((sumseries/4)*3))];
    u=np.mean(series);iqr=q3-q1;out=iqr*1.5;upper=q3+out;lower=q1-out;lower=max(0,lower);
    first=series[0];last=series[las]
    if 'AF' in name:
        print("kind,name,q1,q2,q3,iqr,out,upper,lower,min,max,mean")
        print("%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f"
              % (kind,name,q1,q2,q3,iqr,out,upper,lower,first,last,u))
    else:
        #print("kind,name,q1,q2,q3,iqr,out,upper,lower,min,max,mean")
        print("%s,%s,%.0f,%.0f,%.0f,%.0f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f"
              % (kind,name,q1,q2,q3,iqr,out,upper,lower,first,last,u))

#def instead_of_iquart(meh):

def feed_dataframe(kind,df):
    count=len(df);set1=['AF','VD','DP','RD'];set2=['normalAF','normalVD','normalDP','normalRD']
    set3=['tumorAF','tumorVD','tumorDP','tumorRD'];set4=set1
    if kind == 'snps_normal':
        set4=set2
    if kind == 'snps_tumor':
        set4=set3
    AFcol=df.sort_values(by=[set4[0]])
    VDcol=df.sort_values(by=[set4[1]])
    DPcol=df.sort_values(by=[set4[2]])
    RDcol=df.sort_values(by=[set4[3]])
    iquart(set4[0],AFcol[set4[0]].values,count,kind)
    iquart(set4[1],VDcol[set4[1]].values,count,kind)
    iquart(set4[2],DPcol[set4[2]].values,count,kind)
    iquart(set4[3],RDcol[set4[3]].values,count,kind)
def get_rest_dataframe(rest1):
    sett=[];t=0
    for s in rest1['rsid']:
        if s==None:
            sett.append(t)
        t+=1
    rest_index = pd.DataFrame({'index': sett})
    rest = rest1.loc[rest1['index'].isin(set(rest_index['index']))]
    return rest


def cli_interface():
    parser = argparse.ArgumentParser(description='VCF Filter.')
    parser.add_argument('--input', dest='input_filename', nargs='+', type=str, help='input VCF files', required=True)
    parser.add_argument('--af', dest='af_filter', type=float, help='AF Filter >=', required=True)
    parser.add_argument('--dp', dest='dp_filter', type=int, help='DP Filter >=', required=True)
    args = parser.parse_args()
    return args

os.chdir("/home/kevin/mrd_wes/results/2018-05-09/set")#/bedtools_isec")
directory =("/home/kevin/mrd_wes/results/2018-05-09/set")#/bedtools_isec")

files=sorted(os.listdir(directory))
asum=bsum=csum=dsum=fsum=c=d=0
np.warnings.filterwarnings('ignore')#if np.isnan(np.sum(out_vec)):for MRD-16-028 long scalar
all_tumorAF=[];all_tumorDP=[];all_tumorVD=[];all_rsid=[];all_index=[];all_tumorRD=[]
snps_index=[];snps_tumorAF=[];snps_tumorDP=[];snps_tumorVD=[];snps_normalDP=[];snps_normalAF=[];
snps_rsid=[];snps_normalVD=[];snps_normalRD=[];snps_tumorRD=[];snps_rsample=[];args=[]
snps_CHR=[];snps_POS=[];snps_REF=[];snps_ALT=[]
    #args = cli_interface()
    #files = args.input_filename
    #aft = args.af_filter
    #dpt = args.dp_filter
for filename in files:
    #if filename.endswith(".vcf.gz"):
    if filename.endswith("-ensemble-annotated.vcf.gz"):
        checkname=re.sub("-ensemble-annotated.vcf.gz","",filename)
        args.append(checkname)
        #outname=filename
        #outname=re.sub(".vcf.gz","-snps.vcf",filename)
        vcf_reader = vcf.Reader(open(filename,'rb'))
        #vcf_writer = vcf.Writer(open(outname, 'w'), vcf_reader)
        a=b=f=0
        for r in vcf_reader:
            for rs in r.samples:
                if rs.sample.endswith('aPC'):
                    a+=1
                    #print (r.INFO['CALLERS'],r.CHROM,r.POS,r.REF,r.ALT)
                    fields=get_AF_fields(r.genotype(rs.sample)['AD'],r.genotype(rs.sample)['DP'])
                    AF=fields[0];VD=fields[1];DP=fields[2];RD=fields[3]
                    #if AF>0 and DP > VD:not here
                    all_tumorAF.append(AF);all_tumorDP.append(DP);all_tumorVD.append(VD)
                    all_rsid.append(r.ID);all_index.append(c);all_tumorRD.append(RD)
                    #value = int(round(DP * aft))
                    c+=1

                try:
                    if r.ID.startswith("rs"):
                        snps=get_AF_fields(r.genotype(rs.sample)['AD'],r.genotype(rs.sample)['DP'])
                        AF1=snps[0];VD1=snps[1];DP1=snps[2];RD1=snps[3]
                        if rs.sample.endswith('PB'):
                            snps_normalDP.append(DP1);snps_normalRD.append(RD1)
                            snps_normalAF.append(AF1);snps_normalVD.append(VD1)
                        if rs.sample.endswith('aPC'):
                            b+=1;d+=1
                            co=c-1
                            #e = (co,r.CHROM,r.POS,r.REF,r.ALT)
                            #seq=','.join(map(str,e))
                            snps_index.append(co);snps_CHR.append(r.CHROM);
                            snps_POS.append(r.POS);snps_REF.append(r.REF);snps_ALT.append(r.ALT)
                            snps_rsid.append(r.ID);snps_tumorDP.append(DP1)
                            snps_tumorAF.append(AF1);snps_tumorVD.append(VD1);snps_tumorRD.append(RD1)
                            snps_rsample.append(rs.sample);
                except AttributeError:
                    None

        print(filename, a, b)
        asum = asum + a
        bsum = bsum + b
print(asum, bsum)
snps_k=[]
for q in range(len(snps_normalDP)):
    v=snps_tumorDP[q]-snps_normalDP[q]
    snps_k.append(v)

collect = pd.DataFrame({'index': all_index, 'rsid': all_rsid, 'AF': all_tumorAF, 'VD': all_tumorVD,
                    'RD': all_tumorRD, 'DP': all_tumorDP})

snps = pd.DataFrame({'index': snps_index, 'rsample': snps_rsample, 'chr': snps_CHR, 'pos': snps_POS, 'ref': snps_REF,
                     'alt': snps_ALT, 'rsid': snps_rsid, 'tumorAF': snps_tumorAF, 'tumorRD': snps_tumorRD,
                     'tumorVD': snps_tumorVD, 'tumorDP': snps_tumorDP, 'normalAF': snps_normalAF,
                     'normalRD': snps_normalRD,'normalVD': snps_normalVD, 'normalDP': snps_normalDP, 'k1-k2': snps_k})
data = snps[['index', 'rsample', 'chr', 'pos', 'ref','alt','rsid','tumorAF','normalAF','tumorRD','normalRD',
             'tumorVD','normalVD','tumorDP','normalDP','k1-k2']]

den = {}
for donor in args:
    name = re.sub('$','-aPC',donor)
    den[donor] = data.loc[data['rsample'] == name]
    outnamed=("%s-snps.csv" % donor)
    data_frame = den[donor]
    data_frame.to_csv(outnamed)
#den["MRD-5105"]['pos']
#or print ...
#d = {name: pd.DataFrame() for donor in args}


rest = get_rest_dataframe(collect)
collect = collect.set_index('index')
snps = snps.set_index('index')
rest = rest.set_index('index')

feed_dataframe('all',collect)#you could do on the dataframe directly better surely
feed_dataframe('rest',rest)#
feed_dataframe('snps_normal',snps)#
feed_dataframe('snps_tumor',snps)#

#test=data.set_index(['rsample', 'chr','pos'])
#pewf = df1.loc[df1['Feature_ID'].isin(set(df['Feature'])) & df1['Gene_ID'].isin(set(df['Gene']))]
#df.loc[df['column_name'] == some_value]
#no=0
#for donor in args:
    #string = 'donor'
    #df.loc[df['column_name'] == some_value]
    #frame[no]=dloc[]
    #frame[0]=data.loc[data['rsample'] == donor]
    #print(donor)
    #string +=str(no)
    #string = data.loc[data['rsample'] == donor]
    #no+=1


#print (snps_tumorDP[q], ' - ', snps_normalDP[q], ' = ',v)
#e = (co,r.CHROM,r.POS,r.REF,r.ALT)
#seq=','.join(map(str,e))
#woop=snps_index[13625].split(',')
#','.join(str(i) for i in e)
#snps_index.append(seq)
#suitable for index when iterated together below instead of set of list in snps do other way for speed
#snps_CHR.append(r.CHROM);
#snps_POS.append(r.POS);snps_REF.append(r.REF);snps_ALT.append(r.ALT)
