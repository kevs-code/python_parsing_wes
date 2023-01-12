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

os.chdir("/home/kevin/mrd_wes/results/2018-05-09/set/bedtools_isec")
directory =("/home/kevin/mrd_wes/results/2018-05-09/set/bedtools_isec")

files=sorted(os.listdir(directory))
asum=bsum=csum=dsum=fsum=gsum=hsum=c=d=0
np.warnings.filterwarnings('ignore')#if np.isnan(np.sum(out_vec)):for MRD-16-028 long scalar
all_tumorAF=[];all_tumorDP=[];all_tumorVD=[];all_rsid=[];all_index=[];all_tumorRD=[]
snps_index=[];snps_tumorAF=[];snps_tumorDP=[];snps_tumorVD=[];snps_normalDP=[];snps_normalAF=[];
snps_rsid=[];snps_normalVD=[];snps_normalRD=[];snps_tumorRD=[]
    #args = cli_interface()
    #files = args.input_filename
    #aft = args.af_filter
    #dpt = args.dp_filter
for filename in files:
    #if filename.endswith(".vcf.gz"):
    if filename.endswith("-ensemble-annotated.vcf.gz"):
        outname=filename
        outname=re.sub(".vcf.gz","-snps.vcf",filename)
        vcf_reader = vcf.Reader(open(filename,'rb'))
        vcf_writer = vcf.Writer(open(outname, 'w'), vcf_reader)
        a=b=f=g=h=0
        for r in vcf_reader:
            for rs in r.samples:
                if rs.sample.endswith('aPC'):
                    a+=1
                    if 'mutect2' in r.INFO['CALLERS']:###individual
                        f+=1
                        if r.INFO['CALLERS'][0]=='mutect2':
                            h+=1
                    else:
                        #more=set()
                        #for ann in r.INFO['ANN']:
                            #get=ann.split('|')[3]
                            #more.add(get)
                        #for gene in more:
                        #r.INFO['ANN'][0].split('|')[3]
                        #print (r.INFO['CALLERS'],r.CHROM,r.POS,r.REF,r.ALT)
                        g+=1
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
                            vcf_writer.write_record(r)
                            snps_index.append(co);snps_rsid.append(r.ID);snps_tumorDP.append(DP1)
                            snps_tumorAF.append(AF1);snps_tumorVD.append(VD1);snps_tumorRD.append(RD1)
                except AttributeError:
                    None

        print(filename, a, b ,f, h, g)#f
        #print(outname,b)
        asum = asum + a
        bsum = bsum + b
        fsum = fsum +f
        gsum = gsum +g
        hsum = hsum +h
print(asum, bsum, fsum, hsum, gsum)


collect = pd.DataFrame({'index': all_index, 'rsid': all_rsid, 'AF': all_tumorAF, 'VD': all_tumorVD,
                    'RD': all_tumorRD, 'DP': all_tumorDP})
snps = pd.DataFrame({'index': snps_index, 'rsid': snps_rsid, 'tumorAF': snps_tumorAF, 'tumorRD': snps_tumorRD,
                     'tumorVD': snps_tumorVD, 'tumorDP': snps_tumorDP, 'normalAF': snps_normalAF, 'normalRD': snps_normalRD,
                     'normalVD': snps_normalVD, 'normalDP': snps_normalDP})
rest = get_rest_dataframe(collect)
collect = collect.set_index('index')
snps = snps.set_index('index')
rest = rest.set_index('index')

feed_dataframe('all',collect)#you could do on the dataframe directly better surely
feed_dataframe('rest',rest)#
feed_dataframe('snps_normal',snps)#
feed_dataframe('snps_tumor',snps)#


####some tables ...
df2 = pd.DataFrame({'normalDP': snps.sort_values(by=['normalDP'])['normalDP'].values,'tumorDP': snps.sort_values(by=['tumorDP'])['tumorDP'].values})
plt.plot(range(350), 'b--')
plt.xlim(0, 350)
plt.ylim(0, 350)
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
ax2=plt.scatter(df2['normalDP'],df2['tumorDP'],c='blue',marker=".",linewidths=0.3)
plt.xlabel("DP (normal)")
plt.ylabel("DP (tumor)")
plt.title('snps DP tumor/ normal samples')

#labels = ['AF']
#df = pd.DataFrame({'AF':tumorAF1})
#ax1=df.hist(stacked=True, bins=10, range=(0,1),grid=False)
#plt.xlabel("snps AF tumor")
#plt.ylabel("Mutations")
#plt.title('snps AF tumor samples (e.g., with rsid)')

ax=df2.boxplot(figsize=(15,10),rot=60,fontsize=8)
ax.grid(False)
ax.set_ylim(0,250)

df6 = pd.DataFrame({'normRD': snps.sort_values(by=['normalRD'])['normalRD'].values,'normVD': snps.sort_values(by=['normalVD'])['normalVD'].values,
                    'tumorRD': snps.sort_values(by=['tumorRD'])['tumorRD'].values,'tumorVD': snps.sort_values(by=['tumorVD'])['tumorVD'].values})
ax6=df6.boxplot(figsize=(15,10),rot=60,fontsize=8)
ax6.grid(False)
ax6.set_ylim(0,200)

ax6=df6.boxplot(figsize=(15,10),rot=60,fontsize=8)
ax6.grid(False)

df5 = pd.DataFrame({'normRD': snps.sort_values(by=['normalRD'])['normalRD'].values,'normDP': snps.sort_values(by=['normalDP'])['normalDP'].values})

plt.plot(range(350), 'b--')
plt.xlim(0, 350)
plt.ylim(0, 350)
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
ax3=plt.scatter(df5['normDP'],df5['normRD'],c='blue',marker=".",linewidths=0.3)
plt.xlabel("DP (normal)")
plt.ylabel("RD (normal)")
plt.title('snps DP / RD normal samples')