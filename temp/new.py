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
        VD=ADY[1];DP=np.sum(ADY);AF=VD/DP;ND=DP-VD#or ADY[0]
        if len(ADY)>2:
            print ("error at index %s GT:AD>2" % c)
    else:
        VD=ADY;DP=DPY;AF=VD/DP;ND=DP-VD
    fields = [AF,VD,DP,ND]
    return fields

def iquart(name, series, sumseries):
    las=sumseries-1
    q1=series[int(round(sumseries/4))];q2=series[int(round(sumseries/2))];q3=series[int(round((sumseries/4)*3))];
    u=np.mean(series);iqr=q3-q1;out=iqr*1.5;upper=q3+out;lower=q1-out;lower=max(0,lower);
    first=series[0];last=series[las]
    if (name=='tumorAF'):
        print("%s\nQ1 %.2f Q2 %.2f Q3 %.2f IQR %.2f +/- outlier %.2f upper %.2f lower %.2f min %.2f max %.2f mean %.2f"
              % (name,q1,q2,q3,iqr,out,upper,lower,first,last,u))
        print("name,q1,q2,q3,iqr,out,upper,lower,min,max,mean")
        print("%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f"
              % (name,q1,q2,q3,iqr,out,upper,lower,first,last,u))
    else:
        print("%s\nQ1 %.0f Q2 %.0f Q3 %.0f IQR %.0f +/- outlier %.1f upper %.1f lower %.1f min %.1f max %.1f mean %.1f"
              % (name,q1,q2,q3,iqr,out,upper,lower,first,last,u))
        print("name,q1,q2,q3,iqr,out,upper,lower,min,max,mean")
        print("%s,%.0f,%.0f,%.0f,%.0f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f"
              % (name,q1,q2,q3,iqr,out,upper,lower,first,last,u))

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
    #parser.add_argument('--af', dest='af_filter', type=float, help='AF Filter >=', required=True)
    #parser.add_argument('--dp', dest='dp_filter', type=int, help='DP Filter >=', required=True)
    args = parser.parse_args()
    return args

def main():
    asum=bsum=csum=dsum=c=d=0
    all_tumorAF=[];all_tumorDP=[];all_tumorVD=[];all_rsid=[];all_index=[];snps_index=[]
    snps_tumorAF=[];snps_tumorDP=[];snps_tumorVD=[];snps_normalDP=[];snps_rsid=[]
    args = cli_interface()
    files = args.input_filename
    #aft = args.af_filter
    #dpt = args.dp_filter
    for filename in files:
        if filename.endswith(".vcf.gz"):
            #outname=filename
            #outname=re.sub(".vcf.gz","-filter.vcf",filename)
            vcf_reader = vcf.Reader(open(filename,'rb'))
            #vcf_writer = vcf.Writer(open(outname, 'w'), vcf_reader)
            a=b=0
            for r in vcf_reader:
                for rs in r.samples:
                    if rs.sample.endswith('aPC'):
                        a+=1
                        fields=get_AF_fields(r.genotype(rs.sample)['AD'],r.genotype(rs.sample)['DP'])
                        AF=fields[0];VD=fields[1];DP=fields[2]
                        if AF>0 and DP > VD:
                            all_tumorAF.append(AF);all_tumorDP.append(DP);all_tumorVD.append(VD)
                            all_rsid.append(r.ID);all_index.append(c)
                            #value = int(round(DP * aft))
                            c+=1
                try:
                    if r.ID.startswith("rs"):
                        snps=get_AF_fields(r.genotype(rs.sample)['AD'],r.genotype(rs.sample)['DP'])
                        AF1=snps[0];VD1=snps[1];DP1=snps[2]
                        if AF1>0 and DP1 > VD1:
                            if rs.sample.endswith('PB'):
                                snps_normalDP.append(DP1)
                            if rs.sample.endswith('aPC'):
                                d+=1
                                snps_index.append(c);snps_rsid.append(r.ID);snps_tumorDP.append(DP1)
                                snps_tumorAF.append(AF1);snps_tumorVD.append(VD1)
                except AttributeError:
                    None

            print(filename, a)#filename a
            #print(outname,b)
            asum = asum + a
            bsum = bsum + b
            print(asum, bsum)

    all = pd.DataFrame({'index': all_index, 'rsid': all_rsid, 'AF': all_tumorAF, 'VD': all_tumorVD, 'DP': all_tumorDP})
    snps = pd.DataFrame({'index': snps_index, 'rsid': snps_rsid, 'AF': snps_tumorAF,
                         'VD': snps_tumorVD, 'tumorDP': snps_tumorDP, 'normalDP': snps_normalDP})
    rest = get_rest_dataframe(all)
    alldf = all.set_index('index')
    restdf = rest.set_index('index')

    #all = pd.DataFrame({'index': all_index, 'rsid': all_rsid, 'AF': all_tumorAF, 'VD': all_tumorVD, 'DP': all_tumorDP})
    #meh = sorted(all_tumorAF)
    #peh = sorted(all_tumorDP)
#x1=sorted(normalDP)
#y1=sorted(tumorDP)
#tumorAF1=sorted(tumorAF)

#iqrt('normalDP',x1,c)
#iqrt('tumorDP',y1,c)
#iqrt('tumorAF',tumorAF1,c)


df2 = pd.DataFrame({'x1': x1,'y1': y1})
plt.plot(range(350), 'b--')
plt.xlim(0, 350)
plt.ylim(0, 350)
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
ax2=plt.scatter(df2['x1'],df2['y1'],c='blue',marker=".",linewidths=0.3)
plt.xlabel("DP (normal)")
plt.ylabel("DP (tumor)")
plt.title('snps DP tumor/ normal samples')
labels = ['AF']
df = pd.DataFrame({'AF':tumorAF1})
ax1=df.hist(stacked=True, bins=10, range=(0,1),grid=False)
plt.xlabel("snps AF tumor")
plt.ylabel("Mutations")
plt.title('snps AF tumor samples (e.g., with rsid)')

#data_frame.to_csv(args.output_filename)
                if DP>=dpt and AF>=aft or VD>value:
                            try:
                                vcf_writer.write_record(r)
                                b+=1
                            except AttributeError:
                                None
            print (filename,'records VD>0 =',a,'filtered records =',b)
if __name__ == "__main__":
    main()
