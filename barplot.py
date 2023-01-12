#!/usr/bin/env python
import numpy as np
import os
import pandas as pd
os.chdir("/home/kevin/PycharmProjects/begins/continues")
df = pd.read_csv('23032018_summary_mutations.csv', sep=',')
#,chrom,pos,ref,alt,dp,mq,"mmq,",num_samples,samples,Gene_Name
#donor=df[['chrom', 'pos', 'ref', 'alt']]
Gene_mutation=df[['Gene_Name', 'chrom', 'pos', 'ref', 'alt','num_samples']]
grouped = Gene_mutation.groupby('Gene_Name')
test1=[];#test2=[];
for name, group in grouped:
    lax=group['num_samples'].agg([np.sum, np.mean, np.std])
    seq1=(name,len(group))
    #seq2=(name,lax[0])
    #print (seq2)#barplot gene bt mutations and sum of frequency across patients
    test1.append(seq1)
head1 =['Gene','mutations']
table1 = pd.DataFrame(test1, columns=head1)
plot1=table1.sort_values(by=['mutations'],ascending=False)
fig1= plot1[['Gene','mutations']][:50].plot.bar\
    (x='Gene',y='mutations', title ="Top 50 Number of Mutations per gene", figsize=(15,10),legend=False, fontsize=12)
