#!/usr/bin/env python
import numpy as np
import os
import pandas as pd

os.chdir("/home/kevin/PycharmProjects/begins/continues")
#df = pd.read_csv('20032018_summary_mutations.csv', sep=',')
df = pd.read_csv('high_slicer.csv', sep=',')
Gene_mutation=df[['Gene_Name', 'chrom', 'pos', 'ref', 'alt','samples','num_samples']]
grouped = Gene_mutation.groupby('Gene_Name')
arges = ['4866aPC', 'M-16-031aPC', 'M-17-004aPC', 'M-17-005aPC', 'M-17-021-aPC', 'M-17-024aPC', 'M-17-031aPC', 'M-17-046aPC', 'M-17-052aPC', 'R-17-067-aPC', 'M-17-069aPC','M-17-070aPC', 'M-17-081aPC', '17-108aPC']
test1=[];

def get_sampled_by_sample(ab):
    args1=['0','0','0','0','0','0','0','0','0','0','0','0','0','0']
    i=0
    for so in arges:
        for ro in ab:
            if (so==ro):
                args1[i]=+1
        i+=1
    return args1
look=[]
for name, group in grouped:
    print("hello"+name)
    my_list=[]
    for r in range(len(group['samples'].values)):
        ab=group['samples'].values[r].split(';')
        sampled = get_sampled_by_sample(ab)
        sampled = list(map(int, sampled))
        print(sampled)
        my_list.append(sampled)
    doof=pd.DataFrame(my_list,columns=arges)
    act=doof.sum().values.tolist()
    meep=[name]+act
    look.append(meep)

head2=['gene']+arges

doof2=pd.DataFrame(look,columns=head2)
loft = doof2.set_index('gene')
loft['total'] = loft.sum(axis=1)
loft['donors'] = (loft != 0).astype(int).sum(axis=1)
ok=loft.sort_values(by=['total','donors'],ascending=False)
ok.to_csv('boxplot_again.csv',sep=',')

#(lol == 0).astype(int).sum(axis=1)#

cols_to_use = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
df2 = pd.read_csv('boxplot_again.csv', usecols= cols_to_use, sep=',')
lol=df2[:50]
#(lol == 0).astype(int).sum(axis=1)#
lol = lol.set_index('gene')
ax=lol.T.boxplot(figsize=(15,10),rot=60,fontsize=8)
ax.grid(False)
#ax.right_ax(False)
ax.set_xlabel("Gene")
ax.set_ylabel("Mutation count per donor")
ax.set_title("Top 50 Mutated Genes with High Functional Impact")
#import seaborn as sns
#sns.boxplot(data=lol)