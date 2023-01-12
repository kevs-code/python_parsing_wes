#!/usr/bin/env python
import numpy as np
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

os.chdir("/home/kevin/mrd_wes/processed_data/queries/unfiltered_queries/replicates_as_query/staty/unfiltered")
df = pd.read_csv('varscan.csv', sep=',')
df['sompycmd']
named=[]
for s in df['sompycmd'].values:
    t=re.sub(r'^.*-o MRD','MRD',s)
    t=re.sub('varscan','',t)
    named.append(t)
df['name']=named
dataframe=df[['name','type', 'total.truth', 'total.query', 'tp', 'fp', 'fn', 'recall', 'recall2','precision']]

fscores=[]
for k in range(len(dataframe)):
    sensitivity=dataframe['recall'].values[k]
    ppv=dataframe['precision'].values[k]
    fscore=2*((sensitivity*ppv)/(sensitivity+ppv))
    fscores.append(fscore)
dataframe['fscore']=fscores
dataframe = dataframe.round(2)
dataframe

df1 = dataframe.loc[dataframe['type'] == 'records']
df2 = dataframe.loc[dataframe['type'] == 'indels']
df3 = dataframe.loc[dataframe['type'] == 'SNVs']

#make a plot
#plt.scatter(df1['fscore'])
plt.plot(df1['name'],df1['fscore'],'-s',color='purple')
plt.plot(df2['name'],df2['fscore'],'-s',color='blue')
plt.plot(df3['name'],df3['fscore'],'-s',color='red')
#plt.show()

plt.xlabel("Patient")
plt.ylabel("F$_1$ Score")
plt.title("varscan variant calling F$_1$ Score for aliquots")
purple_line = mlines.Line2D([], [], color='purple', marker='s', label='records')
blue_line = mlines.Line2D([], [], color='blue', marker='s', label='indels')
red_line = mlines.Line2D([], [], color='red', marker='s', label='SNVs')
plt.legend(handles=[purple_line,blue_line,red_line])
plt.show()

#edatadf=datadf[:25

df1=df1.set_index('name')
ax1 = pd.concat([df1['total.truth'].rename('truth'),df1['total.query'].rename('query')]
                ,axis=1).plot.bar(title ="varscan vcf records", figsize=(15,10),legend=True,
                                  color=['g', 'yellow'], rot=60, fontsize=8)
ax1.set_ylabel("Records")
plt.show()

#df1 = pd.read_csv('all_transcripts_for_gene_mutations2.csv', sep=',')
#newdf = df1.loc[df1['Feature_ID'].isin(set(df['Feature']))]#22857
#pewf = df1.loc[df1['Feature_ID'].isin(set(df['Feature'])) & df1['Gene_ID'].isin(set(df['Gene']))]#22857
#df1['poscheck']=df1['pos'].values.astype(int)
#df['position']=df['Start_Position'].values.astype(int)
#df['position']=df['Start_Position'].values.astype(int)
#rewf = df1.loc[df1['Feature_ID'].isin(set(df['Feature'])) & df1['Gene_ID'].isin(set(df['Gene']))
       #        & df1['HGVS.c'].isin(set(df['HGVSc']))]
#trewf = df1.loc[df1['Feature_ID'].isin(set(df['Feature'])) & df1['Gene_ID'].isin(set(df['Gene']))
    #            & df1['pos'].isin(set(df['Start_Position']))]
#trewf2 = df1.loc[df1['Feature_ID'].isin(set(df['Feature'])) & df1['Gene_ID'].isin(set(df['Gene']))
 #                & df1['poscheck'].isin(set(df['position']))]
