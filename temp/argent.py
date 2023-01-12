#!/usr/bin/env python
import numpy as np
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
os.chdir("/home/kevin/mrd_wes/processed_data/queries/unfiltered_queries/replicates_as_query")
df = pd.read_csv('ensemble1.csv', sep=',')
df['sompycmd']
named=[]
for s in df['sompycmd'].values:
    t=re.sub(r'^.*-o MRD','MRD',s)
    t=re.sub('-ensemble','',t)
    named.append(t)
df['name']=named
dataframe=df[['name','type', 'total.truth', 'total.query', 'tp', 'fp', 'fn', 'recall', 'recall2','precision']]

dataframe = dataframe.round(2)
# = pd.DataFrame(dataframe.values, columns=['name','type', 'tot.truth', 'tot.query', 'tp', 'fp', 'fn', 'recall', 'recall2','precision'])
dataframe = dataframe.set_index('name')

#type	total.truth	total.query	tp	fp	fn	unk	ambi	recall	recall_lower	recall_upper	recall2	precision	precision_lower	precision_upper	na	ambiguous	fp.region.size	fp.rate	sompyversion	sompycmd

#df1 = type	total.truth	total.query	tp	fp	fn	recall	recall_lower	recall_upper	recall2	precision	precision_lower	precision_upper	sompycmd


df1 = pd.read_csv('all_transcripts_for_gene_mutations2.csv', sep=',')
newdf = df1.loc[df1['Feature_ID'].isin(set(df['Feature']))]#22857
pewf = df1.loc[df1['Feature_ID'].isin(set(df['Feature'])) & df1['Gene_ID'].isin(set(df['Gene']))]#22857
df1['poscheck']=df1['pos'].values.astype(int)
df['position']=df['Start_Position'].values.astype(int)
df['position']=df['Start_Position'].values.astype(int)
rewf = df1.loc[df1['Feature_ID'].isin(set(df['Feature'])) & df1['Gene_ID'].isin(set(df['Gene'])) & df1['HGVS.c'].isin(set(df['HGVSc']))]
trewf = df1.loc[df1['Feature_ID'].isin(set(df['Feature'])) & df1['Gene_ID'].isin(set(df['Gene'])) & df1['pos'].isin(set(df['Start_Position']))]
trewf2 = df1.loc[df1['Feature_ID'].isin(set(df['Feature'])) & df1['Gene_ID'].isin(set(df['Gene'])) & df1['poscheck'].isin(set(df['position']))]