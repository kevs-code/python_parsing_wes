#!/usr/bin/env python
import numpy as np
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

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

plt.xlabel("Donor Query")
plt.ylabel("F$_1$ Score")
plt.title("ensemble variant calling F$_1$ Score for aliquots aPC truth aP2 query")
purple_line = mlines.Line2D([], [], color='purple', marker='s', label='records')
blue_line = mlines.Line2D([], [], color='blue', marker='s', label='indels')
red_line = mlines.Line2D([], [], color='red', marker='s', label='SNVs')
plt.legend(handles=[purple_line,blue_line,red_line])
plt.show()
data=datframe['total.truth','total.query']

plt.bar(data)
plt.show()

fig1= df1[['total.truth','total_query']].plot.bar\
    (x='meh',y='records', title ="Top e", figsize=(15,10),legend=False, fontsize=12)




Which is pretty much the same as

plt.plot(dates, values, 'or-')
plt.show()