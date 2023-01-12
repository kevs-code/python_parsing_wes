#!/usr/bin/env python
import numpy as np
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

os.chdir("/home/kevin/mrd_wes/processed_data/queries/unfiltered_queries/replicates_as_query/staty/filtered")
df = pd.read_csv('ensemble.csv', sep=',')
df['sompycmd']
named=[]
for s in df['sompycmd'].values:
    t=re.sub(r'^.*-o MRD','MRD',s)
    t=re.sub('ensemble-annotated-filter','af=',t)
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
#df1 = dataframe.loc[dataframe['name'].str.endswith('0.1')]
#df2 = dataframe.loc[dataframe['name'].str.endswith('0.2')]
#df3 = dataframe.loc[dataframe['name'].str.endswith('0.3')]
#df4 = dataframe.loc[dataframe['name'].str.endswith('0.4')]
#df5 = dataframe.loc[dataframe['name'].str.endswith('0.5')]
#plt.scatter(df1['recall'].rename('sensitivity'),df1['precision'].rename('ppv'))


#scatter
#precision=ppv
#recall=sensitivity
df1 = dataframe.loc[dataframe['type'] == 'records']
df2 = dataframe.loc[dataframe['type'] == 'indels']
df3 = dataframe.loc[dataframe['type'] == 'SNVs']

dfr1 = df1.loc[df1['name'].str.endswith('0.1')]
dfr2 = df1.loc[df1['name'].str.endswith('0.2')]
dfr3 = df1.loc[df1['name'].str.endswith('0.3')]
dfr4 = df1.loc[df1['name'].str.endswith('0.4')]
dfr5 = df1.loc[df1['name'].str.endswith('0.5')]

dfi1 = df2.loc[df2['name'].str.endswith('0.1')]
dfi2 = df2.loc[df2['name'].str.endswith('0.2')]
dfi3 = df2.loc[df2['name'].str.endswith('0.3')]
dfi4 = df2.loc[df2['name'].str.endswith('0.4')]
dfi5 = df2.loc[df2['name'].str.endswith('0.5')]

dfs1 = df3.loc[df3['name'].str.endswith('0.1')]
dfs2 = df3.loc[df3['name'].str.endswith('0.2')]
dfs3 = df3.loc[df3['name'].str.endswith('0.3')]
dfs4 = df3.loc[df3['name'].str.endswith('0.4')]
dfs5 = df3.loc[df3['name'].str.endswith('0.5')]

mp1=plt.scatter(dfr1['recall'],dfr1['precision'],color='purple',marker='o')
mp2=plt.scatter(dfi1['recall'],dfi1['precision'],color='blue',marker='s')
mp3=plt.scatter(dfs1['recall'],dfs1['precision'],color='red',marker='D')

plt.xlabel("sensitivity")
plt.ylabel("ppv")
plt.title("ensemble AF=0.1")
plt.legend([mp1,mp2,mp3],["records","indels","SNV"])
#purple_line = mlines.Line2D([], [], color='purple', marker='s', label='records')
#blue_line = mlines.Line2D([], [], color='blue', marker='s', label='indels')
#red_line = mlines.Line2D([], [], color='red', marker='s', label='SNVs')
#plt.legend(handles=[purple_line,blue_line,red_line])
plt.show()

dfr1=dfr1.set_index('name')
ax1 = pd.concat([dfr1['total.truth'].rename('truth'),dfr1['total.query'].rename('query')]
                ,axis=1).plot.bar(title ="AF 0.1 ensemble vcf records", figsize=(15,10),legend=True,
                                  color=['g', 'yellow'], rot=60, fontsize=8)
ax1.set_ylabel("Records")
plt.show()

mp1=plt.scatter(dfr2['recall'],dfr2['precision'],color='purple',marker='o')
mp2=plt.scatter(dfi2['recall'],dfi2['precision'],color='blue',marker='s')
mp3=plt.scatter(dfs2['recall'],dfs2['precision'],color='red',marker='D')

plt.xlabel("sensitivity")
plt.ylabel("ppv")
plt.title("ensemble AF=0.2")
plt.legend([mp1,mp2,mp3],["records","indels","SNV"])
plt.show()

dfr2=dfr2.set_index('name')
ax1 = pd.concat([dfr2['total.truth'].rename('truth'),dfr2['total.query'].rename('query')]
                ,axis=1).plot.bar(title ="AF 0.2 ensemble vcf records", figsize=(15,10),legend=True,
                                  color=['g', 'yellow'], rot=60, fontsize=8)
ax1.set_ylabel("Records")
plt.show()

mp1=plt.scatter(dfr3['recall'],dfr3['precision'],color='purple',marker='o')
mp2=plt.scatter(dfi3['recall'],dfi3['precision'],color='blue',marker='s')
mp3=plt.scatter(dfs3['recall'],dfs3['precision'],color='red',marker='D')

plt.xlabel("sensitivity")
plt.ylabel("ppv")
plt.title("ensemble AF=0.3")
plt.legend([mp1,mp2,mp3],["records","indels","SNV"])
plt.show()

dfr3=dfr3.set_index('name')
ax1 = pd.concat([dfr3['total.truth'].rename('truth'),dfr3['total.query'].rename('query')]
                ,axis=1).plot.bar(title ="AF 0.3 ensemble vcf records", figsize=(15,10),legend=True,
                                  color=['g', 'yellow'], rot=60, fontsize=8)
ax1.set_ylabel("Records")
plt.show()

mp1=plt.scatter(dfr4['recall'],dfr4['precision'],color='purple',marker='o')
mp2=plt.scatter(dfi4['recall'],dfi4['precision'],color='blue',marker='s')
mp3=plt.scatter(dfs4['recall'],dfs4['precision'],color='red',marker='D')

plt.xlabel("sensitivity")
plt.ylabel("ppv")
plt.title("ensemble AF=0.4")
plt.legend([mp1,mp2,mp3],["records","indels","SNV"])
plt.show()

dfr4=dfr4.set_index('name')
ax1 = pd.concat([dfr4['total.truth'].rename('truth'),dfr4['total.query'].rename('query')]
                ,axis=1).plot.bar(title ="AF 0.4 ensemble vcf records", figsize=(15,10),legend=True,
                                  color=['g', 'yellow'], rot=60, fontsize=8)
ax1.set_ylabel("Records")
plt.show()

mp1=plt.scatter(dfr5['recall'],dfr5['precision'],color='purple',marker='o')
mp2=plt.scatter(dfi5['recall'],dfi5['precision'],color='blue',marker='s')
mp3=plt.scatter(dfs5['recall'],dfs5['precision'],color='red',marker='D')


plt.xlabel("sensitivity")
plt.ylabel("ppv")
plt.title("ensemble AF=0.5")
plt.legend([mp1,mp2,mp3],["records","indels","SNV"])
plt.show()

dfr5=dfr5.set_index('name')
ax1 = pd.concat([dfr5['total.truth'].rename('truth'),dfr5['total.query'].rename('query')]
                ,axis=1).plot.bar(title ="AF 0.5 ensemble vcf records", figsize=(15,10),legend=True,
                                  color=['g', 'yellow'], rot=60, fontsize=8)
ax1.set_ylabel("Records")
plt.show()
###################################################################################
df1=df1.set_index('name')
ax1 = pd.concat([df1['total.truth'].rename('truth'),df1['total.query'].rename('query')]
                ,axis=1).plot.bar(title ="varscan vcf records", figsize=(15,10),legend=True,
                                  color=['g', 'yellow'], rot=60, fontsize=8)
ax1.set_ylabel("Records")
plt.show()

#df1 = dataframe.loc[dataframe['name'#plt.show()].str.endswith('0.1')]
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