#!/usr/bin/env python
import numpy as np
import os
import pandas as pd
os.chdir("/home/kevin/PycharmProjects/begins/continues/control2/AF_DP1")

highdf = pd.read_csv('high_barplot.csv', sep=',')
moddf = pd.read_csv('mod_barplot.csv', sep=',')
lowdf = pd.read_csv('low_barplot.csv', sep=',')


data=pd.merge(highdf, moddf, on="Gene",how='left')
datadf=pd.merge(data, lowdf, on="Gene",how='left')
datadf=datadf.fillna(value=0)
datadf = datadf.set_index('Gene')
#different sorts###########################################
totaldf=datadf.copy()
totaldf['total'] = totaldf.sum(axis=1)
totaldf=totaldf.sort_values(by=['total'],ascending=False)

highmoddf=datadf.copy()
highmoddf['total'] = highmoddf['donors_x']+highmoddf['donors_y']
highmoddf=highmoddf.sort_values(by=['total'],ascending=False)
edatadf=datadf[:25]
fdatadf=highmoddf[:25]
gdatadf=totaldf[:25]
ax1=pd.concat([edatadf['donors_x'].rename('high'),edatadf['donors_y'].rename('moderate'),edatadf['donors'].rename('low')],axis=1).plot.bar(title ="Top 25 Genes by donors with a least one high functional impact mutation", figsize=(15,10),legend=True, color=['r', 'yellow', 'g'], rot=60, fontsize=8)
ax1.set_ylabel("Donors with at least 1 mutation")
ax2=pd.concat([fdatadf['donors_x'].rename('high'),fdatadf['donors_y'].rename('moderate'),fdatadf['donors'].rename('low')],axis=1).plot.bar(title ="Top 25 Genes by sum of High and Moderate functional impact mutations occuring at least once in donors", figsize=(15,10),legend=True, color=['r', 'yellow', 'g'], rot=60, fontsize=8)
ax2.set_ylabel("Donors with at least 1 mutation")
ax3=pd.concat([gdatadf['donors_x'].rename('high'),gdatadf['donors_y'].rename('moderate'),gdatadf['donors'].rename('low')],axis=1).plot.bar(title ="Top 25 Genes by donors with a least one mutation", figsize=(15,10),legend=True, color=['r', 'yellow', 'g'], rot=60, fontsize=8)
ax3.set_ylabel("Donors with at least 1 mutation")

