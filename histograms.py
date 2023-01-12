#!/usr/bin/env python
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
os.chdir("/home/kevin/mrd_wes/results/2018-05-09")
#df = pd.read_csv('23032018_summary_mutations.csv', sep=',')
#df = pd.read_csv('newcol_test.csv', sep=',')
df = pd.read_csv('nodup_mutation_impact_effects_sort.csv', sep=',')
#df4=df[['dp','mq','mmq']]
df1=df['dp']
df2=df['mq']
df3=df['mmq']
df4=df['af']
df5=df['vd']
#df2=df2.fillna(value=0)
#df3=df3.fillna(value=0)
plt.figure();

#ax1=df1.hist(stacked=True, bins=[0,20,40,60,80,100,120,140,160,180,200])

ax1=df1.hist(stacked=True, bins=10, range=(0,400),grid=False)
ax1.set_xlabel("DP in width 40 bins")
ax1.set_ylabel("Frequency")
plt.title('Total Depth (DP) of Reads')


#marker='.',markersize=10

#plt.figure();
#ax1=df1.hist(stacked=True, bins=[0,30,800])
#ax1.set_xlim((0, 60))

plt.figure();
ax1=df1.hist(stacked=True, bins=[0,40,800],grid=False)
ax1.set_xlim((0, 80))
plt.xticks([], [])
ax1.set_xlabel("DP below 40                                DP above 40")
ax1.set_ylabel("Frequency")
plt.title('Total Depth (DP) of Reads')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

#self.axes.set_xticks(np.arange(0,6,1))
plt.figure();
ax1=df1.hist(stacked=True, bins=[0,40,800])
ax1.set_xlim((0, 80))

plt.figure();
ax2=df2.plot.hist(stacked=True, bins=[0,10,20,30,40,50,60])
ax2.set_xlabel("MQ")
ax2.set_ylabel("Frequency")
plt.title('Mean Map Quality (MQ)')

plt.figure();
ax3=df3.plot.hist(stacked=True, bins=[0,10,20,30,40,50,60])
ax3.set_xlabel("MMQ")
ax3.set_ylabel("Frequency")
plt.title('Median Map Quality (MMQ)')

plt.figure();
ax4=df4.plot.hist(stacked=True, bins=10)#, range=(0,1))
ax4.set_xlabel("AF")
ax4.set_ylabel("Frequency")
plt.title('Allelle Frequency (AF)')

plt.figure();
ax5=df5.plot.hist(stacked=True, bins=5, range=(0,50))
ax5.set_xlabel("VD")
ax5.set_ylabel("Frequency")
plt.title('Variant Depth (VD)')

