

import pandas as pd
import numpy as np

df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-gorGor6/PCS_hg38gorGor6.csv") 
df2 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-mm39/PCS_hg38mm39.csv") 
df3 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-galVar1/PCS_hg38galVar1.csv") 
df4 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-micMur3/PCS_hg38micMur3.csv") 
df5 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-ponAbe3/PCS_hg38ponAbe3.csv") 
df6 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-rheMac10/PCS_hg38rheMac10.csv") 
df7 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-nomLeu3/PCS_hg38nomLeu3.csv") 
df8 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-saiBol1/PCS_hg38saiBol1.csv") 

df1['Query Species'] = 'Gorilla'
df2['Query Species'] = 'Mouse'
df3['Query Species'] = 'Malayan Lemur'
df4['Query Species'] = 'Mouse Lemur'
df5['Query Species'] = 'Orangutan'
df6['Query Species'] = 'Rhesus'
df7['Query Species'] = 'Gibbon'
df8['Query Species'] = 'Squirrel Monkey'

df = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8])

df = df[df['Number of PCS of Length L'] != 0]

df['Number of PCS of Length L (log10)'] = np.log10(df['Number of PCS of Length L'])

df = df[df['Number of PCS of Length L (log10)'] != 0]

import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
import seaborn as sns

#sns.catplot(x='species', y='Number of PCS of Length L', data=df, kind = 'violin')
log_df = df

fig, ax = plt.subplots(ncols=2, figsize=(10, 5), sharey=True)
#plt.xticks(rotation=45)

sns.violinplot(x= 'Query Species', y='Number of PCS of Length L (log10)', data=log_df, ax=ax[0])
sns.swarmplot(x= 'Query Species', y='Number of PCS of Length L (log10)', data=log_df, s=3, ax=ax[1])

ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=90, fontsize=12)
ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=90, fontsize=12)

ax[0].set_yticklabels(ax[0].get_yticklabels(), fontsize=12)
ax[1].set_yticklabels(ax[1].get_yticklabels(), fontsize=12)

ax[0].set_ylabel(ax[0].get_ylabel(), fontsize=15)
ax[1].set_ylabel(ax[1].get_ylabel(), fontsize=15)
ax[0].set_xlabel(ax[0].get_xlabel(), fontsize=15)
ax[1].set_xlabel(ax[1].get_xlabel(), fontsize=15)

ax[0].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
ymin, ymax = ax[0].get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
ax[0].yaxis.set_ticks(tick_range)
ax[0].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)

plt.tight_layout()
plt.savefig('/bucket/.mabuya/MillerU/Abrar/PCS_comparison/All_PCS_violin_swarm_splitplot.png')


