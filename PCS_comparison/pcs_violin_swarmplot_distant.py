

import pandas as pd
import numpy as np

df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-calJac4/PCS_hg38calJac4.csv") 
df2 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-danRer10/PCS_hg38danRer10.csv") 
df3 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-equCab3/PCS_hg38equCab3.csv") 
df4 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-galGal6/PCS_hg38galGal6.csv") 
df5 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-ornAna2/PCS_hg38ornAna2.csv") 
df6 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-oryCun2/PCS_hg38oryCun2.csv") 
df7 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-petMar3/PCS_hg38petMar3.csv") 
df8 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-tarSyr2/PCS_hg38tarSyr2.csv") 
df9 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-thaSir1/PCS_hg38thaSir1.csv") 
df10 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-triMan1/PCS_hg38triMan1.csv") 
df11 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-xenTro10/PCS_hg38xenTro10.csv") 

df1['Query Species'] = 'Marmoset'
df2['Query Species'] = 'Zebrafish'
df3['Query Species'] = 'Horse'
df4['Query Species'] = 'Chicken'
df5['Query Species'] = 'Platypus'
df6['Query Species'] = 'Rabbit'
df7['Query Species'] = 'Lamprey'
df8['Query Species'] = 'Tarsier'
df9['Query Species'] = 'Garter Snake'
df10['Query Species'] = 'Manatee'
df11['Query Species'] = 'X. tropicalis'

df = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11])

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

ax[0].set_ylabel(ax[0].get_ylabel(), fontsize=14)
ax[1].set_ylabel(ax[1].get_ylabel(), fontsize=14)
ax[0].set_xlabel(ax[0].get_xlabel(), fontsize=14)
ax[1].set_xlabel(ax[1].get_xlabel(), fontsize=14)

ax[0].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
ymin, ymax = ax[0].get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
ax[0].yaxis.set_ticks(tick_range)
ax[0].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)

plt.tight_layout()
plt.savefig('/bucket/.mabuya/MillerU/Abrar/PCS_comparison/All_PCS_violin_swarmplot_distant.png')


