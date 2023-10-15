

import pandas as pd
import numpy as np

#df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38gorGor6/All_bin_15kb/quantile_kurtosis_100pc.csv") 
#df2 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38mm39/All_bin_15kb/quantile_kurtosis_100pc.csv") 
#df3 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38galVar1/All_bin_15kb/quantile_kurtosis_100pc.csv") 
#df4 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38micMur3/All_bin_15kb/quantile_kurtosis_100pc.csv") 
#df5 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38ponAbe3/All_bin_15kb/quantile_kurtosis_100pc.csv") 
#df6 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38rheMac10/All_bin_15kb/quantile_kurtosis_100pc.csv") 
#df7 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38nomLeu3/All_bin_15kb/quantile_kurtosis_100pc.csv") 
#df8 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38saiBol1/All_bin_15kb/quantile_kurtosis_100pc.csv") 

#df1['Query Species'] = 'Gorilla'
#df2['Query Species'] = 'Mouse'
#df3['Query Species'] = 'Malayan Lemur'
#df4['Query Species'] = 'Mouse Lemur'
#df5['Query Species'] = 'Orangutan'
#df6['Query Species'] = 'Rhesus'
#df7['Query Species'] = 'Gibbon'
#df8['Query Species'] = 'Squirrel Monkey'

df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38calJac4/All_bin_15kb/quantile_kurtosis_100pc.csv") 
df2 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38danRer10/All_bin_15kb/quantile_kurtosis_100pc.csv") 
df3 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38equCab3/All_bin_15kb/quantile_kurtosis_100pc.csv") 
df4 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38galGal6/All_bin_15kb/quantile_kurtosis_100pc.csv") 
df5 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38ornAna2/All_bin_15kb/quantile_kurtosis_100pc.csv") 
df6 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38oryCun2/All_bin_15kb/quantile_kurtosis_100pc.csv") 
df7 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38petMar3/All_bin_15kb/quantile_kurtosis_100pc.csv") 
df8 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38tarSyr2/All_bin_15kb/quantile_kurtosis_100pc.csv") 
df9 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38thaSir1/All_bin_15kb/quantile_kurtosis_100pc.csv") 
df10 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38triMan1/All_bin_15kb/quantile_kurtosis_100pc.csv") 
df11 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/hg38xenTro10/All_bin_15kb/quantile_kurtosis_100pc.csv") 

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

df = df.rename(columns={'0': 'Kurtosis'})

import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
import seaborn as sns

fig, ax = plt.subplots(ncols=2, figsize=(10, 5), sharey=True)

sns.violinplot(x= 'Query Species', y='Kurtosis', data=df, ax=ax[0])
sns.stripplot(x= 'Query Species', y='Kurtosis', data=df, ax=ax[1])

ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=90, fontsize=11)
ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=90, fontsize=11)

ax[0].set_yticklabels(ax[0].get_yticklabels(), fontsize=11)
ax[1].set_yticklabels(ax[1].get_yticklabels(), fontsize=11)

ax[0].set_ylabel(ax[0].get_ylabel(), fontsize=13)
ax[1].set_ylabel(ax[1].get_ylabel(), fontsize=13)
ax[0].set_xlabel(ax[0].get_xlabel(), fontsize=13)
ax[1].set_xlabel(ax[1].get_xlabel(), fontsize=13)

ax[0].yaxis.set_major_formatter(mticker.StrMethodFormatter("${{{x:.0f}}}$"))
ymin, ymax = ax[0].get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
ax[0].yaxis.set_ticks(tick_range, minor=True)

plt.tight_layout()
plt.savefig('/bucket/.mabuya/MillerU/Abrar/PCS_comparison/All_kurtosis_distant_violin_strip_plot_15kb.png')


