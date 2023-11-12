

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

############ PLOTS ##################
cutoff = 100
name = 'Cat'
df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/idea3/PCS_hg38-mm39-felCat9_merged_extended_L{}.csv".format(cutoff))

# Importing necessary libraries


# Preparing the data for the plot
x = [i for i in df1.index][:]
y = df1['Number of PCS of Length L'][:]


# Resizing the figure
plt.figure(figsize=[10, 7])

# Plotting the graph with Log ticks at x and y axis using loglog
plt.loglog(x, y, '.r', linewidth=2, label='Merged PCS of Human/Mouse and Human/{0}'.format(name))
plt.title('PCS Length Distribution of Merged PCS of hg38/mm39_felCat9 (log-log) for L>={0}'.format(cutoff), fontsize=14)

plt.xlabel('Length (log10)', fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plt.legend(title="Pairwise Distants Used:\n Human/Mouse, Human/Cat", title_fontsize='large', fontsize=13, markerscale=3, loc='upper left')
plt.savefig('/bucket/MillerU/Abrar/PCS/idea3/PCS_hg38mm39felCat9_loglog_L{}.png'.format(cutoff))





