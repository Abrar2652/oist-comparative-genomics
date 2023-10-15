

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

############ PLOTS ##################
cutoff = 235
name = 'Orangutan'
df4 = pd.read_csv("/bucket/MillerU/Abrar/PCS/PCS_hg38-gg6-ponAbe3_{}.csv".format(cutoff))

# Importing necessary libraries


# Preparing the data for the plot

x = [i for i in df4.index][1:]
y = df4['Number of PCS of Length L'][1:]

# Resizing the figure
plt.figure(figsize=[10, 7])
plt.ylim(0.5, 10**8)
# Plotting the graph with Log ticks at x and y axis using loglog
plt.loglog(x, y, '.r', linewidth=2, 
           label='Human/Gorilla_{} whole genome L-mers'.format(name))
plt.title('Length distributions of Perfectly Conserved\n L-mers from hg38/gorGor6 (L>=1)_ponAbe3 (L>={}) genome alignment (log-log)'.format(cutoff), fontsize=14)
plt.xlabel('Length (log10)', fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plt.legend(fontsize=15, markerscale=3, loc='upper right')
plt.savefig('/bucket/MillerU/Abrar/PCS/PCS_hg38gg6ponAbe3_loglog_L{}.png'.format(cutoff))



######## semilogy ######

plt.figure(figsize=[10, 7])
plt.ylim(0.5, 10**8)
plt.semilogy(x , y, '*m', linewidth=2, base=10, label='Human/Gorilla_{} whole genome L-mers'.format(name))
plt.xlabel("Length (L)", fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)


plt.title("Length distributions of Perfectly Conserved\n L-mers from hg38/gorGor6 (L>=1)_ponAbe3 (L>={}) genome alignment (semi-log)".format(cutoff), fontsize=14)

plt.legend(fontsize=15, markerscale=3, loc='upper right')

plt.savefig('/bucket/MillerU/Abrar/PCS/PCS_hg38gg6ponAbe3_semilog_L{}.png'.format(cutoff))


