#!/usr/bin/env python
# coding: utf-8

# In[206]:


import os
import glob
import pandas as pd
#os.chdir("D:/Research/OIST/Project/Human-Gorilla/All_bin_30kb")
#os.chdir("D:/Research/OIST/Project/Human-Mouse/Human Mouse 2009/All_bin_30kb")

# In[222]:

import pandas as pd
nash_kurtosis = pd.read_csv("/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_30kb/kurtosis_100pc_500col.csv") #kurtosis_100pc_500col #Only working with 30kb
nash_kurtosis = nash_kurtosis[nash_kurtosis['0'] > 0] #nash_kurtosis['0'] != 0 
stat = pd.DataFrame({'counts': nash_kurtosis['0'].value_counts().sort_index()})

# In[224]:

st = 'Quantile Kurtosis' #stat to calculate
st2 = 'PCS distribution having Heavy Tail in Quantile Kurtosis'
species1 = 'hg19'
species2 = 'gorGor3'
window_size = 30 #in kb

############################ SCATTER PLOTS ##########################################

Outdir = ('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_{}kb/').format(window_size)

# In[8]:

import numpy as np
import matplotlib.pyplot as plt

x = stat.index
y = stat.values

plt.figure(figsize = [15, 10])

plt.plot(x, y, '.b', linewidth=2, label='Frequency of {}'.format(st))
plt.title('{} Distribution\n{}/{}'.format(st, species1, species2), fontsize=15)
plt.xlabel('{}'.format(st), fontsize=13)
plt.ylabel('Counts', fontsize=13)

plt.legend()

plt.savefig((Outdir + 'Scatterplots/{}_100pc_500col_linear.png').format(st))
#plt.savefig((Outdir + 'Scatterplots/{}_100pc_linear.png').format(st))


# In[9]:

import numpy as np
import matplotlib.pyplot as plt

x = stat.index
y = stat.values

plt.figure(figsize = [15, 10])

plt.loglog(x ,y, '.b',linewidth=2, label='Frequency (log10) of {} (log10)'.format(st))
plt.title('{} Distribution (log-log)\n{}/{}'.format(st, species1, species2), fontsize=15)
plt.xlabel('{} (log10)'.format(st), fontsize=13)
plt.ylabel('Counts (log10)', fontsize=13)

plt.legend()

plt.savefig((Outdir + 'Scatterplots/{}_100pc_500col_loglog.png').format(st))
#plt.savefig((Outdir + 'Scatterplots/{}_100pc_loglog.png').format(st))


# In[10]:

import numpy as np
import matplotlib.pyplot as plt

x = stat.index
y = stat.values

plt.figure(figsize = [15, 10])

plt.semilogy(x , y, '.b', linewidth=2, base=10, label='Frequency (log10) of {}'.format(st))
plt.title('{} Distribution semilog(y)\n{}/{}'.format(st, species1, species2), fontsize=15)
plt.xlabel('{}'.format(st), fontsize=13)
plt.ylabel('Counts (log10)', fontsize=13)

plt.legend()

plt.savefig((Outdir + 'Scatterplots/{}_100pc_500col_semilogy.png').format(st))
#plt.savefig((Outdir + 'Scatterplots/{}_100pc_semilogy.png').format(st))



################################ HISTOGRAMS ##########################################


# In[246]:


# Import the libraries
import matplotlib.pyplot as plt

plt.figure(figsize=[15,10])
# matplotlib histogram
bin = 300
plt.hist(nash_kurtosis['0'], color = 'red', edgecolor='black',
         bins = bin)

# Add labels
plt.title("{} Distribution \n{}/{} \nNumber of bins = {}".format(st, species1, species2, bin), fontsize=15)
plt.xlabel('{}'.format(st), fontsize=13)
plt.ylabel('Counts', fontsize=13)

plt.savefig((Outdir + 'Histograms/{}_100pc_500col_{}.png').format(st,bin)) #gorilla
#plt.savefig((Outdir + 'Histograms/{}_100pc_{}.png').format(st,bin)) #others


# In[277]:


# Import the libraries
import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=[10,8])

# histogram on log scale. 
# Use non-equal bin sizes, such that they look equal on log scale.
def plot_loghist(x, bin):
    hist, bins = np.histogram(x, bins=bin)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    plt.hist(x, color = 'blue', edgecolor='black', bins=logbins)
    plt.xscale('log')
    plt.yscale('log')
    # Add labels
    plt.title("{} Distribution (log-log plot) \n{}/{} \nNumber of bins = {}".format(st, species1, species2, bin), fontsize = 15)
    plt.xlabel('{} (log10)'.format(st), fontsize = 13)
    plt.ylabel('Counts (log10)', fontsize = 13)
    plt.savefig((Outdir + 'Histograms/{}_100pc_500col_loglog_{}.png').format(st,bin)) #gorilla
    #plt.savefig((Outdir + 'Histograms/{}_100pc_loglog_{}.png').format(st,bin)) #others
    
plot_loghist(nash_kurtosis['0'], 200)


# In[263]:


# Import the libraries
import matplotlib.pyplot as plt

plt.figure(figsize=[15,10])
# matplotlib histogram
bins = 200
plt.hist(nash_kurtosis['0'], color = 'red', edgecolor='black',
         bins = bins)
plt.yscale('log')
# Add labels
plt.title("{} Distribution (semilog) \n{}/{} \nNumber of bins = {}".format(st, species1, species2, bins), fontsize=15)
plt.xlabel('{}'.format(st), fontsize=13)
plt.ylabel('Counts (log10)', fontsize=13)

plt.savefig((Outdir + 'Histograms/{}_100pc_500col_semilog_{}.png').format(st,bins)) #gorilla
#plt.savefig((Outdir + 'Histograms/{}_100pc_semilog_{}.png').format(st,bins)) #others








