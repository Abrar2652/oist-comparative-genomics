#!/usr/bin/env python
# coding: utf-8

## Red in case of Mouse, Blue in case of Gorilla
import os
import glob
import pandas as pd

lines=[]
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
    with open("/bucket/MillerU/Abrar/hg38-saiBol1_outputs/hg38-saiBol1_tiled_pcs_15kb/chr{0}_15kb_min_window_tiled_pcs.txt".format(i)) as file:
        while True:
            # Get next line from file
            line = file.readline().rstrip('\n')
            # if line is empty
            # end of file is reached
            if not line:
                break
            elif line==("\n"):
                continue
            #print(line)
            lines.append([line])
lines_frame = pd.DataFrame(lines, columns=['PCS'])

import pandas as pd
# kurtosis_100pc_500col.csv" in case of gorGor3
#nash_kurtosis = pd.read_csv("/bucket/MillerU/Abrar/Moment_val_bins/hg38saiBol1/All_bin_15kb/kurtosis_100pc_500col.csv") #in the Moment_val_bins folder, the kurtosis files are located
nash_kurtosis = pd.read_csv("/bucket/MillerU/Abrar/Moment_val_bins/hg38saiBol1/All_bin_15kb/quantile_kurtosis_100pc.csv") 
nash_kurtosis = nash_kurtosis[nash_kurtosis['0'] > 0] #nash_kurtosis['0'] != 0  #Manually convert excel to numeric (upto 9 decimal places) if any error is encountered

# In[92]:

st2 = 'PCS distribution (with exons) having Heavy Tail (right hand side cluster) in Quantile Kurtosis'
species1 = 'hg38'
species2 = 'saiBol1'
percent = '100' #change accordingly

# In[40]:

cutoff = 6.25
backtrack_index = nash_kurtosis[nash_kurtosis['0'] > cutoff].index.tolist() #take decision while selecting this 'cutoff' #high kurtosis
#backtrack_index = nash_kurtosis[nash_kurtosis['0'] < cutoff].index.tolist() #take decision while selecting this 'cutoff' #low kurtosis

# In[78]:

tail_pcs=[]
for i in backtrack_index:
    tail_pcs.append(lines_frame.iloc[i].values.tolist())
flat_tail_pcs = [item for sublist in tail_pcs for item in sublist]
final_tail_pcs = [int(i) for i in ','.join(flat_tail_pcs).split(',')]


# In[87]:
serial_count=[]
for i in range(0, int(max(final_tail_pcs)+5)):
    serial_count.append(final_tail_pcs.count(i))


lmer = [x for x in range(int(max(final_tail_pcs)+5))]
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}        
df = pd.DataFrame(dict) 
df.to_csv(r'/bucket/MillerU/Abrar/PCS_heavy_tailed/high_kurtosis_PCS_counts_hg38saiBol1.csv', index = False) 


# Kurtosis verification
import numpy as np
verify_new = [[int(x) for x in sublst[0].split(',')] for sublst in tail_pcs]
# for filtered and vanilla hggg, 60 (L>=20), otherwise 56 (full); for filtered hgmm, 6 (filtered L>=10), 11 (vanilla L>=20), otherwise 3 (full)
denominator = float(11) # IQR
temp = []
for sublist in verify_new:
    if (denominator!=0) and i:
        temp.append(round(float((np.quantile(sublist, 0.99)-np.quantile(sublist, 0.01))/denominator), 1)) # Nash's Kurtosis Definition
    elif not i:
        temp.append(round(float(0), 2))
    else:
        temp.append(round(float(0), 2))


# This loop will print False if there is a single mismatch in the Kurtosis values of each window

# high kurtosis
for ind, val in enumerate(temp):
    if (val != round(nash_kurtosis[nash_kurtosis['0'] > cutoff]['0'], 1).values[ind]):
        print(val,'\t',round(nash_kurtosis[nash_kurtosis['0'] > cutoff]['0'], 1).values[ind])

'''
#low kurtosis
for ind, val in enumerate(temp):
    if (val != round(nash_kurtosis[nash_kurtosis['0'] < cutoff]['0'], 1).values[ind]):
        print(val,'\t',round(nash_kurtosis[nash_kurtosis['0'] < cutoff]['0'], 1).values[ind])       
      
'''

'''
# In[93]:

import numpy as np
import matplotlib.pyplot as plt

x = [x for x in range(int(max(final_tail_pcs)+5))][10:]
y = serial_count[10:]

plt.figure(figsize = [15, 10])

plt.plot(x ,y, '.g',linewidth=2, label='Frequency of PCS')
plt.title('{} Distribution with Cutoff = {} (linear plot)\n{}/{}'.format(st2,cutoff,species1,species2), fontsize=15)
plt.xlabel('PCS', fontsize=13)
plt.ylabel('Counts', fontsize=13)

plt.legend()
plt.savefig('/bucket/MillerU/Abrar/PCS_heavy_tailed/' + species1 + species2 + '_All_' + percent + 'pc_cutoff' + str(cutoff) + '_tiled_pcs_linear.png')

# In[94]:


import numpy as np
import matplotlib.pyplot as plt

x = [x for x in range(int(max(final_tail_pcs)+5))][10:]
y = serial_count[10:]

plt.figure(figsize = [15, 10])

plt.loglog(x ,y, '.g',linewidth=2, label='Frequency (log10) of PCS (log10)')
plt.title('{} Distribution with Cutoff = {} (log-log plot)\n{}/{}'.format(st2,cutoff,species1,species2), fontsize=15)
plt.xlabel('PCS (log10)', fontsize=13)
plt.ylabel('Counts (log10)', fontsize=13)

plt.legend()

plt.savefig('/bucket/MillerU/Abrar/PCS_heavy_tailed/' + species1 + species2 + '_All_' + percent + 'pc_cutoff' + str(cutoff) + '_tiled_pcs_loglog.png')


# In[96]:


import numpy as np
import matplotlib.pyplot as plt

x = [x for x in range(int(max(final_tail_pcs)+5))][10:]
y = serial_count[10:]

plt.figure(figsize = [15, 10])

plt.semilogy(x , y, '.g', linewidth=2, base=10, label='Frequency (log10) of PCS')
plt.title('{} Distribution with Cutoff = {} (semilog(y) plot)\n{}/{}'.format(st2,cutoff,species1,species2), fontsize=15)
plt.xlabel('PCS', fontsize=13)
plt.ylabel('Counts (log10)', fontsize=13)

plt.legend()

plt.savefig('/bucket/MillerU/Abrar/PCS_heavy_tailed/' + species1 + species2 + '_All_' + percent + 'pc_cutoff' + str(cutoff) + '_tiled_pcs_semilogy.png')

'''
