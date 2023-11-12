#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import glob
import pandas as pd
os.chdir("/bucket/MillerU/Abrar/Nash_wrong_PCS/hg19-gorGor3/all/")

# In[2]:

extension = 'tsv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]

# In[4]:

import natsort
all_filenames = natsort.natsorted(all_filenames, reverse=False)

# In[5]:

df1 = pd.concat( [pd.read_csv(f, sep='\t') for f in all_filenames[0:50]], sort=False)#1->50
df2 = pd.concat( [pd.read_csv(f, sep='\t') for f in all_filenames[50:83]], sort=False)#51->83

# In[16]:

total_rows = len(df1) + len(df2)
print(total_rows)

# In[20]:

seq_count = pd.concat([df1['first.width'], df2['first.width']])


# In[21]:

#print(seq_count)

# In[23]:

del(df1)
del(df2)

# In[24]:

serial_count = seq_count.value_counts()

'''
import numpy as np
seq_count = seq_count['first.width'].values.tolist()
serial_count = []
for i in range(int(max(seq_count))+5):
    serial_count.append(seq_count.count(i))
'''

#print(serial_count)

# In[30]:

'''
import numpy as np
import matplotlib.pyplot as plt

x = serial_count.index
y = serial_count

plt.figure(figsize = [8, 6])

plt.loglog(x ,y, '.r',linewidth=2, label='hg19/gorGor3 L-mers')
plt.title('Length distributions of perfectly conserved L-mers from Human/Gorilla genome alignment for L>=1 Nash\'s log-log plot', fontsize=15)
plt.xlabel('Length (L)', fontsize=13)
plt.ylabel('Number of L-mers', fontsize=13)

plt.legend()

plt.savefig('/bucket/MillerU/Abrar/PCS/hg19-gorGor3/Nash_PCS_hg19gorGor3.png')
'''

# In[27]:

#lmer = [x for x in range(int(max(seq_count)+5))]
lmer = serial_count.index
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
df = pd.DataFrame(dict) 
df.to_csv(r'/bucket/MillerU/Abrar/PCS/hg19-gorGor3/Nash_PCS_calculation_hg19gorGor3.csv', index=False)






