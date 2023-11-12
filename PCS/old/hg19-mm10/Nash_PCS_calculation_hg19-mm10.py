#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import glob
import pandas as pd
os.chdir("/bucket/MillerU/Abrar/Nash_wrong_PCS/hg19-mm10/all/")

# In[2]:

extension = 'tsv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]

# In[4]:

import natsort
all_filenames = natsort.natsorted(all_filenames, reverse=False)

# In[5]:

df1 = pd.concat( [pd.read_csv(f, sep='\t') for f in all_filenames[0:50]], sort=False)#1->50
df2 = pd.concat( [pd.read_csv(f, sep='\t') for f in all_filenames[50:100]], sort=False)#51->100
df3 = pd.concat( [pd.read_csv(f, sep='\t') for f in all_filenames[100:150]], sort=False)#101->150

# In[14]:

df4 = pd.concat( [pd.read_csv(f, sep='\t') for f in all_filenames[150:200]], sort=False)#151->200
df5 = pd.concat( [pd.read_csv(f, sep='\t') for f in all_filenames[200:250]], sort=False)#201->250
df6 = pd.concat( [pd.read_csv(f, sep='\t') for f in all_filenames[250:300]], sort=False)#251->300
df7 = pd.concat( [pd.read_csv(f, sep='\t') for f in all_filenames[300:364]], sort=False)#301->364


# In[16]:

total_rows = len(df1) + len(df2) + len(df3) + len(df4) + len(df5) + len(df6) + len(df7)
print(total_rows)

# In[20]:

seq_count = pd.concat([df1['first.width'], df2['first.width'], df3['first.width'], df4['first.width'], df5['first.width'], df6['first.width'], df7['first.width']])

# In[21]:

print(seq_count)


# In[23]:

del(df1)
del(df2)
del(df3)
del(df4)
del(df5)
del(df6)
del(df7)

# In[24]:
serial_count = seq_count.value_counts()

'''
import numpy as np
seq_count = seq_count['first.width'].values.tolist()
serial_count = []
for i in range(int(max(seq_count))+5):
    serial_count.append(seq_count.count(i))
'''

print(serial_count)


# In[30]:

'''
import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0, int(max(seq_count)+5), 1)
y = serial_count

plt.figure(figsize = [8, 6])

plt.loglog(x ,y, '.r',linewidth=2, label='hg19/mm10 L-mers')
plt.title('Length distributions of perfectly conserved L-mers from Human/Mouse genome alignment for L>=1 Nash\'s log-log plot', fontsize=15)
plt.xlabel('Length (L)', fontsize=13)
plt.ylabel('Number of L-mers', fontsize=13)

plt.legend()

plt.savefig('/bucket/MillerU/Abrar/PCS/hg19-mm10/Nash_PCS_hg19mm10.png')
'''


# In[27]:


#lmer = [x for x in range(int(max(seq_count)+5))]
# dictionary of lists  
lmer = serial_count.index
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
df = pd.DataFrame(dict) 
df.to_csv(r'/bucket/MillerU/Abrar/PCS/Nash_PCS_calculation_hg19mm10.csv', index=False)







