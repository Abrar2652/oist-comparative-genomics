#!/usr/bin/env python
# coding: utf-8


# In[45]:

#pip3 install bx-python
import bx.align.axt
import os
count = 0
#file = open("chr13.hg19.mm10.net.axt")
lines=[]
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']: # reading all genome AXT files
    with open('/bucket/MillerU/Abrar/axtNet/hg19gorGor3/chr{0}.hg19.gorGor3.net.axt'.format(i)) as file:
        while True: 
            # Get next line from file
            line = bx.align.axt.readline(file, skip_blank=False) # function for read the necessary AXT lines

            # if line is empty
            # end of file is reached
            if not line:
                break
            elif line==("\n"):
                continue
            count += 1
            # Keep only those lines containing the sequences and don't include unnecessary lines
            # (count%3-1) generates a sequence of [0,1,-1,0,1,-1,0,1,-1,...] for [1,2,3,4,5,...]
            # The first line of AXT file is a header info, 2nd and 3rd line containes pairwise sequences. 
            # If we need every 2nd and 3rd line then (count%3-1)!=0 will give us the desired sequences
            if (count%3-1)!=0:
                #print("Line{}: {}".format(count, line.strip()))
                lines.append(line)

#while True: 
    # Get next line from file
#    line = bx.align.axt.readline(file, skip_blank=False)
 
    # if line is empty
    # end of file is reached
#    if not line:
#        break
#    elif line==("\n"):
#        continue
#    count += 1
#    if (count%3-1)!=0:
       # print("Line{}: {}".format(count, line.strip()))
#        lines.append(line)
        


# In[46]:


import pandas as pd
lines_frame = pd.DataFrame(lines, columns=['Genome Sequence'])
lines_frame


# In[47]:


# concatenates all the human genomes obtained from the even rows of the AXT file's dataframe and removes newlines
#lines_frame.iloc[::2,:]['Genome Sequence'].str.cat(sep='')
#len(lines_frame.iloc[::2,:]['Genome Sequence'].str.cat(sep=''))
hum = ''.join(lines_frame.iloc[::2,:]['Genome Sequence']).replace("\n", "")
mouse = ''.join(lines_frame.iloc[1::2,:]['Genome Sequence']).replace("\n", "")


# In[48]:
# Change the bin_size as you wish
bin_size=30000/4

import math
# If we divide the total number of bases by the bin_size, we'll get the number of bins and ceiling function will maximize it
num_bins = math.ceil(len(hum)/bin_size)
print(num_bins)


# In[49]:


# Dividing the bases into 30kb equal bins/windows
# chunk function divides all the first bases into the equal number of windows except the last one 
# The last one contains the remaining bases, so no bases left out for binning
# In this way both sequences get comparibility because of having equal sizes
def chunk(in_string,num_chunks):
    chunk_size = len(in_string)//num_chunks
    if len(in_string) % num_chunks: chunk_size += 1
    iterator = iter(in_string)
    for _ in range(num_chunks):
        accumulator = list()
        for _ in range(chunk_size):
            try: accumulator.append(next(iterator))
            except StopIteration: break
        yield ''.join(accumulator)
hum_bin = list(chunk(hum, num_bins))
mouse_bin = list(chunk(mouse, num_bins))
total_bin = len(hum_bin) + len(mouse_bin)
print(total_bin)



# In[52]:


import numpy as np
res = np.arange(0, total_bin)
df = pd.DataFrame(res, columns=['Genome Sequence'])
# Each binned hum  and mouse sequences are stored in the even and odd rows respectively like earlier format
df.loc[0::2,"Genome Sequence"] = hum_bin #even rows
df.loc[1::2,"Genome Sequence"] = mouse_bin #odd rows



# In[54]:

import re
def compare_strings(df, column, j):
    # First, I removed the split... it is already an array
    str1 = df["Genome Sequence"][j] #even rows
    str2 = df["Genome Sequence"][j+1] #odd rows

    #then creating a new variable to store the result after  
    #comparing the strings. You note that I added result2 because 
    #if string 2 is longer than string 1 then you have extra characters 
    #in result 2, if string 1 is  longer then the result you want to take 
    #a look at is result 2

    result1 = ''
    result2 = ''

    #handle the case where one string is longer than the other
    maxlen=len(str2) if len(str1)<len(str2) else len(str1)

    #loop through the characters
    for i in range(maxlen):
      #use a slice rather than index in case one string longer than other
        letter1=str1[i:i+1]
        letter2=str2[i:i+1]
        #create string with differences
        # Compares and keeps only uppercase A, T, C, G. Otherwise, discards other letters and lowercase letters as well.
        if (letter1.islower() and letter2.islower() and letter1==letter2) or (((letter1 not in ['A','T','C','G']) or (letter2 not in ['A','T','C','G'])) and letter1==letter2)  or ((letter1 != letter2) or (letter1=='-') or (letter2=='-')):
            result1+='\n'
            result2+='\n'
        if ((letter1 == letter2) and letter1.isupper() and letter2.isupper() and (letter1!='-') and (letter2!='-') and 
            letter1 in ['A','T','C','G'] and letter2 in ['A','T','C','G']):
            result1+=letter1
            result2+=letter2

    word=re.sub(r'\n+', '\n', result1).strip()
    return word


ultimate = []
total_count = []
for i in range(int(len(df))-1):
    seq_count = []
    if i%2==0:
        for x in (compare_strings(df, "Genome Sequence", i).split()):
            #print(x)
            seq_count.append(len(x))
            total_count.append(len(x))
        ultimate.append(seq_count)


# In[56]:

ultimate2 = []
count = 1
for i in range(int(len(df))-1): # iterates through each sequences and compares each pair
    seq_count2 = []
    if i%2==0:
        seq_count2.append("bin_num = {}, window_size = {}".format(count, bin_size)) # Stores the bin number and window size information
        for x in (compare_strings(df, "Genome Sequence", i).split()):
            #print(x)
            seq_count2.append(len(x))
        count+=1
        ultimate2.append(seq_count2)

# In[57]:


#setting fisher=False in the above code
#does the calculation of the Pearson’s definition of 
#kurtosis where the kurtosis value for normal distribution = 3.
from scipy.stats import kurtosis
from scipy.stats import skew
from scipy.stats import moment
kurt_val = []
mean_val = []
variance_val = []
skew_val = []
sixth_mom = []
denominator = float(np.quantile(total_count, 0.75)-np.quantile(total_count, 0.25))
print("Denominator: ",denominator)
for i in [[0] if not arr else arr for arr in ultimate]:
    if (denominator!=0) and i:
        x=float((np.quantile(i, 0.99)-np.quantile(i, 0.01))/denominator) # Nash's Kurtosis Definition
    elif not i:
        x=float(0)
    else:
        x=float(0)
    kurt_val.append(x)
    mean_val.append(np.mean(i))
    variance_val.append(np.var(i))
    skew_val.append(skew(i))
    sixth_mom.append(moment(i, moment=6))


# In[58]:

from scipy.stats import kurtosis
kurt_val2 = []
kurt_val3 = []
kurt_val4 = []
kurt_val5 = []
for i in [[0] if not arr else arr for arr in ultimate]:
    kurt_val2.append(kurtosis(i, bias=False, fisher = True)) #Fisher's definition of kurtosis
    kurt_val3.append(kurtosis(i, bias=False, fisher = False)) #Pearson's definition of kurtosis
    kurt_val4.append(kurtosis(i, bias=True, fisher = True)) #Fisher's definition of kurtosis with bias
    kurt_val5.append(kurtosis(i, bias=True, fisher = False)) #Pearson's definition of kurtosis with bias
# If bias is False then the kurtosis is calculated using k statistics to eliminate bias coming from biased moment estimators


# In[59]:


# dictionary of lists 
dict1 = {'Kurtosis values (Nash)': kurt_val}  
kurt1 = pd.DataFrame(dict1) 
dict2 = {'Mean values': mean_val}  
mean1 = pd.DataFrame(dict2) 
dict3 = {'Variance values': variance_val}  
variance1 = pd.DataFrame(dict3) 
dict4 = {'Skewness values': skew_val}  
skew1 = pd.DataFrame(dict4) 
dict5 = {'PCS lengths': ultimate2}  
pcs1 = pd.DataFrame(dict5) 
dict6 = {'Kurtosis values (Fisher)': kurt_val2}  
kurt2 = pd.DataFrame(dict6) 
dict7 = {'Kurtosis values (Pearson)': kurt_val3}  
kurt3 = pd.DataFrame(dict7)
dict8 = {'Sixth moment values': sixth_mom}  
sixth = pd.DataFrame(dict8)
dict9 = {'Kurtosis values (Fisher) biased': kurt_val4}  
kurt4 = pd.DataFrame(dict9) 
dict10 = {'Kurtosis values (Pearson) biased': kurt_val5}  
kurt5 = pd.DataFrame(dict10)


# saving the dataframe to csv files
kurt1.to_csv(r'/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Kurtosis_val_Nash_hg19gorGor3_all.csv') 
kurt2.to_csv(r'/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Kurtosis_val_Fisher_hg19gorGor3_all.csv') 
kurt3.to_csv(r'/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Kurtosis_val_Pearson_hg19gorGor3_all.csv') 
kurt4.to_csv(r'/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Kurtosis_val_Fisher_biased_hg19gorGor3_all.csv') 
kurt5.to_csv(r'/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Kurtosis_val_Pearson_biased_hg19gorGor3_all.csv') 
mean1.to_csv(r'/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Mean_val_hg19gorGor3_all.csv') 
variance1.to_csv(r'/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Variance_val_hg19gorGor3_all.csv') 
skew1.to_csv(r'/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Skewness_val_hg19gorGor3_all.csv') 
pcs1.to_csv(r'/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/PCS_val_hg19gorGor3_all.csv') 
sixth.to_csv(r'/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Sixth_moment_val_hg19gorGor3_all.csv') 

# In[60]:

#get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#plt.rcParams.update({'figure.max_open_warning': 0})

x = np.linspace( -5, 5, 10000 )
y = kurt_val

plt.figure(figsize=[7, 5])

plt.hist(y, bins='auto')
plt.title('Kurtosis density (Nash\'s equation)')
plt.xlabel('Kurtosis values', fontsize=13)
plt.ylabel('Counts', fontsize=13)

plt.tight_layout()
plt.savefig('/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Nash_kurt.png')


# In[61]:


norm = pd.DataFrame({'Kurtosis values (Nash)':kurt_val})
norm.plot(kind = 'density')
plt.savefig('/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Density_plot.png')


# In[62]:


x = np.linspace( -5, 5, 10000 )
y = kurt_val2 

import matplotlib.pyplot as plt
#plt.rcParams.update({'figure.max_open_warning': 0})

plt.figure(figsize=[7, 5])

plt.hist(y, bins='auto')
plt.title('Kurtosis density (Fisher)')
plt.xlabel('Kurtosis values', fontsize=13)
plt.ylabel('Counts', fontsize=13)
plt.tight_layout()
chrom = '13'
plt.savefig('/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Fisher_kurt.png')


# In[63]:


x = np.linspace( -5, 5, 10000 )
y = kurt_val3  # normal distribution

import matplotlib.pyplot as plt
#plt.rcParams.update({'figure.max_open_warning': 0})

plt.figure(figsize=[7, 5])

plt.hist(y, bins='auto')
plt.title('Kurtosis density (Pearson)')
plt.xlabel('Kurtosis values', fontsize=13)
plt.ylabel('Counts', fontsize=13)
plt.tight_layout()
plt.savefig('/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Pearson_kurt.png')



x = np.linspace( -5, 5, 10000 )
y = kurt_val4  # normal distribution

import matplotlib.pyplot as plt
#plt.rcParams.update({'figure.max_open_warning': 0})

plt.figure(figsize=[7, 5])

plt.hist(y, bins='auto')
plt.title('Kurtosis density (Fisher) biased')
plt.xlabel('Kurtosis values (biased)', fontsize=13)
plt.ylabel('Counts', fontsize=13)
plt.tight_layout()
plt.savefig('/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Fisher_biased_kurt.png')




x = np.linspace( -5, 5, 10000 )
y = kurt_val5  # normal distribution

import matplotlib.pyplot as plt
#plt.rcParams.update({'figure.max_open_warning': 0})

plt.figure(figsize=[7, 5])

plt.hist(y, bins='auto')
plt.title('Kurtosis density (Pearson) biased')
plt.xlabel('Kurtosis values (biased)', fontsize=13)
plt.ylabel('Counts', fontsize=13)
plt.tight_layout()
plt.savefig('/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Pearson_biased_kurt.png')






# In[64]:


x = np.linspace( -5, 5, 10000 )
y = mean_val  

import matplotlib.pyplot as plt
#plt.rcParams.update({'figure.max_open_warning': 0})

plt.figure(figsize=[7, 5])

plt.hist(y, bins='auto')
plt.title('Mean density')
plt.xlabel('Mean values', fontsize=13)
plt.ylabel('Counts', fontsize=13)
plt.tight_layout()
plt.savefig('/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Mean.png')


# In[65]:


x = np.linspace( -5, 5, 10000 )
y = variance_val 

import matplotlib.pyplot as plt
#plt.rcParams.update({'figure.max_open_warning': 0})

plt.figure(figsize=[7, 5])

plt.hist(y, bins='auto')
plt.title('Variance density')
plt.xlabel('Variance values', fontsize=13)
plt.ylabel('Counts', fontsize=13)
plt.tight_layout()
plt.savefig('/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Variance.png')


# In[66]:


x = np.linspace( -5, 5, 10000 )
y = skew_val

import matplotlib.pyplot as plt
#plt.rcParams.update({'figure.max_open_warning': 0})

plt.figure(figsize=[7, 5])

plt.hist(y, bins='auto')
plt.title('Skewness density')
plt.xlabel('Skewness values', fontsize=13)
plt.ylabel('Counts', fontsize=13)
plt.tight_layout()
plt.savefig('/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_7500/Skewness.png')




