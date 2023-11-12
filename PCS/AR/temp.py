#!/usr/bin/env python
# coding: utf-8

# ## Md. Abrar Jahin
# ### Research Intern, Miller Unit, OIST

#get_ipython().system('ls')
genome = 'gorGor6'

# hg/gg IGS has unique largest 3327 sequences if we sort them in descending order
# hg/mm IGS has unique largest 2299 sequences if we sort them in descending order

import pandas as pd
import numpy as np
new = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/indel_terminated/hg38_{0}_all_vanilla_indel_terminated_seqs.tsv".format(genome),sep='\t')
new = new.sort_values('pcs.width',ascending=False)
df = new[['first.sequence','second.sequence']]
lines = []
for i,j in zip(df['first.sequence'].astype(str),df['second.sequence'].astype(str)):
    lines.append(i)
    lines.append(j)
lines_frame = pd.DataFrame(lines, columns=['Genome Sequence'])


import os
import re
from itertools import starmap
import itertools

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
        if (letter1.islower() and letter2.islower() and letter1==letter2):
            result1+='\n'
            result2+='\n'
        if ((letter1 not in ['A','T','C','G']) or (letter2 not in ['A','T','C','G'])) and letter1==letter2:
            result1+='\n'
            result2+='\n'
        if ((letter1 == letter2) and letter1.isupper() and letter2.isupper() and (letter1!='-') and (letter2!='-') and 
            letter1 in ['A','T','C','G'] and letter2 in ['A','T','C','G']):
            result1+=letter1
            result2+=letter2
        if ((letter1 != letter2) or (letter1=='-') or (letter2=='-')):
            result1+='\n'
            result2+='\n'

    word=re.sub(r'\n+', '\n', result1).strip()
    return word
    

# seq_count stores the PCS lengthss
len_of_segment = new['pcs.width'].values.tolist()
ultimate = []

for i in range(int(len(lines_frame))-1):  #Iterates through the whole genome
    seq_count = []
    if i%2==0: # only even rows are input, but inside the function both even and odd lines are compared
        for x in (compare_strings(lines_frame, "Genome Sequence", i).split()):
            seq_count.append(len(x)) # append pcs lengths 
        ultimate.append(seq_count) # handle null pcs, append pcs lengths for each sequence
    
cum_pcs_per_igs = list(itertools.accumulate(ultimate)) # cumulatively concatenates pcs of sequence to the pcs of the next sequence
# the sequences is sorted in descending order
# only "cumulative" PCS counts are reported, starting from the largest IGS 
# then largest IGS + second largest IGS
# then largest IGS + second largest IGS + third largest IGS and so on
cum_pcs_frame = pd.DataFrame(list(zip(len_of_segment, cum_pcs_per_igs)), columns=['length_of_each_sequence','cumulative_pcs_descending_order'])
cum_pcs_frame.to_csv(r'/bucket/MillerU/Abrar/PCS/AR/cumulative_vanilla_pcs_for_each_IGS_whole_genome/cumulative_vanilla_pcs_in_vanilla_IGS_hg38{0}_all.csv'.format(genome), index = False) 

# RUN THE FOLLOWING CODE ONLY FOR THE ALL_MERGED FILE
dct1 = {}
dct2 = {}
for i,k in enumerate(cum_pcs_per_igs):
  #print(k)
  serial_count=[]
  for j in range(0, int(max(k)+5)):
      serial_count.append(k.count(j))

  dct1['s{0}'.format(i)] = [x for x in range(int(max(k)+5))]
  dct2['t{0}'.format(i)] = serial_count

### Function for generating unique colors for distinct plots ######
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib._cm
#https://matplotlib.org/stable/gallery/color/colormap_reference.html
#https://stackoverflow.com/questions/4971269/how-to-pick-a-new-color-for-each-plotted-line-within-a-figure-in-matplotlib
def get_cycle(cmap, N=None, use_index="auto"):
    if isinstance(cmap, str):
        if use_index == "auto":
            if cmap in ['Pastel1', 'Pastel2', 'Paired', 'Accent',
                        'Dark2', 'Set1', 'Set2', 'Set3',
                        'tab10', 'tab20', 'tab20b', 'tab20c']:
                use_index=True
            else:
                use_index=False
        cmap = plt.get_cmap(cmap)
    if not N:
        N = cmap.N
    if use_index=="auto":
        if cmap.N > 100:
            use_index=False
        elif isinstance(cmap, LinearSegmentedColormap):
            use_index=False
        elif isinstance(cmap, ListedColormap):
            use_index=True
    if use_index:
        ind = np.arange(int(N)) % cmap.N
        return plt.cycler("color",cmap(ind))
    else:
        colors = cmap(np.linspace(0,1,N))
        return plt.cycler("color",colors)

from matplotlib import pyplot as plt
from matplotlib.colors import hsv_to_rgb
from cycler import cycler
from matplotlib.cm import get_cmap
# Resizing the figure 
plt.figure(figsize=[10, 7])

with plt.rc_context({"axes.prop_cycle": get_cycle('gist_rainbow',len(cum_pcs_per_igs))}):
  for i,k in enumerate(cum_pcs_per_igs):
    plt.loglog(dct1['s{0}'.format(i)], dct2['t{0}'.format(i)], marker='.', alpha=0.5, linewidth=2, linestyle='None', label = "{0}".format(i+1))

# This IGS denotes the whole genome indel terminated sequences (vanilla+lowercase)
plt.title('Cumulative Vanilla PCS Length Distribution Comparison\n For Each Vanilla IGS in hg38/{0}'.format(genome), fontsize=14)

plt.xlabel('PCS Length (log10)', fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
#plt.legend(title="IGS in descending order: ".format(i), title_fontsize='large', ncol=5, fontsize=9, markerscale=1.5)
plt.savefig('/bucket/MillerU/Abrar/PCS/AR/cumulative_vanilla_pcs_for_each_IGS_whole_genome/cumulative_vanilla_pcs_in_vanilla_IGS_hg38{0}_all.png'.format(genome))























