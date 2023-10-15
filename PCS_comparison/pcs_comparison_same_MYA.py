

from pandas._libs.algos import diff_2d
import os
import glob
import matplotlib.pyplot as plt
import csv
import pandas as pd
import scipy
import lmfit

###### malayan lemur ######
df3 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-equCab3/PCS_hg38equCab3.csv") 
x3 = [i for i in df3.index][1:]
y3 = df3['Number of PCS of Length L'][1:]

###### mouse lemur ######
df4 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-felCat9/PCS_hg38felCat9.csv") 
x4 = [i for i in df4.index][1:]
y4 = df4['Number of PCS of Length L'][1:]

###### squirrel monkey ######
df8 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-oviAri4/PCS_hg38oviAri4.csv") 
x8 = [i for i in df8.index][1:]
y8 = df8['Number of PCS of Length L'][1:]

###### squirrel monkey ######
df10 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-susScr11/PCS_hg38susScr11.csv") 
x10 = [i for i in df10.index][1:]
y10 = df10['Number of PCS of Length L'][1:]


# Resizing the figure
plt.figure(figsize=[10, 7])
plt.xlim(0.7,10**3)
plt.ylim(0.5,10**8)
# Plotting the graph with Log ticks at x and y axis using loglog

plt.loglog(x3, y3, 'crimson', marker='.', linewidth=2, linestyle='None', label='hg38/equCab3')
plt.loglog(x4, y4, 'yellow', marker='.', linewidth=2, linestyle='None', label='hg38/felCat9')
plt.loglog(x8, y8, 'olive', marker='.', linewidth=2, linestyle='None', label='hg38/oviAri4')
plt.loglog(x10, y10, 'teal', marker='.', linewidth=2, linestyle='None', label='hg38/susScr11')


plt.title('PCS Length Distribution Comparison Among Similar Distant Genomes', fontsize=15)
plt.xlabel('Length (log10)', fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.legend(fontsize=15, markerscale=3, loc='upper right')

plt.savefig('/bucket/.mabuya/MillerU/Abrar/PCS_comparison/All_PCS_comparison_same_MYA.png')























