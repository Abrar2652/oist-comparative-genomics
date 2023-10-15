#!pip install lmfit
from pandas._libs.algos import diff_2d
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import scipy

###### Mouse Lemur #######
df1 = pd.read_csv("/content/PCS_hg38micMur3.csv") 
x1 = [i for i in df1.index][1:]
y1 = df1['Number of PCS of Length L'][1:]

###### Zebrafish ######
df2 = pd.read_csv("/content/PCS_hg38danRer10.csv") 
x2 = [i for i in df2.index][1:]
y2 = df2['Number of PCS of Length L'][1:]

###### Horse ######
df3 = pd.read_csv("/content/PCS_hg38equCab3.csv") 
x3 = [i for i in df3.index][1:]
y3 = df3['Number of PCS of Length L'][1:]

###### Chicken ######
df5 = pd.read_csv("/content/PCS_hg38galGal6.csv") 
x5 = [i for i in df5.index][1:]
y5 = df5['Number of PCS of Length L'][1:]

###### Platypus ######
df6 = pd.read_csv("/content/PCS_hg38ornAna2.csv") 
x6 = [i for i in df6.index][1:]
y6 = df6['Number of PCS of Length L'][1:]

###### Mouse ######
df7 = pd.read_csv("/content/PCS_hg38mm39.csv") 
x7 = [i for i in df7.index][1:]
y7 = df7['Number of PCS of Length L'][1:]

###### Lamprey ######
df9 = pd.read_csv("/content/PCS_hg38petMar3.csv") 
x9 = [i for i in df9.index][1:]
y9 = df9['Number of PCS of Length L'][1:]

###### Tarsier ######
df11 = pd.read_csv("/content/PCS_hg38tarSyr2.csv") 
x11 = [i for i in df11.index][1:]
y11 = df11['Number of PCS of Length L'][1:]

###### Garter Snake ######
df12 = pd.read_csv("/content/PCS_hg38thaSir1.csv") 
x12 = [i for i in df12.index][1:]
y12 = df12['Number of PCS of Length L'][1:]

###### Manatee ######
df13 = pd.read_csv("/content/PCS_hg38triMan1.csv") 
x13 = [i for i in df13.index][1:]
y13 = df13['Number of PCS of Length L'][1:]

###### X. tropicalis ######
df14 = pd.read_csv("/content/PCS_hg38xenTro10.csv") 
x14 = [i for i in df14.index][1:]
y14 = df14['Number of PCS of Length L'][1:]

###### Malayan Mouse Lemur ######
df15 = pd.read_csv("/content/PCS_hg38galVar1.csv") 
x15 = [i for i in df15.index][1:]
y15 = df15['Number of PCS of Length L'][1:]


# Resizing the figure
plt.figure(figsize=[10, 7])
plt.xlim(0.7,10**4)
plt.ylim(0.5,10**9)
# Plotting the graph with Log ticks at x and y axis using loglog
plt.loglog(x1, y1, 'blue', marker='.', linewidth=2, linestyle='None', label='hg38/micMur3: 155')
plt.loglog(x2, y2, 'red', marker='.', linewidth=2, linestyle='None', label='hg38/danRer10: 35')
plt.loglog(x3, y3, 'crimson', marker='.', linewidth=2, linestyle='None', label='hg38/equCab3: 135')
plt.loglog(x5, y5, 'darkorange', marker='.', linewidth=2, linestyle='None', label='hg38/galGal6: 69')
plt.loglog(x6, y6, 'black', marker='.', linewidth=2, linestyle='None', label='hg38/ornAna2: 80')
plt.loglog(x7, y7, 'green', marker='.', linewidth=2, linestyle='None', label='hg38/mm39: 105')
plt.loglog(x9, y9, 'gold', marker='.', linewidth=2, linestyle='None', label='hg38/petMar3: 30')
plt.loglog(x11, y11, 'turquoise', marker='.', linewidth=2, linestyle='None', label='hg38/tarSyr2: 146')
plt.loglog(x12, y12, 'indigo', marker='.', linewidth=2, linestyle='None', label='hg38/thaSir1: 55')
plt.loglog(x13, y13, 'deeppink', marker='.', linewidth=2, linestyle='None', label='hg38/triMan1: 120')
plt.loglog(x14, y14, 'darkslategray', marker='.', linewidth=2, linestyle='None', label='hg38/xenTro10: 48')
plt.loglog(x15, y15, 'darkslategray', marker='.', linewidth=2, linestyle='None', label='hg38/galVar1: 140')


plt.hlines(y=10**2.3, lw=1,xmin=0, xmax=10**4, colors='r')
plt.axvline(x=30, lw=1, ymin=0, ymax=10**9, color='r')
plt.axvline(x=35, lw=1,ymin=0, ymax=10**9, color='r')
plt.axvline(x=48, lw=1,ymin=0, ymax=10**9, color='r')
plt.axvline(x=55, lw=1,ymin=0, ymax=10**9, color='r')
plt.axvline(x=69, lw=1,ymin=0, ymax=10**9, color='r')
plt.axvline(x=80, lw=1,ymin=0, ymax=10**9, color='r')
plt.axvline(x=120, lw=1,ymin=0, ymax=10**9, color='r')
plt.axvline(x=135, lw=1,ymin=0, ymax=10**9, color='r')
plt.axvline(x=105, lw=1,ymin=0, ymax=10**9, color='r')
plt.axvline(x=146, lw=1,ymin=0, ymax=10**9, color='r')
plt.axvline(x=155, lw=1,ymin=0, ymax=10**9, color='r')
plt.axvline(x=140, lw=1,ymin=0, ymax=10**9, color='r')


plt.title('PCS Length Distribution Comparison Among Distant Genomes', fontsize=15)
plt.xlabel('Length (log10)', fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.legend(fontsize=15, markerscale=3, loc='upper right')

plt.savefig('All_PCS_displacements_distant.png')

plt.show()
