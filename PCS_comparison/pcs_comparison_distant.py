

from pandas._libs.algos import diff_2d
import os
import glob
import matplotlib.pyplot as plt
import csv
import pandas as pd
import scipy
import lmfit

###### gorilla #######
df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-calJac4/PCS_hg38calJac4.csv") 
x1 = [i for i in df1.index][1:]
y1 = df1['Number of PCS of Length L'][1:]

###### mouse ######
df2 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-danRer10/PCS_hg38danRer10.csv") 
x2 = [i for i in df2.index][1:]
y2 = df2['Number of PCS of Length L'][1:]

###### malayan lemur ######
df3 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-equCab3/PCS_hg38equCab3.csv") 
x3 = [i for i in df3.index][1:]
y3 = df3['Number of PCS of Length L'][1:]

###### mouse lemur ######
df4 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-felCat9/PCS_hg38felCat9.csv") 
x4 = [i for i in df4.index][1:]
y4 = df4['Number of PCS of Length L'][1:]

###### orangutan ######
df5 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-galGal6/PCS_hg38galGal6.csv") 
x5 = [i for i in df5.index][1:]
y5 = df5['Number of PCS of Length L'][1:]

###### rhesus ######
df6 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-ornAna2/PCS_hg38ornAna2.csv") 
x6 = [i for i in df6.index][1:]
y6 = df6['Number of PCS of Length L'][1:]

###### gibbon ######
df7 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-oryCun2/PCS_hg38oryCun2.csv") 
x7 = [i for i in df7.index][1:]
y7 = df7['Number of PCS of Length L'][1:]

###### squirrel monkey ######
df8 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-oviAri4/PCS_hg38oviAri4.csv") 
x8 = [i for i in df8.index][1:]
y8 = df8['Number of PCS of Length L'][1:]

###### squirrel monkey ######
df9 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-petMar3/PCS_hg38petMar3.csv") 
x9 = [i for i in df9.index][1:]
y9 = df9['Number of PCS of Length L'][1:]

###### squirrel monkey ######
df10 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-susScr11/PCS_hg38susScr11.csv") 
x10 = [i for i in df10.index][1:]
y10 = df10['Number of PCS of Length L'][1:]

###### squirrel monkey ######
df11 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-tarSyr2/PCS_hg38tarSyr2.csv") 
x11 = [i for i in df11.index][1:]
y11 = df11['Number of PCS of Length L'][1:]

###### squirrel monkey ######
df12 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-thaSir1/PCS_hg38thaSir1.csv") 
x12 = [i for i in df12.index][1:]
y12 = df12['Number of PCS of Length L'][1:]

###### squirrel monkey ######
df13 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-triMan1/PCS_hg38triMan1.csv") 
x13 = [i for i in df13.index][1:]
y13 = df13['Number of PCS of Length L'][1:]

###### squirrel monkey ######
df14 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-xenTro10/PCS_hg38xenTro10.csv") 
x14 = [i for i in df14.index][1:]
y14 = df14['Number of PCS of Length L'][1:]


# Resizing the figure
plt.figure(figsize=[10, 7])
plt.xlim(0.7,10**4)
plt.ylim(0.5,10**9)
# Plotting the graph with Log ticks at x and y axis using loglog
plt.loglog(x1, y1, 'blue', marker='.', linewidth=2, linestyle='None', label='hg38/calJac4')
plt.loglog(x2, y2, 'red', marker='.', linewidth=2, linestyle='None', label='hg38/danRer10')
plt.loglog(x3, y3, 'crimson', marker='.', linewidth=2, linestyle='None', label='hg38/equCab3')
plt.loglog(x4, y4, 'firebrick', marker='.', linewidth=2, linestyle='None', label='hg38/felCat9')
plt.loglog(x5, y5, 'darkorange', marker='.', linewidth=2, linestyle='None', label='hg38/galGal6')
plt.loglog(x6, y6, 'black', marker='.', linewidth=2, linestyle='None', label='hg38/ornAna2')
plt.loglog(x7, y7, 'green', marker='.', linewidth=2, linestyle='None', label='hg38/oryCun2')
plt.loglog(x8, y8, 'olive', marker='.', linewidth=2, linestyle='None', label='hg38/oviAri4')
plt.loglog(x9, y9, 'gold', marker='.', linewidth=2, linestyle='None', label='hg38/petMar3')
plt.loglog(x10, y10, 'teal', marker='.', linewidth=2, linestyle='None', label='hg38/susScr11')
plt.loglog(x11, y11, 'turquoise', marker='.', linewidth=2, linestyle='None', label='hg38/tarSyr2')
plt.loglog(x12, y12, 'indigo', marker='.', linewidth=2, linestyle='None', label='hg38/thaSir1')
plt.loglog(x13, y13, 'deeppink', marker='.', linewidth=2, linestyle='None', label='hg38/triMan1')
plt.loglog(x14, y14, 'darkslategray', marker='.', linewidth=2, linestyle='None', label='hg38/xenTro10')

plt.title('PCS Length Distribution Comparison Among Distant Genomes', fontsize=15)
plt.xlabel('Length (log10)', fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.legend(fontsize=15, markerscale=3, loc='upper right')

plt.savefig('/bucket/.mabuya/MillerU/Abrar/PCS_comparison/All_PCS_comparison_distant.png')























