

from pandas._libs.algos import diff_2d
import os
import glob
import matplotlib.pyplot as plt
import csv
import pandas as pd
import scipy
import lmfit

###### gorilla #######
df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-gorGor6/PCS_hg38gorGor6.csv") 
x1 = [i for i in df1.index][1:]
y1 = df1['Number of PCS of Length L'][1:]

###### mouse ######
df2 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-mm39/PCS_hg38mm39.csv") 
x2 = [i for i in df2.index][1:]
y2 = df2['Number of PCS of Length L'][1:]

###### malayan lemur ######
df3 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-galVar1/PCS_hg38galVar1.csv") 
x3 = [i for i in df3.index][1:]
y3 = df3['Number of PCS of Length L'][1:]

###### mouse lemur ######
df4 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-micMur3/PCS_hg38micMur3.csv") 
x4 = [i for i in df4.index][1:]
y4 = df4['Number of PCS of Length L'][1:]

###### orangutan ######
df5 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-ponAbe3/PCS_hg38ponAbe3.csv") 
x5 = [i for i in df5.index][1:]
y5 = df5['Number of PCS of Length L'][1:]

###### rhesus ######
df6 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-rheMac10/PCS_hg38rheMac10.csv") 
x6 = [i for i in df6.index][1:]
y6 = df6['Number of PCS of Length L'][1:]

###### gibbon ######
df7 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-nomLeu3/PCS_hg38nomLeu3.csv") 
x7 = [i for i in df7.index][1:]
y7 = df7['Number of PCS of Length L'][1:]

###### squirrel monkey ######
df8 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-saiBol1/PCS_hg38saiBol1.csv") 
x8 = [i for i in df8.index][1:]
y8 = df8['Number of PCS of Length L'][1:]

# Resizing the figure
plt.figure(figsize=[10, 7])

# Plotting the graph with Log ticks at x and y axis using loglog
plt.loglog(x1, y1, 'blue', marker='.', linewidth=2, linestyle='None', label='Log-Log (hg38/gorGor6)')
plt.loglog(x2, y2, 'red', marker='.', linewidth=2, linestyle='None', label='Log-Log (hg38/mm39)')
plt.loglog(x3, y3, 'gold', marker='.', linewidth=2, linestyle='None', label='Log-Log (hg38/galVar1)')
plt.loglog(x4, y4, 'green', marker='.', linewidth=2, linestyle='None', label='Log-Log (hg38/micMur3)')
plt.loglog(x5, y5, 'darkorange', marker='.', linewidth=2, linestyle='None', label='Log-Log (hg38/ponAbe3)')
plt.loglog(x6, y6, 'black', marker='.', linewidth=2, linestyle='None', label='Log-Log (hg38/rheMac10)')
plt.loglog(x7, y7, 'firebrick', marker='.', linewidth=2, linestyle='None', label='Log-Log (hg38/nomLeu3)')
plt.loglog(x8, y8, 'olive', marker='.', linewidth=2, linestyle='None', label='Log-Log (hg38/saiBol1)')

plt.title('PCS Length Distribution Comparison', fontsize=15)
plt.xlabel('Length (log10)', fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.legend(fontsize=15, markerscale=3, loc='lower left')

plt.savefig('/bucket/.mabuya/MillerU/Abrar/PCS_comparison/All_PCS_comparison_primate_distant.png')























